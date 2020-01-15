from autoqchem.gaussian_log_extractor import *
from autoqchem.helper_classes import *

conv = pybel.ob.OBConversion()


def input_to_OBMol(input, input_format, input_type) -> pybel.ob.OBMol:
    """create OBMol from input (type and format as supported in the enumerators)"""

    mol = pybel.ob.OBMol()
    conv.SetInFormat(input_format)

    if input_type == input_types.file:
        conv.ReadFile(mol, input)
    elif input_type == input_types.string:
        conv.ReadString(mol, input)

    return mol


def OBMol_to_string(mol, format) -> str:
    """from OBMol to a string format"""

    conv.SetOutFormat(format)
    return conv.WriteString(mol).strip()


def OBMol_to_file(mol, format, target_path) -> str:
    """from OBMol to a string format"""

    conv.SetOutFormat(format)
    return conv.WriteFile(mol, target_path)


def OBMol_to_geom_df(mol) -> pd.DataFrame:
    """extract geometry dataframe from OBMol"""

    array = [[pybel.ob.GetSymbol(a.GetAtomicNum()),
              a.GetIsotope(), a.x(), a.y(), a.z()] for a in pybel.ob.OBMolAtomIter(mol)]
    df = pd.DataFrame(array, columns=['Atom', 'Isotope'] + list('XYZ'))
    df['Atom'] = df['Atom'] + df['Isotope'].astype(str).where(df['Isotope'] > 0, '')
    return df.drop('Isotope', axis=1)


def OBMol_from_done_slurm_job(job) -> pybel.ob.OBMol:
    """create OBMol from a finished gaussian job"""

    assert job.status == slurm_status.done
    le = gaussian_log_extractor(f"{job.directory}/{job.base_name}.log")
    return OBMol_from_can_and_geom(job.can, le.geom)


def OBMol_from_can_and_geom(can, geom) -> pybel.ob.OBMol:
    """create OBMol from canonical smiles and geometry dataframe"""

    # create OBMol from can
    mol = input_to_OBMol(can, input_type=input_types.string, input_format="can")
    mol.AddHydrogens()

    # adjust geometry
    for atom in pybel.ob.OBMolAtomIter(mol):
        pos = geom.iloc[atom.GetIdx() - 1]
        atom.SetVector(pos.X, pos.Y, pos.Z)

    return mol


def deduplicate_list_of_OBMols(mols, symmetry) -> list:
    """Filter conformers based on RMSD,
    the function is slow if symmetry=True in OBAlign setup"""

    # safety check, assert all mols convert to the same canonical smiles
    assert (len(set(OBMol_to_string(mol, "can") for mol in mols)) == 1)

    # trivial case
    if len(mols) < 2:
        return []

    alignment = pybel.ob.OBAlign(True, symmetry)  # alignment class from OB
    RMSD_threshold = config['gaussian']['conformer_RMSD_threshold']

    dupliacte_indices = []
    for i in range(len(mols) - 1):
        alignment.SetRefMol(mols[i])
        for j in range(i + 1, len(mols)):
            alignment.SetTargetMol(mols[j])
            alignment.Align()
            if alignment.GetRMSD() < RMSD_threshold and i not in dupliacte_indices:
                dupliacte_indices.append(j)

    return dupliacte_indices
