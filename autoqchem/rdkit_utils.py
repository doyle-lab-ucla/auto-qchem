import os
import numpy as np

from ipywidgets import interact, fixed

import py3Dmol
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem

from .molecule import GetSymbol
from .gaussian_log_extractor import gaussian_log_extractor


def extract_from_rdmol(mol: Chem.Mol) -> tuple([list, np.ndarray, np.ndarray, np.ndarray]):
    """Extract information from RDKit Mol object with conformers."""

    # Take out elements, coordinates and connectivity matrix
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    charges = np.array([atom.GetFormalCharge() for atom in mol.GetAtoms()])

    conformer_coordinates = []
    for conformer in mol.GetConformers():
        coordinates = conformer.GetPositions()
        conformer_coordinates.append(coordinates)
    conformer_coordinates = np.array(conformer_coordinates)

    connectivity_matrix = Chem.GetAdjacencyMatrix(mol, useBO=True)

    return elements, conformer_coordinates, connectivity_matrix, charges


def get_rdkit_mol(elements: list, conformer_coordinates: np.ndarray,
                  connectivity_matrix: np.ndarray, charges: np.ndarray) -> Chem.Mol:
    _RDKIT_BOND_TYPES = {
        1.0: Chem.BondType.SINGLE,
        1.5: Chem.BondType.AROMATIC,
        2.0: Chem.BondType.DOUBLE,
        3.0: Chem.BondType.TRIPLE,
        4.0: Chem.BondType.QUADRUPLE,
        5.0: Chem.BondType.QUINTUPLE,
        6.0: Chem.BondType.HEXTUPLE,
    }

    mol = Chem.RWMol()

    # Add atoms
    for element, charge in zip(elements, charges):
        atom = Chem.Atom(element)
        atom.SetFormalCharge(int(charge))
        mol.AddAtom(atom)

    # Add bonds
    for i, j in zip(*np.tril_indices_from(connectivity_matrix)):
        if i != j:
            bo = connectivity_matrix[i, j]
            if bo != 0:
                bond_type = _RDKIT_BOND_TYPES[float(bo)]
                mol.AddBond(int(i), int(j), bond_type)

    # Add conformers
    add_conformers_to_rdmol(mol, conformer_coordinates)

    mol = mol.GetMol()
    Chem.SanitizeMol(mol)

    return mol


def add_conformers_to_rdmol(mol: Chem.Mol, conformer_coordinates: np.ndarray) -> None:
    """Add conformers to RDKit Mol object.
    Args:
        mol: RDKit mol object
        conformer_coordinates: Conformer coordinates (Å)
    """
    conformer_coordinates = np.array(conformer_coordinates)
    if len(conformer_coordinates.shape) == 2:
        conformer_coordinates.reshape(-1, conformer_coordinates.shape[0], 3)

    for coordinates in conformer_coordinates:
        conformer = Chem.Conformer()
        for i, coord in enumerate(coordinates):
            point = Geometry.Point3D(*coord)
            conformer.SetAtomPosition(i, point)
        mol.AddConformer(conformer, assignId=True)


def graph_conf(m, confId=0):
    mb = Chem.MolToMolBlock(m, confId=confId)
    p = py3Dmol.view(width=400, height=400)
    p.removeAllModels()
    p.addModel(mb, 'sdf')
    p.setStyle({'stick': {}})
    p.setBackgroundColor('0xeeeeee')
    p.zoomTo()
    return p


def draw(mol: Chem.Mol, energies: list) -> interact:
    """Make a drawing of all conformers in 3d"""

    # emin = min(energies)
    # energies = energies - emin
    # p = py3Dmol.view(width=400,height=400)
    # return interact(drawit, m=fixed(mol), p=fixed(p),
    #                confId=(0, mol.GetNumConformers()-1))


def generate_conformations_from_rdkit(smiles, num_conf, rdkit_ff='MMFF94', n_threads=os.cpu_count() - 1):
    """Generate conformational ensemble using a method of choice"""

    # initialize rdmol
    rdmol = Chem.AddHs(Chem.MolFromSmiles(smiles))

    params = AllChem.EmbedParameters()
    params.useSymmetryForPruning = True
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    params.ETversion = 2
    params.pruneRmsThresh = 0.35
    params.numThreads = n_threads

    # embed and optimized conformers
    AllChem.EmbedMultipleConfs(rdmol, num_conf, params)
    if rdkit_ff == "MMFF94":
        AllChem.MMFFOptimizeMoleculeConfs(rdmol, mmffVariant="MMFF94", numThreads=n_threads)
    elif rdkit_ff == "MMFF94s":
        AllChem.MMFFOptimizeMoleculeConfs(rdmol, mmffVariant="MMFF94s", numThreads=n_threads)
    elif rdkit_ff == "UFF":
        AllChem.UFFOptimizeMoleculeConfs(rdmol, numThreads=n_threads)

    elements, conformer_coordinates, connectivity_matrix, charges = extract_from_rdmol(rdmol)

    return elements, conformer_coordinates, connectivity_matrix, charges


def get_light_and_heavy_elements(mol: Chem.Mol, max_light_atomic_number: int) -> tuple:
    """Group molecule elements into light and heavy."""

    atomic_nums = set(atom.GetAtomicNum() for atom in mol.GetAtoms())
    light_elements = [GetSymbol(n) for n in atomic_nums if n <= max_light_atomic_number]
    heavy_elements = [GetSymbol(n) for n in atomic_nums if n > max_light_atomic_number]
    return light_elements, heavy_elements


def rdmol_from_slurm_jobs(jobs, postDFT=True, ordered=False):
    """Create an rdkit molecule from a set of finished slurm jobs"""

    # check that these are jobs for the same molecule
    assert len(set(j.inchi for j in jobs)) == 1
    if postDFT == True:
        # check that the jobs are in done status
        assert all(j.status.value == slurm_status.done.value for j in jobs)

    elements, connectivity_matrix, charges = jobs[0].elements, jobs[0].connectivity_matrix, jobs[0].charges
    conformer_coordinates = []
    energies = []
    for j in jobs:
        if postDFT == True:

            le = gaussian_log_extractor(f"{j.directory}/{j.base_name}.log")
            le.check_for_exceptions()
            le.get_atom_labels()

            # verify that the labels are in the same order in gaussian after running it
            assert tuple(le.labels) == tuple(elements)

            le.get_geometry()
            conformer_coordinates.append(le.geom[list('XYZ')].values)

            le._get_freq_part_descriptors()
            energies.append(le.descriptors['G'])

        else:
            with open(f"{j.directory}/{j.base_name}.gjf") as f:
                text = f.read()
                geom = re.findall(f"{j.base_name}\n\n\d\s\d\n(.*?)\n\n\n--Link1--", text, re.DOTALL)[0]
                geom = map(str.split, map(str.strip, geom.splitlines()))
                geom = np.array(list(geom))
                geom = geom[:, 1:].astype(float)
                conformer_coordinates.append(geom)

    rdmol = get_rdkit_mol(elements, conformer_coordinates, connectivity_matrix, charges)
    if not postDFT:
        props = AllChem.MMFFGetMoleculeProperties(rdmol)
        energies = [AllChem.MMFFGetMoleculeForceField(rdmol, props, confId=i).CalcEnergy()
                    for i in range(rdmol.GetNumConformers())]

    if ordered:
        pass

    return rdmol, energies


def get_rmsd_rdkit(rdmol):
    """Calculate RMSD row-wise with RDKit."""

    # Construct atom list for rmsd: do not include hydrogen
    atom_ids = [atom.GetIdx() for atom in rdmol.GetAtoms() if atom.GetAtomicNum() != 1]

    # Calculated RMSD row-wise with RDKit
    rmsds = []
    conformers = list(rdmol.GetConformers())
    for i in range(len(conformers)):
        ref_mol = Chem.Mol(rdmol)
        ref_conformer = conformers[i]
        ref_mol.RemoveAllConformers()
        ref_mol.AddConformer(ref_conformer)
        for j in range(len(conformers)):
            conformer = conformers[j]
            ref_mol.AddConformer(conformer)
        rmsds_row = []
        AllChem.AlignMolConformers(ref_mol, atomIds=atom_ids, RMSlist=rmsds_row)
        rmsds.append(rmsds_row)
    rmsds = np.array(rmsds)

    return rmsds


def prune_rmsds(rdmol, thres):
    """Get a list of conformer indices to keep"""

    rmsds = get_rmsd_rdkit(rdmol)

    working_array = rmsds
    candidates = np.array(range(rdmol.GetNumConformers()))

    keep_list = []
    while len(working_array) > 0:
        keeper = candidates[0]
        keep_list.append(keeper)
        rmsd = working_array[0]
        mask = rmsd > thres
        candidates = candidates[mask]
        working_array = working_array[mask, :][:, mask]

    return keep_list
