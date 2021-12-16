from autoqchem.helper_classes import *
from autoqchem.molecule import pybel, GetSymbol

conv = pybel.ob.OBConversion()


def generate_conformations_from_openbabel(smiles, num_conf, ob_gen3D_option='best') -> tuple:
    """

    :param smiles: SMILES string
    :type smiles: str
    :param num_conf: maximum number of conformations to generate
    :type num_conf: int
    :param ob_gen3D_option: option to use with openbabel gen3D for search of initial geometry
    :type ob_gen3D_option: str
    :return: tuple(elements, conformer_coordinates, connectivity_matrix, charges)
    """

    # initialize obmol
    obmol = pybel.readstring('smi', smiles).OBMol
    obmol.AddHydrogens()

    # initial geometry
    gen3D = pybel.ob.OBOp.FindType("gen3D")
    gen3D.Do(obmol, ob_gen3D_option)

    # conf search
    confSearch = pybel.ob.OBConformerSearch()
    confSearch.Setup(obmol, num_conf)
    confSearch.Search()
    confSearch.GetConformers(obmol)

    elements, conformer_coordinates, connectivity_matrix, charges = extract_from_obmol(obmol)

    return elements, conformer_coordinates, connectivity_matrix, charges


def extract_from_obmol(mol) -> tuple:
    """Extract information from Openbabel OBMol object with conformers.

    :param mol: pybel.ob.OBMol object
    :type mol: pybel.ob.OBMol
    :return: tuple(elements, conformer_coordinates, connectivity_matrix, charges)
    """

    py_mol = pybel.Molecule(mol)
    elements = [GetSymbol(atom.atomicnum) for atom in py_mol.atoms]
    charges = np.array([atom.formalcharge for atom in py_mol.atoms])

    n_atoms = len(py_mol.atoms)
    connectivity_matrix = np.zeros((n_atoms, n_atoms))
    for bond in pybel.ob.OBMolBondIter(mol):
        i = bond.GetBeginAtomIdx() - 1
        j = bond.GetEndAtomIdx() - 1
        bo = bond.GetBondOrder()
        connectivity_matrix[i, j] = bo
        connectivity_matrix[j, i] = bo

    # Retrieve conformer coordinates
    conformer_coordinates = []
    for i in range(mol.NumConformers()):
        mol.SetConformer(i)
        coordinates = np.array([atom.coords for atom in py_mol.atoms])
        conformer_coordinates.append(coordinates)
    conformer_coordinates = np.array(conformer_coordinates)

    return elements, conformer_coordinates, connectivity_matrix, charges
