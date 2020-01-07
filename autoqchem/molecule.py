import hashlib

import numpy as np
import rdkit.Chem
from scipy.spatial.distance import cdist

from autoqchem.gaussian_input_generator import *
from autoqchem.openbabel_conversions import *

logger = logging.getLogger(__name__)


class molecule(object):
    """Class that holds a single molecule information"""

    def __init__(self, input, input_format='smi', input_type=input_types.string, ):
        """initialize the molecule"""

        # read the molecule
        self.mol = input_to_OBMol(input, input_format, input_type)

        # get canonical smiles of this molecule
        self.name = OBMol_to_string(self.mol, "can")
        logger.info(f"Initializing molecule with canonical smiles: {self.name}")

        # reload the molecule using its canonical smiles
        self.mol = input_to_OBMol(self.name, "can", input_types.string)

        # group elements into light and heavy
        self.__get_light_and_heavy_elements()

        # create a unique name for files and directories, aka filesystem name (use stoichiometric formula)
        # add 4 hash digits of its canonical smiles in case of collisions of formulas
        self.fs_name = f"{self.mol.GetFormula()}_{hashlib.md5(self.name.encode()).hexdigest()[:4]}"

        # generate initial geometry and conformations
        self.__generate_geometry(config['babel']['gen3D_option'])
        self.__generate_conformers(config['babel']['num_conformers'])

        # find central atoms
        self.__find_central_atoms()

        # extra steps for molecules with multiple fragments
        if len(self.centers) == 2:
            logger.info(f"Molecule has 2 non-bonded fragments")
            # adjust distance between fragments (in-case it's not enough already)
            self.__adjust_geometries(config['babel']['min_dist_for_salts'])
        elif len(self.centers) > 2:
            message = f"Molecule has {len(self.centers)} non-bonded fragments. Only up to 2 are supported"
            logger.error(message)
            raise Exception(message)

    def __generate_geometry(self, option):
        """generate initial geometry using get3D option"""

        logger.info(f"Creating initial geometry")
        gen3D = pybel.ob.OBOp.FindType("gen3D")
        gen3D.Do(self.mol, option)
        logger.info(f"Initial geometry created successfully.")

    def __find_central_atoms(self):
        """First find molecular fragments (if any) and group atoms
        by the fragment they belong to, then, find central atoms
        for all fragments (closest to centroids)"""

        # use RDKit to get molecular fragments
        rdkit_mol = rdkit.Chem.MolFromSmiles(self.name)  # get rdkit mol from canonical
        rdkit_mol = rdkit.Chem.AddHs(rdkit_mol)  # add hydrogens
        fragments = rdkit.Chem.GetMolFrags(rdkit_mol)
        self.fragments_dict = {i: j for j, frag in enumerate(fragments) for i in frag}

        # find centers for each fragment
        geom = self.get_initial_geometry()
        geom['Fragment'] = geom.index.map(self.fragments_dict)
        self.centers = geom.groupby('Fragment')[['X', 'Y', 'Z']].apply(
            lambda frag: frag.sub(frag.mean()).pow(2).sum(1).idxmin()
        ).values

    def __generate_conformers(self, num_conformers):
        """generate conformations with GA algorithm"""

        confSearch = pybel.ob.OBConformerSearch()
        confSearch.Setup(self.mol, num_conformers)
        confSearch.Search()
        confSearch.GetConformers(self.mol)

        logger.info(f"Conformer Search generated {self.mol.NumConformers()} conformations of {self.name} molecule")

    def __get_light_and_heavy_elements(self):
        """group elements into light and heavy for this molecule"""

        max_light_z = config['gaussian']['max_light_atomic_number']
        atomic_nums = set(atom.GetAtomicNum() for atom in pybel.ob.OBMolAtomIter(self.mol))
        self.light_elements = [pybel.ob.GetSymbol(n) for n in atomic_nums if n <= max_light_z]
        self.heavy_elements = [pybel.ob.GetSymbol(n) for n in atomic_nums if n > max_light_z]

    def __adjust_geometries(self, min_dist):
        """adjust geometries such that the minimum separation
        between any atoms for fragments is at least 'min_dist' Angstroms"""

        for conf_id in range(self.mol.NumConformers()):
            self.mol.SetConformer(conf_id)
            geom = OBMol_to_geom_df(self.mol)

            geom['Fragment'] = geom.index.map(self.fragments_dict)
            geom = geom.set_index(['Fragment', 'Atom'])
            v = geom.iloc[self.centers].diff().dropna().values  # separate along the center-center axis
            v = v / np.linalg.norm(v)

            init_mdist = cdist(geom.loc[0], geom.loc[1]).min()
            counter = 0
            mdist = init_mdist
            while mdist < min_dist:
                geom.loc[1] = geom.loc[1].values + 0.1 * v
                mdist = cdist(geom.loc[0], geom.loc[1]).min()
                counter += 1

            if counter > 0:
                logger.info(f"Conformation {conf_id}: repositioned molecular fragments from"
                            f" {init_mdist:.2f} separation to "
                            f"{mdist:.2f} separation in {counter} iterations.")

            # reposition the atoms in the OBMol object
            for atom in pybel.ob.OBMolAtomIter(self.mol):
                pos = geom.iloc[atom.GetIdx() - 1]
                atom.SetVector(pos.X, pos.Y, pos.Z)

    def get_initial_geometry(self, conformer_num=0) -> pd.DataFrame:
        """get coordinates dataframe for a given conformer"""

        self.mol.SetConformer(conformer_num)
        return OBMol_to_geom_df(self.mol)

    def draw(self, conformer_num=0) -> pybel.Molecule:
        """draw molecule"""

        self.mol.SetConformer(conformer_num)
        return pybel.Molecule(self.mol)
