import hashlib

try:
    from openbabel import pybel  # openbabel 3.0.0

    GetSymbol = pybel.ob.GetSymbol
    GetVdwRad = pybel.ob.GetVdwRad
except ImportError:
    import pybel  # openbabel 2.4

    table = pybel.ob.OBElementTable()
    GetSymbol = table.GetSymbol
    GetVdwRad = table.GetVdwRad

from autoqchem.gaussian_input_generator import *
from autoqchem.openbabel_functions import *

logger = logging.getLogger(__name__)


class molecule(object):
    """Wrapper class for openbabel.OBMol class"""

    def __init__(self,
                 input,
                 isotopes_as_labels=False,
                 max_num_conformers=1,
                 input_type="string",
                 input_format='smi',
                 gen3D_option='best',
                 min_fragment_dist=2,
                 ):
        """
        Initialize the OBMol molecule, generate initial 3D geometry and do conformer search.

        :param input: string or file path
        :param isotopes_as_labels: a flag to treat isotopes in the smiles string as atomic labels
        :param max_num_conformers: maximum number of conformers to generate, use 1 (default) for no conformer search, \
        a reasonable number of conformations to use is 30. Conformations are generated with genetic algorithm \
        using pybel.ob.OBConformerSearch class.
        :param input_type: "string" or "file", in line with the input
        :param input_format: any format supported by OpenBabel, e.g. 'smi', 'cdx', 'pdb', etc.
        :param gen3D_option: "best", "medium", "fast" or "gen2D" (no 3D generation)
        :param min_fragment_dist: minimum distance between molecular fragments for salts, no more than 2 \
        molecular fragments are supported
        """

        # read the molecule
        self.mol = input_to_OBMol(input, input_type, input_format)

        # get canonical smiles of this molecule
        self.can = OBMol_to_string(self.mol, "can")
        logger.info(f"Initializing molecule with canonical smiles: {self.can}")

        # reload molecule from can
        self.mol = input_to_OBMol(self.can, "string", "can")

        # create a unique name for files and directories, aka filesystem name (use stoichiometric formula)
        # add 4 hash digits of its canonical smiles in case of collisions of formulas
        self.fs_name = f"{self.mol.GetFormula()}_{hashlib.md5(self.can.encode()).hexdigest()[:4]}"

        # TODO a chemical name would be good to have for lookups (smiles strings are long and hard to handle)

        # generate initial geometry and conformations
        self._generate_geometry(gen3D_option)
        self._generate_conformers(max_num_conformers)
        self.max_num_conformers = max_num_conformers

        # find central atoms
        self._find_central_atoms()

        # save isotopes as labels flag
        self.isotopes_as_labels = isotopes_as_labels

        # extra steps for molecules with multiple fragments
        if len(self.centers) == 2:
            logger.info(f"Molecule has 2 non-bonded fragments")
            # adjust distance between fragments (in-case it's not enough already)
            self._adjust_geometries(min_fragment_dist)
        elif len(self.centers) > 2:
            message = f"Molecule has {len(self.centers)} non-bonded fragments. Only up to 2 are supported"
            logger.error(message)
            raise Exception(message)

    def get_geometry(self, conformer_num=0) -> pd.DataFrame:
        """Get coordinates DataFrame for a given conformer.

        :param conformer_num: conformer number, if larger than number of conformations available, the last is returned
        :return: pandas.core.frame.DataFrame
        """

        self.mol.SetConformer(conformer_num)
        array = [[GetSymbol(a.GetAtomicNum()),
                  a.GetIsotope(), a.x(), a.y(), a.z()] for a in pybel.ob.OBMolAtomIter(self.mol)]
        return pd.DataFrame(array, columns=['Atom', 'Isotope'] + list('XYZ'))

    def draw(self, conformer_num=0, ipython_3d=True) -> None:
        """Draw a depiction of the molecule for a given conformer.

        :param conformer_num: conformer number, if larger than number of conformations available, the last is returned
        :param ipython_3d: if True 3-dim rotatable graphics is generated
        """

        self.mol.SetConformer(conformer_num)
        pybel.ipython_3d = ipython_3d

        try:
            from IPython.core.display import display
            display(pybel.Molecule(self.mol))
        except ImportError as e:
            logger.warning(f"Import Error: {e}")

    def _generate_geometry(self, gen3D_option) -> None:
        """Generate initial geometry of the molecule.

        :param gen3D_option: "best", "medium", "fast" or "gen2D" (no 3D generation)
        """

        logger.info(f"Creating initial geometry with option '{gen3D_option}'.")
        if gen3D_option == "gen2D":
            self.mol.AddHydrogens()
            gen2D = pybel.ob.OBOp.FindType("gen2D")
            gen2D.Do(self.mol)
        else:
            gen3D = pybel.ob.OBOp.FindType("gen3D")
            gen3D.Do(self.mol, gen3D_option)
        logger.info(f"Initial geometry created successfully.")

    def _find_central_atoms(self) -> None:
        """Find central atoms for molecular fragments of the molecule."""

        fragments = [[atom.GetId() for atom in pybel.ob.OBMolAtomIter(part)] for part in self.mol.Separate()]
        self.fragments_dict = {i: j for j, frag in enumerate(fragments) for i in frag}

        # find centers for each fragment
        geom = self.get_geometry()
        geom['Fragment'] = geom.index.map(self.fragments_dict)
        self.centers = geom.groupby('Fragment')[['X', 'Y', 'Z']].apply(
            lambda frag: frag.sub(frag.mean()).pow(2).sum(1).idxmin()
        ).values

    def _generate_conformers(self, num_conformers) -> None:
        """Generate conformations with genetic algorithm using pybel.ob.OBConformerSearch class.

        :param num_conformers: maximum number of conformers to generate
        """

        # safety check
        assert num_conformers > 0

        # skip trivial case
        if num_conformers < 2:
            return

        confSearch = pybel.ob.OBConformerSearch()
        confSearch.Setup(self.mol, num_conformers)
        confSearch.Search()
        confSearch.GetConformers(self.mol)

        logger.info(f"Conformer Search generated {self.mol.NumConformers()} conformations of {self.can} molecule")

    def get_light_and_heavy_elements(self, max_light_atomic_number) -> tuple:
        """Group molecule elements into light and heavy.

        :param max_light_atomic_number: maximum atomic number classified as light
        :return: light, heavy element lists
        """

        atomic_nums = set(atom.GetAtomicNum() for atom in pybel.ob.OBMolAtomIter(self.mol))
        light_elements = [GetSymbol(n) for n in atomic_nums if n <= max_light_atomic_number]
        heavy_elements = [GetSymbol(n) for n in atomic_nums if n > max_light_atomic_number]
        return light_elements, heavy_elements

    def _adjust_geometries(self, min_fragment_dist) -> None:
        """Adjust molecule fragment geometries such that the minimum separation
        between any atom between fragments is at least 'min_fragment_dist' in Angstroms.

        :param min_fragment_dist: minimum distance between molecular fragments for salts, no more than 2 \
        molecular fragments are supported
        """

        for conf_id in range(self.mol.NumConformers()):
            geom = self.get_geometry(conf_id)

            geom['Fragment'] = geom.index.map(self.fragments_dict)
            geom = geom.set_index(['Fragment', 'Atom'])
            v = geom.iloc[self.centers].diff().dropna().values  # separate along the center-center axis
            v = v / np.linalg.norm(v)

            init_mdist = cdist(geom.loc[0], geom.loc[1]).min()
            counter = 0
            mdist = init_mdist
            while mdist < min_fragment_dist:
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
