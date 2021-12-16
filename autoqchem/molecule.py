import logging

try:
    from openbabel import pybel  # openbabel 3.0.0

    GetSymbol = pybel.ob.GetSymbol
    GetVdwRad = pybel.ob.GetVdwRad
except ImportError:
    import pybel  # openbabel 2.4

    table = pybel.ob.OBElementTable()
    GetSymbol = table.GetSymbol
    GetVdwRad = table.GetVdwRad

from rdkit.Chem import Descriptors

from autoqchem.openbabel_utils import *
from autoqchem.rdkit_utils import *

logger = logging.getLogger(__name__)


class molecule(object):
    """Wrapper class for molecule"""

    def __init__(self, smiles, num_conf, engine='rdkit', rdkit_ff='MMFF94', ob_gen3D_option='best') -> None:
        """Initialize the molecule with a conformational ensemble

        :arg smiles: SMILES string
        :type smiles: str
        :arg num_conf: maximum number of conformations to generate
        :type num_conf: int
        :param engine: conformation search engine (rdkit or openbabel)
        :type engine: str
        :param rdkit_ff: rdkit supported force-field to use when engine='rdkit'
        :type rdkit_ff: str
        :param ob_gen3D_option: option to use with openbabel gen3D for search of initial geometry
        :type ob_gen3D_option: str
        """

        # run conformer generation
        if engine == 'rdkit':
            self.elements, \
            self.conformer_coordinates, \
            self.connectivity_matrix, \
            self.charges = generate_conformations_from_rdkit(smiles=smiles, num_conf=num_conf,
                                                             rdkit_ff=rdkit_ff)
        elif engine == 'openbabel':
            self.elements, \
            self.conformer_coordinates, \
            self.connectivity_matrix, \
            self.charges = generate_conformations_from_openbabel(smiles=smiles, num_conf=num_conf,
                                                                 ob_gen3D_option=ob_gen3D_option)

        # make an internal rdkit mol from the conformational ensemble
        self.mol = get_rdkit_mol(self.elements, self.conformer_coordinates,
                                 self.connectivity_matrix, self.charges)

        # fetch identification variables
        self.can = smiles
        self.inchi = Chem.MolToInchi(self.mol)
        self.inchikey = Chem.MolToInchiKey(self.mol)

        # add configuration info
        self.max_num_conformers = num_conf
        self.conformer_engine = engine

        # add charge and spin
        self.charge = sum(self.charges)
        self.spin = Descriptors.NumRadicalElectrons(self.mol) + 1
