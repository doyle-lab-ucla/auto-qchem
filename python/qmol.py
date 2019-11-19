import os

from openbabel import pybel
from rdkit.Chem import GetPeriodicTable

from python.enums import *
from python.gaussian_file_generator import *

class qmol(object):
    """Class that holds a single molecule information"""

    pt = GetPeriodicTable()

    def __init__(self, input_string, input_type=input_types.smiles, n_conformers=None):
        """initialize the molecule"""

        self.name = input_string
        self.directory = f"{config.data_dir}/{self.name}"
        self.mol = pybel.readstring(input_type.value, input_string)
        self.mol.addh()
        self.mol.make3D()  # make use of pybel api (harder to do with OBBuilder directly)
        self.mol.localopt()

        if n_conformers is not None:
            self.__generate_conformers(n_conformers)

    def __generate_conformers(self, n_conformers):
        """generate conformations"""

        confSearch = pybel.ob.OBConformerSearch()
        confSearch.Setup(self.mol.OBMol, n_conformers)
        confSearch.Search()
        confSearch.GetConformers(self.mol.OBMol)

    def __group_light_heavy_elements(self):
        """group elements into light and heavy for this molecule"""

        atomic_nums = set(atom.atomicnum for atom in self.mol.atoms)
        heavy_elements = [self.pt.GetElementSymbol(n) for n in atomic_nums if n > config.max_light_atomic_number]
        light_elements = [self.pt.GetElementSymbol(n) for n in atomic_nums if n <= config.max_light_atomic_number]
        return heavy_elements, light_elements

    def create_gaussian_files(self, workflow=gaussian_workflows.equilibrium):
        """write gaussian files for each conformation with defined options"""

        # make directory for files
        os.makedirs(self.directory, exist_ok=True)

        # group elements by heavy and light
        heavy_elements, light_elements = self.__group_light_heavy_elements()
        has_heavy = len(heavy_elements) > 0

        # configure gaussian output writer
        gfg = gaussian_file_generator(has_heavy, workflow, self.directory)

        # resources configuration
        n_processors = max(20, len(self.mol.atoms) // config.atoms_per_processor)
        RAM = n_processors * config.ram_per_processor
        resource_block = f"%nprocshared={n_processors}\n%Mem={RAM}GB\n"

        conv = pybel.ob.OBConversion()
        conv.SetOutFormat("xyz")

        for conf_id in range(self.mol.OBMol.NumConformers()):
            # set conformer
            self.mol.OBMol.SetConformer(conf_id)
            conf_name = f"{self.name}_conf_{conf_id}"

            # coordinates block, use babel to convert and skip first 2 lines
            coords_string = conv.WriteString(self.mol.OBMol)
            coords_block = "\n".join(coords_string.split("\n")[2:])

            gfg.generate_file(conf_name,
                              resource_block,
                              coords_block,
                              light_elements,
                              heavy_elements,
                              self.mol.charge,
                              self.mol.spin)

    def cleanup_gaussian_files(self):
        """cleanup gaussian input/output files"""
        for file in os.listdir(self.directory):
            os.remove(os.path.join(self.directory, file))
        os.rmdir(self.directory)