from autoqchem.helper_classes import *
from autoqchem.helper_functions import *

logger = logging.getLogger(__name__)


class gaussian_input_generator(object):
    """generator of gaussian input files"""

    def __init__(self, molecule):
        """configure for this molecule"""

        self.directory = f"{config['gaussian']['work_dir']}/{molecule.fs_name}"
        self.molecule = molecule

        if self.molecule.heavy_elements:
            self.basis_set = config['gaussian']['generic_basis_set']
        else:
            self.basis_set = config['gaussian']['light_basis_set']

    def create_gaussian_files(self, workflow=gaussian_workflows.equilibrium):
        """write gaussian files for each conformation with defined options"""

        # prepare directory for gaussian files
        cleanup_directory_files(self.directory, types=["gjf"])
        os.makedirs(self.directory, exist_ok=True)

        # resources configuration
        n_processors = min(config['gaussian']['max_processors'],
                           self.molecule.mol.NumAtoms() // config['gaussian']['atoms_per_processor'])
        ram = n_processors * config['gaussian']['ram_per_processor']
        resource_block = f"%nprocshared={n_processors}\n%Mem={ram}GB\n"

        logger.info(f"Generating Gaussian input files for {self.molecule.mol.NumConformers()} conformations.")

        for conf_id in range(self.molecule.mol.NumConformers()):
            # set conformer
            conf_name = f"{self.molecule.name}_conf_{conf_id}"
            fs_conf_name = f"{self.molecule.fs_name}_conf_{conf_id}"

            # coordinates block
            geom_np_array = self.molecule.get_initial_geometry(conf_id).astype(str).values
            coords_block = "\n".join(map(" ".join, geom_np_array))

            # create the gaussian input file
            self.__generate_file(workflow,
                                 conf_name,
                                 fs_conf_name,
                                 resource_block,
                                 coords_block,
                                 self.molecule.light_elements,
                                 self.molecule.heavy_elements,
                                 self.molecule.mol.GetTotalCharge(),
                                 self.molecule.mol.GetTotalSpinMultiplicity())

    def __generate_file(self, workflow, name, fs_name, resource_block, coords_block,
                        light_elements, heavy_elements, charge, multiplicity):

        heavy_block = ""
        if self.molecule.heavy_elements:
            heavy_block += f"{' '.join(light_elements + ['0'])}\n"
            heavy_block += f"{config['gaussian']['light_basis_set']}\n****\n"
            heavy_block += f"{' '.join(heavy_elements + ['0'])}\n"
            heavy_block += f"{config['gaussian']['heavy_basis_set']}\n****\n"
            heavy_block += f"\n"
            heavy_block += f"{' '.join(heavy_elements + ['0'])}\n"
            heavy_block += f"{config['gaussian']['heavy_basis_set']}\n"

        output = ""

        # loop through the tasks in the workflow and create input file
        for i, task in enumerate(workflow.value):

            if i == 0:  # first task is special, coordinates follow
                output += resource_block
                output += f"%Chk={fs_name}_{i}.chk\n"
                output += f"# {task['route'](self.basis_set)}\n\n"
                output += f"{name}\n\n"
                output += f"{charge} {multiplicity}\n"
                output += f"{coords_block.strip()}\n"
                output += f"\n"
            else:
                output += "\n--Link1--\n"
                output += resource_block
                output += f"%Oldchk={fs_name}_{i - 1}.chk\n"
                output += f"%Chk={fs_name}_{i}.chk\n"
                output += f"# {task['route'](self.basis_set)}\n"
                output += f"\n"

            if self.molecule.heavy_elements:
                output += heavy_block

        output += f"\n\n"

        file_path = f"{self.directory}/{fs_name}.gjf"
        with open(file_path, "w") as file:
            file.write(output)

        logger.info(f"Generated a Gaussian input file in {file_path}")
