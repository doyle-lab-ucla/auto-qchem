from autoqchem.helper_classes import *
from autoqchem.helper_functions import *

logger = logging.getLogger(__name__)


class gaussian_input_generator(object):
    """Generator of gaussian input files class"""

    def __init__(self, molecule, workflow_type, directory):
        """Initialize input generator for a given molecule.

        :param molecule: molecule object
        :type molecule: molecule
        :param workflow_type: Gaussian workflow type, allowed types are: 'equilibrium' or 'transition_state'
        :type workflow_type: str
        :param directory: local directory to store input files
        :type directory: str
        """

        self.directory = directory
        self.molecule = molecule

        if self.molecule.heavy_elements:
            basis_set = config['gaussian']['generic_basis_set']
        else:
            basis_set = config['gaussian']['light_basis_set']

        # create gaussian tasks tuple
        theory = config['gaussian']['theory']

        if workflow_type == "equilibrium":
            self.tasks = (
                f"opt=CalcFc {theory}/{basis_set} scf=xqc",
                f"freq {theory}/{basis_set} volume NMR pop=NPA density=current Geom=AllCheck Guess=Read",
                f"TD(NStates=10, Root=1) {theory}/{basis_set} volume pop=NPA density=current Geom=AllCheck Guess=Read"
            )
        elif workflow_type == "transition_state":
            self.tasks = (
                f"opt=(calcfc,ts,noeigentest) scf=xqc {theory}/{basis_set}",
                f"freq {theory}/{basis_set} volume NMR pop=NPA density=current Geom=AllCheck Guess=Read"
            )
        elif workflow_type == "test":
            self.tasks = (
                f"{theory}/{basis_set}"
            )
        else:
            raise ValueError(f"Not supported gaussian job type {workflow_type}. "
                             f"Allowed types are: equilibrium, transition_state.")

    def create_gaussian_files(self) -> None:
        """Create the actual gaussian files for each conformer of the molecule."""

        # prepare directory for gaussian files
        cleanup_directory_files(self.directory, types=["gjf"])
        os.makedirs(self.directory, exist_ok=True)

        # resources configuration
        n_processors = min(config['slurm']['max_processors'],
                           self.molecule.mol.NumAtoms() // config['slurm']['atoms_per_processor'])
        ram = n_processors * config['slurm']['ram_per_processor']
        resource_block = f"%nprocshared={n_processors}\n%Mem={ram}GB\n"

        logger.info(f"Generating Gaussian input files for {self.molecule.mol.NumConformers()} conformations.")

        for conf_id in range(self.molecule.mol.NumConformers()):
            # set conformer
            conf_name = f"{self.molecule.can}_conf_{conf_id}"
            fs_conf_name = f"{self.molecule.fs_name}_conf_{conf_id}"

            # coordinates block
            geom_df = self.molecule.get_geometry(conf_id)
            geom_df['Iso_String'] = geom_df['Isotope'].astype(str).map(lambda s: f"(Iso={s})")
            geom_df['Atom'] += geom_df['Iso_String'].where(geom_df['Isotope'] > 0, '')
            geom_np_array = geom_df[['Atom', 'X', 'Y', 'Z']].astype(str).values
            coords_block = "\n".join(map(" ".join, geom_np_array))

            # create the gaussian input file
            self._generate_file(self.tasks,
                                conf_name,
                                fs_conf_name,
                                resource_block,
                                coords_block,
                                self.molecule.light_elements,
                                self.molecule.heavy_elements,
                                self.molecule.mol.GetTotalCharge(),
                                self.molecule.mol.GetTotalSpinMultiplicity())

    def _generate_file(self, tasks, name, fs_name, resource_block, coords_block,
                       light_elements, heavy_elements, charge, multiplicity) -> None:
        """

        :param tasks: tuple of Gaussian tasks
        :param name:  conformation name
        :param fs_name: filesystem conformation name
        :param resource_block: resource block for the Gaussian input file
        :param coords_block: coordinates block for the Gaussian input file
        :param light_elements: list of light elements of the molecule
        :param heavy_elements: list of heavy elements of the molecule
        :param charge: molecule charge
        :param multiplicity: molecule multiplicity
        """

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
        for i, task in enumerate(tasks):

            if i == 0:  # first task is special, coordinates follow
                output += resource_block
                output += f"%Chk={fs_name}_{i}.chk\n"
                output += f"# {task}\n\n"
                output += f"{name}\n\n"
                output += f"{charge} {multiplicity}\n"
                output += f"{coords_block.strip()}\n"
                output += f"\n"
            else:
                output += "\n--Link1--\n"
                output += resource_block
                output += f"%Oldchk={fs_name}_{i - 1}.chk\n"
                output += f"%Chk={fs_name}_{i}.chk\n"
                output += f"# {task}\n"
                output += f"\n"

            if self.molecule.heavy_elements:
                output += heavy_block

        output += f"\n\n"

        file_path = f"{self.directory}/{fs_name}.gjf"
        with open(file_path, "w") as file:
            file.write(output)

        logger.debug(f"Generated a Gaussian input file in {file_path}")
