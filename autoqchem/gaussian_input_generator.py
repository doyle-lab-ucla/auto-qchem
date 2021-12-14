from autoqchem.helper_classes import *
from autoqchem.helper_functions import *
from autoqchem.rdkit_utils import *

logger = logging.getLogger(__name__)


class gaussian_input_generator(object):
    """Generator of gaussian input files class"""

    def __init__(self, molecule, workflow_type, directory, theory, solvent, light_basis_set,
                 heavy_basis_set, generic_basis_set, max_light_atomic_number):
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
        # group elements into light and heavy
        light_elements, heavy_elements = get_light_and_heavy_elements(molecule.mol, max_light_atomic_number)
        self.heavy_block = ""

        if heavy_elements:
            basis_set = generic_basis_set
            self.heavy_block += f"{' '.join(light_elements + ['0'])}\n"
            self.heavy_block += f"{light_basis_set}\n****\n"
            self.heavy_block += f"{' '.join(heavy_elements + ['0'])}\n"
            self.heavy_block += f"{heavy_basis_set}\n****\n"
            self.heavy_block += f"\n"
            self.heavy_block += f"{' '.join(heavy_elements + ['0'])}\n"
            self.heavy_block += f"{heavy_basis_set}\n"
        else:
            basis_set = light_basis_set

        if workflow_type == "equilibrium":
            self.tasks = (
                f"opt=CalcFc {theory}/{basis_set} SCRF=(Solvent={solvent}) scf=xqc ",
                f"freq {theory}/{basis_set} SCRF=(Solvent={solvent}) volume NMR pop=NPA density=current Geom=AllCheck Guess=Read",
                f"TD(NStates=10, Root=1) {theory}/{basis_set} SCRF=(Solvent={solvent}) volume pop=NPA density=current Geom=AllCheck Guess=Read"
            )
        elif workflow_type == "transition_state":
            self.tasks = (
                f"opt=(calcfc,ts,noeigentest) scf=xqc {theory}/{basis_set} SCRF=(Solvent={solvent})",
                f"freq {theory}/{basis_set} SCRF=(Solvent={solvent}) volume NMR pop=NPA density=current Geom=AllCheck Guess=Read"
            )
        elif workflow_type == "test":
            self.tasks = (
                f"Opt B3LYP/6-31G** SCRF=(Solvent=TetraHydroFuran) EmpiricalDispersion=GD3",
                f"Freq B3LYP/6-31G** volume NMR pop=NPA density=current Geom=AllCheck Guess=Read "
                f"SCRF=(Solvent=TetraHydroFuran) EmpiricalDispersion=GD3"
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
        n_processors = max(1, min(config['slurm']['max_processors'],
                                  self.molecule.mol.GetNumAtoms() // config['slurm']['atoms_per_processor']))
        ram = n_processors * config['slurm']['ram_per_processor']
        resource_block = f"%nprocshared={n_processors}\n%Mem={ram}GB\n"

        logger.info(f"Generating Gaussian input files for {self.molecule.mol.GetNumConformers()} conformations.")

        for conf_id, conf_coord in enumerate(self.molecule.conformer_coordinates):
            # set conformer
            conf_name = f"{self.molecule.inchikey}_conf_{conf_id}"

            # coordinates block
            geom_np_array = np.concatenate((np.array([self.molecule.elements]).T, conf_coord), axis=1)
            coords_block = "\n".join(map(" ".join, geom_np_array))

            # create the gaussian input file
            self._generate_file(self.tasks,
                                conf_name,
                                resource_block,
                                coords_block,
                                self.molecule.charge,
                                self.molecule.spin)

    def _generate_file(self, tasks, name, resource_block, coords_block, charge, multiplicity) -> None:
        """

        :param tasks: tuple of Gaussian tasks
        :param name:  conformation name
        :param resource_block: resource block for the Gaussian input file
        :param coords_block: coordinates block for the Gaussian input file
        :param light_elements: list of light elements of the molecule
        :param heavy_elements: list of heavy elements of the molecule
        :param charge: molecule charge
        :param multiplicity: molecule multiplicity
        """

        output = ""

        # loop through the tasks in the workflow and create input file
        for i, task in enumerate(tasks):
            if i == 0:  # first task is special, coordinates follow
                output += resource_block
                output += f"%Chk={name}_{i}.chk\n"
                output += f"# {task}\n\n"
                output += f"{name}\n\n"
                output += f"{charge} {multiplicity}\n"
                output += f"{coords_block.strip()}\n"
                output += f"\n"
            else:
                output += "\n--Link1--\n"
                output += resource_block
                output += f"%Oldchk={name}_{i - 1}.chk\n"
                output += f"%Chk={name}_{i}.chk\n"
                output += f"# {task}\n"
                output += f"\n"

            output += self.heavy_block  # this is an empty string if no heavy elements are in the molecule

        output += f"\n\n"

        file_path = f"{self.directory}/{name}.gjf"
        with open(file_path, "w") as file:
            file.write(output)

        logger.debug(f"Generated a Gaussian input file in {file_path}")
