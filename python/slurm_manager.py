import re
import glob

from python.gaussian_file_generator import  *
from python.util import *

logger = logging.getLogger(__name__)


class slurm_manager(object):

    jobs = {}  # bag of jobs under management

    def create_jobs_for_molecule(self, molecule, workflow=gaussian_workflows.equilibrium):
        """generate slurm jobs for a given molecule"""

        directory = f"{config['gaussian']['work_dir']}/{molecule.fs_name}"
        gfg = gaussian_file_generator(molecule)

        if glob.glob(f"{directory}/*.gjf"):
            if yes_or_no(f"Molecule already has gaussian (gjf) input files created "
                         f"in directory: {directory}.\n"
                      f"Shall I use these files?"):
                logger.info(f"Using pre-existing gaussian (gjf) files.")
            else:
                gfg.create_gaussian_files(workflow)
        else:
            gfg.create_gaussian_files(workflow)

        self.__create_jobs_from_directory(directory)


    def __create_jobs_from_directory(self, directory):
        """generate slurm for all the available gaussian files in a directory"""

        # fetch gaussian files
        gjf_files = [f for f in os.listdir(directory) if f.endswith("gjf")]

        # cleanup directory from slurm .sh files
        cleanup_directory_files(directory, types=['sh'])

        for gjf_file in gjf_files:
            self.create_job_from_gaussian_file(os.path.join(directory, gjf_file), directory)


    def create_job_from_gaussian_file(self, gjf_file_path, directory):
        """generate a single slurm submission file based off gaussian input file"""

        gjf_file_name = os.path.basename(gjf_file_path)
        base_file_name = os.path.splitext(gjf_file_name)[0]

        # get information from gaussian file needed for submission
        with open(gjf_file_path) as f:
            file_string = f.read()

        n_processors = re.search("nprocshared=(.*?)\n", file_string).group(1)
        mem = re.search("Mem=(.*?)\n", file_string).group(1)
        checkpoints = re.findall("Chk=(.*?)\n", file_string)    


        output = ""
        output += f"#!/bin/bash\n"
        output += f"#SBATCH -N 1\n"
        output += f"#SBATCH --ntasks-per-node={n_processors}\n"
        output += f"#SBATCH -t {config['slurm']['wall_time']}\n"
        output += f"#SBATCH -C haswell\n\n"
        output += f"mkdir -p /scratch/$USER/job.$$\n"
        output += f"export GAUSS_SCRDIR=/scratch/$USER/job.$$\n\n"
        output += f"g16 {gjf_file_name}\n\n"
        output += f"echo 'Job Complete' > {base_file_name}.done\n\n"
        output += f"rm {' '.join(checkpoints)}"

        file_path = f"{directory}/{base_file_name}.sh"
        with open(file_path, "w") as f:
            f.write(output)
        convert_crlf_to_lf(file_path)

        self.jobs[file_path] = "Created"
        logger.info(f"Created a Slurm job file in {file_path}")

        self.jobs[file_path] = slurm_status.created