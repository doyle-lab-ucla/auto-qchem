import re
import glob
import pickle
from dataclasses import dataclass

import pandas as pd

from autoqchem.gaussian_file_generator import *
from autoqchem.util import *

logger = logging.getLogger(__name__)


@dataclass
class slurm_job:
    """data class for slurm job"""

    job_id: int
    path: str
    name: str
    checkpoints: list
    status: slurm_status
    n_submissions: int = 0
    n_success_steps: int = 0


class slurm_manager(object):

    cache_file = f"{config['gaussian']['work_dir']}/slurm_manager.pkl"

    def __init__(self):
        """load jobs under management in constructor"""

        self.jobs = {} # jobs under management
        if os.path.exists(self.cache_file):
            with open(self.cache_file, 'rb') as cf:
                self.jobs = pickle.load(cf)
        self.connection = None

    def __del__(self):
        """save jobs under management on destruction"""

        self.disconnect()
        self.__cache()

    def __cache(self):
        """save jobs under management"""

        with open(self.cache_file, 'wb') as cf:
            pickle.dump(self.jobs, cf)

    def connect(self):
        """connect to remote server"""

        self.host = config['slurm']['host']
        self.user = config['slurm']['user']
        logger.info(f"Connecting to {self.host} as {self.user}")
        self.remote_dir = f"/home/{self.user}/gaussian"
        self.connection = ssh_connect(self.host, self.user)
        self.connection.run(f"mkdir -p {self.remote_dir}")
        logger.info(f"Connected to {self.host} as {self.user}.")

    def disconnect(self):
        """disconnect from remote server"""

        if self.connection is not None:
            self.connection.close()
            self.connection = None

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
        self.__cache()

    def __create_jobs_from_directory(self, directory):
        """generate slurm for all the available gaussian files in a directory"""

        # fetch gaussian files
        gjf_files = [f for f in os.listdir(directory) if f.endswith("gjf")]

        # cleanup directory from slurm .sh files
        cleanup_directory_files(directory, types=['sh'])

        for gjf_file in gjf_files:
            self.__create_job_from_gaussian_file(os.path.join(directory, gjf_file), directory)

    def __create_job_from_gaussian_file(self, gjf_file_path, directory):
        """generate a single slurm submission file based on the gaussian input file"""

        gjf_file_name = os.path.basename(gjf_file_path)
        base_file_name = os.path.splitext(gjf_file_name)[0]

        # get information from gaussian file needed for submission
        with open(gjf_file_path) as f:
            file_string = f.read()

        n_processors = re.search("nprocshared=(.*?)\n", file_string).group(1)
        checkpoints = re.findall("Chk=(.*?)\n", file_string)

        output = ""
        output += f"#!/bin/bash\n"
        output += f"#SBATCH -N 1\n" \
                  f"#SBATCH --ntasks-per-node={n_processors}\n" \
                  f"#SBATCH -t {config['slurm']['wall_time']}\n" \
                  f"#SBATCH -C haswell\n\n"
        output += f"input={base_file_name}\n\n"
        output += f"# create scratch directory for the job\n" \
                  f"export GAUSS_SCRDIR=/scratch/${{USER}}/${{SLURM_JOB_ID}}\n" \
                  f"tempdir=${{GAUSS_SCRDIR}}\n" \
                  f"mkdir -p ${{tempdir}}\n\n"
        output += f"# copy input file to scratch directory\n" \
                  f"cp ${{SLURM_SUBMIT_DIR}}/${{input}}.gjf ${{tempdir}}\n\n"
        output += f"# run the code \n" \
                  f"cd ${{tempdir}}\n" \
                  f"g16 ${{input}}.gjf\n\n"
        output += f"# copy output\n" \
                  f"cp ${{input}}.log ${{SLURM_SUBMIT_DIR}}"

        sh_file_path = f"{directory}/{base_file_name}.sh"
        with open(sh_file_path, "w") as f:
            f.write(output)
        convert_crlf_to_lf(sh_file_path)

        job = slurm_job(None,  # job_id
                        directory,  # sh_path
                        base_file_name,  # gjf_path
                        checkpoints,  # checkpoints
                        slurm_status.created)  # status

        self.jobs[base_file_name] = job  # keep dictionary with name as key
        logger.info(f"Created a Slurm job file in {sh_file_path}")

    def get_jobs(self, status):
        """get a list of jobs by status"""

        return {name: job for name, job in self.jobs.items() if job.status == status}

    def get_job_stats(self):
        """count jobs with each status under management"""

        return pd.DataFrame([v.__dict__ for _, v in self.jobs.items()]).groupby('status').size()

    def submit_jobs(self, status=slurm_status.created):
        """submit jobs of a given status"""

        jobs_bag = self.get_jobs(status)
        logger.info(f"Submitting {len(jobs_bag)} jobs with status {status}")
        self.__submit_jobs_from_jobs_bag(self.get_jobs(status))

    def __submit_jobs_from_jobs_bag(self, jobs):
        """submit jobs"""

        # check if connected
        if self.connection is None:
            self.connect()

        # check if jobs are in status created or failed
        for name, job in jobs.items():
            # copy .sh and .gjf file to remote_dir
            sh_file = self.connection.put(f"{job.path}/{job.name}.sh", self.remote_dir)
            gjf_file = self.connection.put(f"{job.path}/{job.name}.gjf", self.remote_dir)

            with self.connection.cd(self.remote_dir):
                ret = self.connection.run(f"sbatch {sh_file.remote}", hide=True)
                job.job_id = re.search("job\s*(\d+)\n", ret.stdout).group(1)
                job.status = slurm_status.submitted
                job.n_submissions = job.n_submissions + 1
                logger.info(f"Submitted job {name}, job_id: {job.job_id}.")

        self.__cache()

    def squeue(self):
        """run squeue command"""

        # check if connected
        if self.connection is None:
            self.connect()
        self.connection.run(f"squeue -u {self.user}")

    def check_jobs(self):
        """check for finished jobs"""

        ids_to_check = [j.job_id for _, j in self.get_jobs(slurm_status.submitted).items()]
        if not ids_to_check:
            logger.info(f"There are no submitted jobs to check.")
            return

        # check if connected
        if self.connection is None:
                self.connect()
        ret = self.connection.run(f"squeue -j {','.join(ids_to_check)} -o %A,%T", hide=True)
        ids_running = [s.split(',')[0] for s in ret.stdout.splitlines()[1:]]
        ids_finished = [id for id in ids_to_check if id not in ids_running]

        logger.info(f"There are {len(ids_running)} running jobs, {len(ids_finished)} finished jobs.")

        # get finished jobs
        finished_jobs = {name: job for name, job in self.jobs.items() if job.job_id in ids_finished}
        failed_jobs, success_jobs = {}, {}

        if finished_jobs:
            logger.info(f"Checking the termination of finished jobs.")
            for name, job in finished_jobs.items():
                log_file = self.connection.get(f"{self.remote_dir}/{job.name}.log",
                                               local=f"{job.path}/{job.name}.log")
                with open(log_file.local) as f:
                    log = f.read()
                job.n_success_steps = len(re.findall("Normal termination", log))

                if job.n_success_steps == len(job.checkpoints):
                    job.status = slurm_status.success
                    success_jobs[name] = job
                else:
                    job.status = slurm_status.failed
                    failed_jobs[name] = job
                    log.warning(f"Job {job.name} has {job.success_steps} out of {len(job.checkpoints)} completed."
                                f"The job needs to be resubmitted.")

                # clean up files on the remote site
                with self.connection.cd(self.remote_dir):
                    self.connection.rm(f"slurm-{job.id}.out")
                    self.connection.rm(f"{job.name}.*")

            self.__cache()
            logger.info(f"{len(success_jobs)} finished successfuly. {len(failed_jobs)} failed.")
