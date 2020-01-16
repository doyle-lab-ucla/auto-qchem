import hashlib
import pickle
from contextlib import suppress

import pymongo

from autoqchem.gaussian_input_generator import *
from autoqchem.helper_functions import *
from autoqchem.openbabel_conversions import *

logger = logging.getLogger(__name__)


class slurm_manager(object):
    cache_file = f"{config['slurm']['work_dir']}/slurm_manager.pkl"

    def __init__(self, user, host):
        """load jobs under management in constructor"""

        self.jobs = {}  # jobs under management

        # load jobs under management from cache_file (suppress exceptions, no file, empty file, etc.)
        with suppress(Exception):
            with open(self.cache_file, 'rb') as cf:
                self.jobs = pickle.load(cf)

        self.host = host
        self.user = user
        self.remote_dir = f"/home/{self.user}/gaussian"
        self.connection = None

    def __cache(self):
        """save jobs under management and cleanup empty directories"""

        with open(self.cache_file, 'wb') as cf:
            pickle.dump(self.jobs, cf)

        cleanup_empty_dirs(config['slurm']['work_dir'])

    def connect(self) -> None:
        """connect to remote server"""

        create_new_connection = False
        # check if connection already exists
        if self.connection is not None:
            # check if it went stale
            if not self.connection.is_connected:
                logger.info(f"Connection got disconnected, reconnecting.")
                self.connection.close()
                self.connection = None
                create_new_connection = True
        else:
            logger.info(f"Creating connection to {self.host} as {self.user}")
            create_new_connection = True
        if create_new_connection:
            self.connection = ssh_connect(self.host, self.user)
            self.connection.run(f"mkdir -p {self.remote_dir}")
            logger.info(f"Connected to {self.host} as {self.user}.")

    def create_jobs_for_molecule(self, molecule, workflow=gaussian_workflows.equilibrium) -> None:
        """generate slurm jobs for a given molecule"""

        molecule_work_dir = f"{config['slurm']['work_dir']}/{molecule.fs_name}"

        # create gaussian files
        gig = gaussian_input_generator(molecule, molecule_work_dir)
        tasks = gig.create_gaussian_files(workflow)

        # create slurm files
        for gjf_file in glob.glob(f"{molecule_work_dir}/*.gjf"):

            base_name = os.path.basename(os.path.splitext(gjf_file)[0])
            slurm_manager.create_slurm_file_from_gaussian_file(base_name, molecule_work_dir)
            # create job structure
            job = slurm_job(can=molecule.can,
                            conformation=int(base_name.split("_conf_")[1]),
                            tasks=tasks,
                            job_id=-1,  # job_id (not assigned yet)
                            directory=gig.directory,  # filesystem path
                            base_name=base_name,  # filesystem basename
                            status=slurm_status.created)  # status

            # create a key for the job
            key = hashlib.md5((job.can + str(job.conformation) + ','.join(map(str, job.tasks))).encode()).hexdigest()

            # check if a job like that already exists
            if key in self.jobs:  # a job like that is already present:
                logger.warning(f"A job with exactly the same parameters, molecule {job.can}, conformation "
                               f"{job.conformation}, workflow {job.tasks} already exists. "
                               f"Not creating a duplicate")
                continue

            self.jobs[key] = job  # add job to bag
        self.__cache()

    @staticmethod
    def create_slurm_file_from_gaussian_file(base_name, directory) -> None:
        """generate a single slurm submission file based on the gaussian input file"""

        # get information from gaussian file needed for submission
        with open(f"{directory}/{base_name}.gjf") as f:
            file_string = f.read()

        n_processors = re.search("nprocshared=(.*?)\n", file_string).group(1)

        output = ""
        output += f"#!/bin/bash\n"
        output += f"#SBATCH -N 1\n" \
                  f"#SBATCH --ntasks-per-node={n_processors}\n" \
                  f"#SBATCH -t {config['slurm']['wall_time']}\n" \
                  f"#SBATCH -C haswell\n\n"
        output += f"input={base_name}\n\n"
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

        sh_file_path = f"{directory}/{base_name}.sh"
        with open(sh_file_path, "w") as f:
            f.write(output)
        convert_crlf_to_lf(sh_file_path)
        logger.debug(f"Created a Slurm job file in {sh_file_path}")

    def get_jobs(self, status) -> dict:
        """get a list of jobs by status"""

        return {name: job for name, job in self.jobs.items() if job.status == status}

    def remove_jobs(self, jobs) -> None:
        """remove jobs"""

        for name, job in jobs.items():
            logger.debug(f"Removing job {name}.")
            os.remove(f"{job.directory}/{job.base_name}.sh")  # slurm file
            os.remove(f"{job.directory}/{job.base_name}.gjf")  # gaussian file
            if os.path.exists(f"{job.directory}/{job.base_name}.log"):
                os.remove(f"{job.directory}/{job.base_name}.log")  # log file
            del self.jobs[name]
        self.__cache()

    def get_job_stats(self, split_by_can=False) -> pd.DataFrame:
        """count jobs with each status under management"""

        df = pd.DataFrame([[v.status.value, v.can] for v in self.jobs.values()], columns=['status', 'can'])
        if split_by_can:
            return df.groupby(['status', 'can']).size().unstack(level=1).fillna(0).astype(int)
        else:
            return df.groupby('status').size().to_frame('jobs')

    def submit_jobs(self) -> None:
        """submit jobs of a given status"""

        jobs = self.get_jobs(slurm_status.created)
        logger.info(f"Submitting {len(jobs)} jobs.")
        self.__submit_jobs_from_jobs_dict(jobs)

    def __submit_jobs_from_jobs_dict(self, jobs) -> None:
        """submit jobs"""

        # check if there are any jobs to be submitted
        if jobs:
            # get or create connection
            self.connect()

            # check if jobs are in status created or failed
            for name, job in jobs.items():
                # copy .sh and .gjf file to remote_dir
                self.connection.put(f"{job.directory}/{job.base_name}.sh", self.remote_dir)
                self.connection.put(f"{job.directory}/{job.base_name}.gjf", self.remote_dir)

                with self.connection.cd(self.remote_dir):
                    ret = self.connection.run(f"sbatch {self.remote_dir}/{job.base_name}.sh", hide=True)
                    job.job_id = re.search("job\s*(\d+)\n", ret.stdout).group(1)
                    job.status = slurm_status.submitted
                    job.n_submissions = job.n_submissions + 1
                    logger.info(f"Submitted job {name}, job_id: {job.job_id}.")

            self.__cache()

    def squeue(self) -> None:
        """run squeue command"""

        self.connect()
        self.connection.run(f"squeue -u {self.user}")

    def retrieve_jobs(self) -> None:
        """retrieve finished jobs"""

        ids_to_check = [j.job_id for j in self.get_jobs(slurm_status.submitted).values()]
        if not ids_to_check:
            logger.info(f"There are no jobs submitted to cluster. Nothing to retrieve.")
            return

        # get or create connection
        self.connect()

        # retrieve job ids that are running on the server
        ret = self.connection.run(f"squeue -u {self.user} -o %A,%T", hide=True)
        user_running_ids = [s.split(',')[0] for s in ret.stdout.splitlines()[1:]]
        running_ids = [id for id in user_running_ids if id in ids_to_check]
        finished_ids = [id for id in ids_to_check if id not in running_ids]

        logger.info(f"There are {len(running_ids)} running jobs, {len(finished_ids)} finished jobs.")

        # get finished jobs
        finished_jobs = {name: job for name, job in self.jobs.items() if job.job_id in finished_ids}
        failed_jobs, done_jobs = {}, {}

        if finished_jobs:
            logger.info(f"Retrieving log files of finished jobs.")
            for name, job in finished_jobs.items():
                try:
                    log_file = self.connection.get(f"{self.remote_dir}/{job.base_name}.log",
                                                   local=f"{job.directory}/{job.base_name}.log")
                    with open(log_file.local) as f:
                        log = f.read()
                    job.n_success_steps = len(re.findall("Normal termination", log))

                    n_steps = len(job.tasks)
                    if job.n_success_steps == n_steps:
                        job.status = slurm_status.done
                        done_jobs[name] = job
                    else:
                        job.status = slurm_status.failed
                        failed_jobs[name] = job
                        logger.warning(f"Job {job.base_name} has {job.n_success_steps} out of {n_steps} completed."
                                       f"The job needs to be resubmitted.")
                except FileNotFoundError:  # TODO rare occurence, not exactly sure what exception is thrown here
                    job.status = slurm_status.failed
                    failed_jobs[name] = job
                    logger.warning(f" Job {job.base_name} run out of time.")

                # clean up files on the remote site
                with self.connection.cd(self.remote_dir):
                    self.connection.run(f"rm slurm-{job.job_id}.out")
                    self.connection.run(f"rm {job.base_name}.*")

            self.__cache()
            logger.info(f"{len(done_jobs)} jobs finished successfully (all Gaussian steps finished normally)."
                        f" {len(failed_jobs)} jobs failed.")

    def resubmit_failed_jobs(self):
        """resubmit jobs that failed"""

        failed_jobs = self.get_jobs(slurm_status.failed)

        if not failed_jobs:
            logger.info("There are no failed jobs to resubmit.")

        for key, job in failed_jobs.items():

            # put a limit on resubmissions
            if job.n_submissions >= 3:
                logger.warning(f"Job {key} has been already failed 3 times, not submitting again.")
                continue

            job_log = f"{job.directory}/{job.base_name}.log"
            job_gjf = f"{job.directory}/{job.base_name}.gjf"
            job_sh = f"{job.directory}/{job.base_name}.sh"

            # option 1 replace geometry
            try:
                le = gaussian_log_extractor(job_log)
                # old coords block
                with open(job_gjf, "r") as f:
                    file_string = f.read()
                old_coords_block = re.search(f"\w+\s+({float_or_int_regex})"
                                             f"\s+({float_or_int_regex})"
                                             f"\s+({float_or_int_regex}).*?\n\n",
                                             file_string, re.DOTALL).group(0)

                # new coords block
                coords = le.geom[list('XYZ')]
                coords.insert(0, 'Atom', le.labels)
                coords_block = "\n".join(map(" ".join, coords.values.astype(str))) + "\n\n"

                # make sure they are the same length and replace
                assert len(old_coords_block.splitlines()) == len(coords_block.splitlines())
                file_string = file_string.replace(old_coords_block, coords_block)
                with open(job_gjf, "w") as f:
                    f.write(file_string)

            # option 2, log file is empty or does not have geometry
            except Exception:

                # open sh file and increase the runtime by factor of 2
                with open(job_sh, "r") as f:
                    file_string = f.read()

                # fetch time
                h, m, s = re.search("#SBATCH -t (\d+):(\d+):(\d+)", file_string).groups()
                h = str((int(h) + 1) * 2 - 1)  # increase time by factor of 2
                # substitute time
                file_string = re.sub("#SBATCH -t (\d+):(\d+):(\d+)", f"#SBATCH -t {h}:{m}:{s}", file_string)

                with open(job_sh, "w") as f:
                    f.write(file_string)
                convert_crlf_to_lf(job_sh)

        self.__submit_jobs_from_jobs_dict(failed_jobs)

    def upload_done_jobs_to_db(self, tag, symmetry=True) -> None:
        """Deduplicate done jobs molecule by molecule. The conformers
         are removed if their RMSD is below thershold.
        This is a postprocessing step."""

        done_jobs = self.get_jobs(slurm_status.done)
        if not done_jobs:
            logger.info("There are no jobs in done status. Exitting.")
            return

        # create db connection
        db = pymongo.MongoClient(config['mongoDB']['host'],
                                 username=config['mongoDB']['user'],
                                 password=config['mongoDB']['password'],
                                 port=config['mongoDB']['port'])

        if not yes_or_no("Please make sure all jobs for given molecules have finished. The conformers will be "
                         "deduplicated and uploaded to the db. Would you like to continue?"):
            return

        logger.info(f"Deduplicating conformers if RMSD < {config['gaussian']['conformer_RMSD_threshold']}.")

        done_jobs_df = pd.DataFrame([job.__dict__ for job in done_jobs.values()], index=done_jobs.keys())

        for (can, tasks), keys in done_jobs_df.groupby(["can", "tasks"]).groups.items():
            mols = [OBMol_from_done_slurm_job(done_jobs[key]) for key in keys]
            duplicates = deduplicate_list_of_OBMols(mols, symmetry)
            logger.info(f"Molecule {can} has {len(duplicates)} / {len(keys)} duplicate conformers.")
            can_keys_to_keep = [key for i, key in enumerate(keys) if i not in duplicates]
            self.upload_can_to_db(db, can, tasks, can_keys_to_keep, tag)

        # cleanup
        db.close()
        if yes_or_no("Ok to remove log files for done jobs?"):
            self.remove_jobs(done_jobs)

    def upload_can_to_db(self, db, can, tasks, keys, tag) -> None:
        """uploading done jobs to database"""

        # make sure tag is a non-empty string
        assert isinstance(tag, str)
        assert len(tag.strip()) > 0

        # loop over the conformers
        conformations = []
        for key in keys:
            # fetch job, verify that there are not can issues (just in case)
            job = self.jobs[key]
            assert job.can == can

            # extract descriptors for this conformer from log file
            log = f"{job.directory}/{job.base_name}.log"
            le = gaussian_log_extractor(log)
            le.extract_descriptors()
            # add descriptors to conformations list
            conformations.append(le.data())

        # compute weights
        free_energies = np.array(
            [Hartree_in_kcal_per_mol * c['descriptors']['G'] for c in conformations])  # in kcal_mol
        free_energies -= free_energies.min()  # to avoid huge exponentials
        weights = np.exp(-free_energies / (k_in_kcal_per_mol_K * T))
        weights /= weights.sum()

        # fetch only significant conformers
        min_weight = config['gaussian']['conformer_weight_threshold']
        big_weights_indices = np.where(weights > min_weight)[0]
        logger.info(f"Conformer weights: {weights}")
        logger.info(f"Keeping only {len(big_weights_indices)} conformers with weights > {min_weight}.")

        # filter out only those conformations
        conformations = np.array(conformations)[big_weights_indices]
        weights = weights[big_weights_indices]
        weights /= weights.sum()

        for weight, conformation in zip(weights, conformations):
            data = {'can': can,
                    'metadata': {'gaussian_config': config['gaussian'], 'gaussian_tasks': tasks, 'tag': tag},
                    'weight': weight}
            # update with descriptors
            data.update(conformation)

            # db insertion
            db['autoqchem']['dft_descriptors'].insert_one(data)
        logger.info(f"Uploaded descriptors to DB for smiles: {can}, number of conformers: {len(conformations)}.")
