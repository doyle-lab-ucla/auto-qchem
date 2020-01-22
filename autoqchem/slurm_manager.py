import hashlib
import pickle
from contextlib import suppress

import appdirs
import pymongo

from autoqchem.gaussian_input_generator import *
from autoqchem.helper_functions import *
from autoqchem.openbabel_functions import *

logger = logging.getLogger(__name__)


class slurm_manager(object):
    """Slurm manager class."""

    def __init__(self, user, host):
        """Initialize slurm manager and load the cache file.

        :param user: username at remote host
        :type user: str
        :param host: remote host name
        :type host: str
        """

        # set workdir and cache file
        self.workdir = appdirs.user_data_dir(appname="autoqchem")
        self.cache_file = os.path.join(self.workdir, "slurm_manager.pkl")
        os.makedirs(self.workdir, exist_ok=True)

        self.jobs = {}  # jobs under management

        # load jobs under management from cache_file (suppress exceptions, no file, empty file, etc.)
        with suppress(Exception):
            with open(self.cache_file, 'rb') as cf:
                self.jobs = pickle.load(cf)

        self.host = host
        self.user = user
        self.remote_dir = f"/home/{self.user}/gaussian"
        self.connection = None

    def connect(self) -> None:
        """Connect to remote host."""

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

    def create_jobs_for_molecule(self, molecule, workflow_type="equilibrium") -> None:
        """Generate slurm jobs for a molecule. Gaussian input files are also generated.

        :param molecule: molecule object
        :type molecule: molecule
        :param workflow_type: Gaussian workflow type, allowed types are: 'equilibrium' or 'transition_state'
        :type workflow_type: str
        """

        # create gaussian files
        molecule_workdir = os.path.join(self.workdir, molecule.fs_name)
        gig = gaussian_input_generator(molecule, workflow_type, molecule_workdir)
        gig.create_gaussian_files()

        # create slurm files
        for gjf_file in glob.glob(f"{molecule_workdir}/*.gjf"):

            base_name = os.path.basename(os.path.splitext(gjf_file)[0])
            slurm_manager._create_slurm_file_from_gaussian_file(base_name, molecule_workdir)
            # create job structure
            job = slurm_job(can=molecule.can,
                            conformation=int(base_name.split("_conf_")[1]),
                            max_num_conformers=gig.molecule.max_num_conformers,
                            tasks=gig.tasks,
                            job_id=-1,  # job_id (not assigned yet)
                            directory=gig.directory,  # filesystem path
                            base_name=base_name,  # filesystem basename
                            status=slurm_status.created,
                            n_submissions=0,
                            n_success_tasks=0)  # status

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

    def submit_jobs(self) -> None:
        """Submit jobs that have status 'created' to remote host."""

        jobs = self.get_jobs(slurm_status.created)
        logger.info(f"Submitting {len(jobs)} jobs.")
        self._submit_jobs_from_jobs_dict(jobs)

    def retrieve_jobs(self) -> None:
        """Retrieve finished jobs from remote host and check which finished succesfully and which failed."""

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
                    job.n_success_tasks = len(re.findall("Normal termination", log))

                    n_steps = len(job.tasks)
                    if job.n_success_tasks == n_steps:
                        job.status = slurm_status.done
                        done_jobs[name] = job
                    else:
                        job.status = slurm_status.failed
                        failed_jobs[name] = job
                        logger.warning(f"Job {job.base_name} has {job.n_success_tasks} out of {n_steps} completed."
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

    def resubmit_failed_jobs(self) -> None:
        """Resubmit jobs that failed. If the job has failed but a log file has been retrieved, then \
        the last geometry will be used for the next submission. If the job didn't finish in time the time required \
        will be increased by a factor of 2. Maximum number of allowed submission of \
        the same job is 3."""

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

                logger.info("Log file retrieved. Substituting last checked geometry in the new input file.")

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

                logger.warning(f"No log file was retrieved from the server. Resubmitting with more time: {h}:{m}:{s}.")

        self._submit_jobs_from_jobs_dict(failed_jobs)

    def upload_done_molecules_to_db(self, tag, RMSD_threshold=0.01, symmetry=True) -> None:
        """Upload done molecules to db. Molecules are considered done when all jobs for a given \
         smiles are in 'done' status. The conformers are deduplicated and uploaded to database using a metadata tag.

        :param tag: metadata tag to use for these molecules in the database
        :type tag: str
        :param RMSD_threshold: RMSD threshold (in Angstroms) to use when deduplicating multiple conformers \
        after Gaussian has found optimal geometry
        :type RMSD_threshold: float
        :param symmetry: if True symmetry is taken into account when comparing molecules in OBAlign(symmetry=True)
        :type symmetry: bool
        """

        done_jobs = self.get_jobs(slurm_status.done)
        if not done_jobs:
            logger.info("There are no jobs in done status. Exitting.")
            return

        # check if there are molecules with all jobs done
        dfj = self.get_job_stats(split_by_can=True)
        dfj_done = dfj[dfj['done'] == dfj.sum(1)]  # only done jobs
        done_cans = dfj_done.index.tolist()

        if not done_cans:
            logger.info("There are no molecules with all jobs done. Exitting.")
            return

        logger.info(f"There are {len(done_cans)} finished molecules {done_cans}.")

        # create db connection
        db = pymongo.MongoClient(config['mongoDB']['host'],
                                 username=config['mongoDB']['user'],
                                 password=config['mongoDB']['password'],
                                 port=config['mongoDB']['port'])

        # select jobs for done molecules
        done_can_jobs = self.get_jobs(can=done_cans)
        jobs_df = pd.DataFrame([job.__dict__ for job in done_can_jobs.values()], index=done_can_jobs.keys())

        logger.debug(f"Deduplicating conformers if RMSD < {RMSD_threshold}.")

        for (can, tasks, max_n_conf), keys in jobs_df.groupby(["can", "tasks", "max_num_conformers"]).groups.items():
            # deduplicate conformers
            mols = [OBMol_from_done_slurm_job(done_jobs[key]) for key in keys]
            duplicates = deduplicate_list_of_OBMols(mols, RMSD_threshold=RMSD_threshold, symmetry=symmetry)
            logger.info(f"Molecule {can} has {len(duplicates)} / {len(keys)} duplicate conformers.")

            # fetch non-duplicate keys
            can_keys_to_keep = [key for i, key in enumerate(keys) if i not in duplicates]
            self._upload_can_to_db(db, can, tasks, can_keys_to_keep, tag, max_n_conf)

        # cleanup
        db.close()
        if yes_or_no("Ok to remove log files for done jobs?"):
            self.remove_jobs(done_can_jobs)

    def _upload_can_to_db(self, db, can, tasks, keys, tag, max_conf) -> None:
        """Uploading single molecule conformers to database.

        :param db: database client
        :type db: pymongo.MongoClient
        :param can: canonical smiles
        :type can: str
        :param tasks: tuple of Gaussian tasks
        :type tasks: tuple
        :param keys: list of keys to the self.jobs dictionary to upload
        :type keys: list
        :param tag: metadata tag
        :type tag: str
        :param max_conf: max number of conformers used for this molecule
        :type max_conf: int
        """

        # check if the tag is properly provided
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
            # add descriptors to conformations list
            conformations.append(le.get_descriptors())

        # compute weights
        free_energies = np.array(
            [Hartree_in_kcal_per_mol * c['descriptors']['G'] for c in conformations])  # in kcal_mol
        free_energies -= free_energies.min()  # to avoid huge exponentials
        weights = np.exp(-free_energies / (k_in_kcal_per_mol_K * T))
        weights /= weights.sum()

        for weight, conformation in zip(weights, conformations):
            data = {'can': can,
                    'metadata': {
                        'gaussian_config': config['gaussian'],
                        'gaussian_tasks': tasks,
                        'tag': tag,
                        'max_num_conformers': max_conf,
                    },
                    'weight': weight}
            # update with descriptors
            data.update(conformation)

            # db insertion
            db['autoqchem']['dft_descriptors'].insert_one(data)
        logger.info(f"Uploaded descriptors to DB for smiles: {can}, number of conformers: {len(conformations)}.")

    def get_jobs(self, status=None, can=None) -> dict:
        """Get a dictionary of jobs, optionally filter by status and canonical smiles.

        :param status: slurm status of the jobs
        :type status: slurm_status
        :param can: canonical smiles of the molecules, single string for one smiles, a list for multiple smiles
        :type can: str or list
        :return: dict
        """

        def match(job, status, can):
            match = True
            if status is not None:
                match = match and job.status.value == status.value
            if can is not None:
                if isinstance(can, str):
                    can = [can]
                match = match and (job.can in can)
            return match

        return {name: job for name, job in self.jobs.items() if match(job, status, can)}

    def get_job_stats(self, split_by_can=False) -> pd.DataFrame:
        """Job stats for jobs currently under management, optionally split by canonical smiles.

        :param split_by_can: if True each canonical smiles will be listed separately
        :type split_by_can: bool
        :return: pandas.core.frame.DataFrame
        """

        df = pd.DataFrame([[v.status.name, v.can] for v in self.jobs.values()], columns=['status', 'can'])
        if split_by_can:
            return df.groupby(['status', 'can']).size().unstack(level=1).fillna(0).astype(int).T
        else:
            return df.groupby('status').size().to_frame('jobs').T

    def remove_jobs(self, jobs) -> None:
        """Remove jobs.

        :param jobs: dictionary of jobs to remove
        :type jobs: dict
        """

        for name, job in jobs.items():
            logger.debug(f"Removing job {name}.")
            os.remove(f"{job.directory}/{job.base_name}.sh")  # slurm file
            os.remove(f"{job.directory}/{job.base_name}.gjf")  # gaussian file
            if os.path.exists(f"{job.directory}/{job.base_name}.log"):
                os.remove(f"{job.directory}/{job.base_name}.log")  # log file
            del self.jobs[name]
        self.__cache()

    def squeue(self) -> None:
        """Run 'squeue -u $user' command on the server."""

        self.connect()
        self.connection.run(f"squeue -u {self.user}")

    def _scancel(self) -> None:
        """Run 'scancel -u $user' command on the server."""

        self.connect()
        self.connection.run(f"scancel -u {self.user}")
        self.remove_jobs(self.get_jobs(status=slurm_status.submitted))

    def _submit_jobs_from_jobs_dict(self, jobs) -> None:
        """Submit jobs to remote host.

        :param jobs: dictionary of jobs to submit
        :type jobs: dict
        """

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

    @staticmethod
    def _create_slurm_file_from_gaussian_file(base_name, directory) -> None:
        """Generate a single slurm submission file based on the Gaussian input file.

        :param base_name: base name of the Gaussian file
        :param directory: directory location of the Gaussian file
        """

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

    def __cache(self) -> None:
        """save jobs under management and cleanup empty directories"""

        with open(self.cache_file, 'wb') as cf:
            pickle.dump(self.jobs, cf)

        cleanup_empty_dirs(self.workdir)
