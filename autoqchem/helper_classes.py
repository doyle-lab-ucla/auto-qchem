import enum
from dataclasses import dataclass

import yaml

config = yaml.safe_load(open("config.yml"))


@enum.unique
class input_types(enum.IntEnum):
    """enumeration of input types"""

    file = enum.auto()
    string = enum.auto()


@enum.unique
class gaussian_workflows(enum.Enum):
    """enumeration of gaussian workflows"""

    equilibrium = [
        {"desc": "Geometry Optimization",
         "route": lambda basis_set: f"opt=CalcFc {config['gaussian']['theory']}/{basis_set} scf=xqc"},
        {"desc": "Frequency Calculation",
         "route": lambda basis_set: f"freq {config['gaussian']['theory']}/{basis_set}"
                                    f" volume NMR pop=NPA density=current Geom=AllCheck Guess=Read"},
        {"desc": "Time Dependent Calcualtion",
         "route": lambda basis_set: f"TD(NStates=10, Root=1) {config['gaussian']['theory']}/{basis_set} "
                                    f"volume pop=NPA density=current Geom=AllCheck Guess=Read"},
    ]
    transition_state = [
        {"desc": "Geometry Optimization",
         "route": lambda basis_set: f"opt=(calcfc,ts,noeigentest) scf=xqc {config['gaussian']['theory']}/{basis_set}"},
        {"desc": "Frequency Calculation",
         "route": lambda basis_set: f"freq {config['gaussian']['theory']}/{basis_set}"
                                    f" volume NMR pop=NPA density=current Geom=AllCheck Guess=Read"},
    ]
    test = [
        {"desc": "Hartree-Fock", "route": lambda basis_set: f"{config['gaussian']['theory']}/{basis_set}"},
    ]


@enum.unique
class slurm_status(enum.Enum):
    """enumeration for slurm job status"""

    created = 'created'  # job files have been created
    submitted = 'submitted'  # jobs has been submitted to slurm
    finished = 'finished'  # job finished running with slurm (failed or success)
    failed = 'failed'  # job failed
    success = 'success'  # job finished successfully and is ready for retrieval


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
