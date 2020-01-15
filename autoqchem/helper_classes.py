import enum
import os
from dataclasses import dataclass

import yaml

config = yaml.safe_load(open(os.path.join(os.path.dirname(__file__), "..", "config.yml")))
k_in_kcal_per_mol_K = 0.0019872041
Hartree_in_kcal_per_mol = 627.5
T = 298


@enum.unique
class input_types(enum.IntEnum):
    """enumeration of input types"""

    file = enum.auto()
    string = enum.auto()


@enum.unique
class gaussian_workflows(enum.Enum):
    """enumeration of gaussian workflows"""

    equilibrium = (
        lambda basis_set: f"opt=CalcFc {config['gaussian']['theory']}/{basis_set} scf=xqc",  # geometry optimization
        lambda basis_set: f"freq {config['gaussian']['theory']}/{basis_set} volume NMR pop=NPA "  # frequency calc
                          f"density=current Geom=AllCheck Guess=Read",
        lambda basis_set: f"TD(NStates=10, Root=1) {config['gaussian']['theory']}/{basis_set} "  # Time-Dependent
                          f"volume pop=NPA density=current Geom=AllCheck Guess=Read",
    )
    transition_state = (
        lambda basis_set: f"opt=(calcfc,ts,noeigentest) scf=xqc {config['gaussian']['theory']}/{basis_set}",  # geometry
        lambda basis_set: f"freq {config['gaussian']['theory']}/{basis_set} volume NMR pop=NPA"  # frequency calc.
                          f" density=current Geom=AllCheck Guess=Read",
    )
    test = (
        lambda basis_set: f"{config['gaussian']['theory']}/{basis_set}",  # Hartree-Fock
    )


@enum.unique
class slurm_status(enum.Enum):
    """enumeration for slurm job status"""

    created = 'created'  # job files have been created
    submitted = 'submitted'  # jobs has been submitted to slurm
    done = 'done'  # job finished successfully (all gaussian steps finished)
    failed = 'failed'  # job failed


@dataclass
class slurm_job:
    """data class for slurm job"""

    # molecule and gaussian config
    can: str
    conformation: int
    tasks: tuple

    # slurm variables
    job_id: int
    directory: str
    base_name: str
    status: slurm_status
    n_submissions: int = 0
    n_success_steps: int = 0
