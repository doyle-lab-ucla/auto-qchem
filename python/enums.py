import enum
from config import *

class input_types(enum.Enum):
    """enumeration of input types"""

    smiles = "smi"
    chemdraw = "cdx"


class gaussian_workflows(enum.Enum):
    """enumeration of gaussian workflows"""

    equilibrium = {
        "opt": lambda basis_set: f"opt=CalcFc {config.theory}/{basis_set} scf=xqc",
        "freq": lambda basis_set: f"freq {config.theory}/{basis_set} volume NMR pop=NPA density=current Geom=AllCheck Guess=Read",
        "td": lambda basis_set: f"TD(NStates=10, Root=1) {config.theory}/{basis_set} volume pop=NPA density=current Geom=AllCheck Guess=Read",
    }
    transition_state = {
        "opt": lambda basis_set: f"opt=(calcfc,ts,noeigentest) scf=xqc {config.theory}/{basis_set}",
        "freq": lambda basis_set: f"freq {config.theory}/{basis_set} volume NMR pop=NPA density=current Geom=AllCheck Guess=Read"
    }