import enum
import yaml

config = yaml.safe_load(open("config.yml"))

class input_types(enum.Enum):
    """enumeration of input types"""

    file = 1
    string = 2


class input_formats(enum.Enum):
    """enumeration of input formats"""

    smiles = "smi"
    canonical = "can"
    chemdraw = "cdx"


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
        {"desc": "Geoemtry Optimization",
         "route": lambda basis_set: f"opt=(calcfc,ts,noeigentest) scf=xqc {config['gaussian']['theory']}/{basis_set}"},
        {"desc": "Frequency Calculation",
         "route": lambda basis_set: f"freq {config['gaussian']['theory']}/{basis_set}"
                                    f" volume NMR pop=NPA density=current Geom=AllCheck Guess=Read"},
    ]
    test = [
        {"desc": "Hartree-Fock", "route": lambda basis_set: f"{config['gaussian']['theory']}/{basis_set}"},
        {"desc": "Frequency Calculation",
         "route": lambda basis_set: f"freq {config['gaussian']['theory']}/{basis_set}"
                                    f" volume NMR pop=NPA density=current Geom=AllCheck Guess=Read"},
        {"desc": "Time Dependent Calcualtion",
         "route": lambda basis_set: f"TD(NStates=10, Root=1) {config['gaussian']['theory']}/{basis_set} "
                                    f"volume pop=NPA density=current Geom=AllCheck Guess=Read"},
    ]

class slurm_status(enum.Enum):
    """enumeration for slurm job status"""

    created = 1
    submitted = 2
    failed = 3
    finished = 4
