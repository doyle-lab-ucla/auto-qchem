from autoqchem.molecule import molecule
from autoqchem.slurm_manager import *

logging.basicConfig(level=logging.INFO)

s1 = "c1c([2c](ccc1)[1P]([3CH]1CCCCC1)"
s2 = "CN(C)c1ccc(cc1)P(C2CCCCC2)C3CCCCC3"
s3 = "n1([2c](ccc1)[1P]([3CH]1CCCCC1)[4CH]1CCCCC1)c1ccccc1"
s4 ="[1P]([2CH3])([4CH3])[3CH3]"
s5 = "c1ccc(cc1)[As](c2ccccc2)c3ccccc3"
s6 = "[1P]([4CH2]C)([2CH2]C)[3CH2]C"
s7 = "[B-](F)(F)(F)F.C1CCC(CC1)[PH+](C2CCCCC2)C3CCCCC3"

sm = slurm_manager()


for s in [s5][:]:
    m = molecule(s, input_types.string, input_formats.smiles)
    sm.create_jobs_for_molecule(m)
