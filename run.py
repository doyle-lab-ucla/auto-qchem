from python.qmol import qmol
from python.enums import *

m=qmol("Cc1ccccc1P(c2ccccc2C)c3ccccc3C", input_types.smiles, n_conformers=5)
m.create_gaussian_files()
m.cleanup_gaussian_files()