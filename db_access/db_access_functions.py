import pymongo
import pandas as pd
import numpy as np
import logging

from rdkit import Chem, Geometry
from rdkit.Chem import rdFMCS

from autoqchem.db_functions import db_select_molecules

df = db_select_molecules(tags=['primary_amides_janssen2'])
print(df[['molecule_id', '_ids']].head(10))