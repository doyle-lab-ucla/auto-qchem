import pymongo
import pandas as pd
import numpy as np
import logging

from rdkit import Chem, Geometry
from rdkit.Chem import rdFMCS

from autoqchem.db_functions import db_select_molecules

# fetch conformer data with dataset tag
def conformers(tags=[]):
    df = db_select_molecules(tags=['primary_amides_janssen2'])
    df = df[['can', 'molecule_id', '_ids', 'weights', 'num_conformers']]

    # export molecules first
    mols = df[['can', 'molecule_id', 'weights']]
    mols.to_csv('mols.csv')

if __name__ == '__main__':
    conformers()