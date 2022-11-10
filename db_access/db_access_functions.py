import pymongo
import pandas as pd
import numpy as np
import logging
import pathlib
from tqdm import tqdm
import json

logger = logging.getLogger(__name__)

from rdkit import Chem, Geometry
from rdkit.Chem import rdFMCS

from autoqchem.db_functions import db_select_molecules, db_connect

Hartree_in_kcal_per_mol = 627.5


# fetch all conformer data with dataset tag

def all_conformer_data(tags=None):
    """

    """
    df = db_select_molecules(tags=tags)
    df = df[['can', 'molecule_id', '_ids', 'weights', 'num_conformers']]

    # export molecules with metadata first
    mols = df[['can', 'molecule_id', 'weights', '_ids']]
    mols['weights'] = mols['weights'].apply(lambda x: str(sorted(x, reverse=True)))  # reverse sort weights: lowest energy = highest weight
    #mols.to_csv('mols.csv')

    mol_cans = df['can'].to_list()
    mols_coll = db_connect('molecules')
    tags_coll = db_connect('tags')

    # loop through molecules with can_smiles
    for c in tqdm(mol_cans):
        if tags:  # recommended; prefilter molecules that belong to a specific tag
            mol_ids = [record['molecule_id'] for record in tags_coll.find({'tag': {'$in': tags}},
                                                                          {'molecule_id': 1, '_id': 0})]
            m = mols_coll.find_one({'can': c, '_id': {"$in": mol_ids}})
        else:
            m = mols_coll.find_one({'can': c})

        if m is None:
            logger.warning(f"Molecule {c} not found.")
            return None, None

        # find all features
        feats_coll = db_connect("qchem_descriptors")
        feats = feats_coll.find({'molecule_id': m['_id']},
                                {'_id': 0, 'descriptors.G': 1, 'labels': 1, 'atom_descriptors': 1,})
        feats = list(feats)

        # sort with energies
        energies = np.array([f['descriptors']['G'] * Hartree_in_kcal_per_mol for f in feats])
        order = energies.argsort()
        energies = energies[order]
        elements = m['elements']

        # fetch atom descriptors and sort with energies
        atom_descriptors = np.array([np.array([f['atom_descriptors']]).T for f in feats])
        atom_descriptors = atom_descriptors[order,:]

        # save individual conformer data
        num_conf = atom_descriptors.shape[0]
        for ii in range(num_conf):
            data = pd.DataFrame.from_dict(atom_descriptors[ii][0])
            data.index = elements
            fp = pathlib.Path.cwd() / 'data' / c
            fp.mkdir(parents=True, exist_ok=True)
            fn = str(ii)+'.csv'
            data.to_csv(fp / fn)


if __name__ == '__main__':
    all_conformer_data(tags=['primary_amides_janssen2'])
