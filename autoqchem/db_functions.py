import logging
import pybel

import numpy as np
import pandas as pd
import pymongo
from bson.objectid import ObjectId

from autoqchem.helper_classes import config
from autoqchem.helper_functions import add_numbers_to_repeated_items

logger = logging.getLogger(__name__)

desc_presets = ['global', 'min_max', 'substructure', 'transitions']
desc_presets_long = ['Global', 'Min Max Atomic', 'Substructure Atomic', "Excited State Transitions"]
conf_options = ['boltzmann', 'max', 'min', 'mean', 'std', 'any']
conf_options_long = ['Boltzman Average', 'Lowest Energy Conformer', 'Highest Energy Conformer', 'Arithmetic Average',
                     'Standard Deviation', 'Random']


def db_connect(collection="molecules") -> pymongo.collection.Collection:
    """Create a connection to the database and return the table (Collection).

    :return: pymongo.collection.Collection
    """

    db = pymongo.MongoClient(config['mongoDB']['host'],
                             username=config['mongoDB']['user'],
                             password=config['mongoDB']['password'],
                             port=config['mongoDB']['port'])

    return db['autoqchem'][collection]


def db_select_molecules(tags=[], substructure="") -> pd.DataFrame:
    """Get a summary frame of molecules in the database

    :param tags: a list of tags of the db records (if multiple an 'OR' is taken)
    :type tags: list
    :param substructure: substructure SMARTS string
    :type substructure: str
    :return: pandas.core.frame.DataFrame
    """

    tags_coll = db_connect('tags')
    mols_coll = db_connect('molecules')
    feats_coll = db_connect('qchem_descriptors')

    tags_cur = tags_coll.find({'tag': {'$in': tags}} if tags else {})
    tags_df = pd.DataFrame(tags_cur)

    mols_cur = mols_coll.find({'_id': {'$in': tags_df.molecule_id.tolist()}})
    mols_df = pd.DataFrame(mols_cur)

    if substructure:
        pattern = pybel.Smarts(substructure)
        mols_df['pybel_mol'] = mols_df['can'].map(lambda can: pybel.readstring("smi", can))
        mols_df = mols_df[mols_df['pybel_mol'].map(lambda mol: bool(pattern.findall(mol)))]
        mols_df = mols_df.drop('pybel_mol', axis=1)

    # merge tags in an outer way
    df = pd.merge(mols_df, tags_df, how='outer', left_on='_id', right_on='molecule_id', suffixes=('', '_tag'))

    # make tags into a list of tags
    df['metadata_str'] = df['metadata'].map(repr)
    grouped = df.groupby(['can', 'metadata_str'])
    # groupby tags
    df = pd.concat([grouped['metadata', 'molecule_id'].first(),
                    grouped['tag'].apply(list)], axis=1).reset_index().drop('metadata_str', axis=1)

    # fetch ids
    df['_ids'] = df['molecule_id'].map(lambda mid: [item['_id'] for item in feats_coll.find(
        {'molecule_id': ObjectId(mid)}, {'_id': 1})
                                                    ])
    df['num_conformers'] = df['_ids'].map(len)
    return df


def db_check_exists(can, gaussian_config, max_num_conformers) -> tuple:
    """Check if a molecule is already present in the database with the same Gaussian config (theory, basis_sets, etc.)

    :param can: canonical smiles
    :type can: str
    :param gaussian_config: gaussian config dictionary
    :type gaussian_config: dict
    :return: exists(bool), list of tags that are associated with the molecule if it exists
    """

    mols_coll = db_connect(collection='molecules')
    tags_coll = db_connect(collection='tags')
    mol_id = mols_coll.find_one({"can": can,
                                 "metadata.gaussian_config": gaussian_config,
                                 "metadata.max_num_conformers": max_num_conformers},
                                {"molecule_id": 1})

    exists, tags = False, []
    if mol_id is not None:
        exists = True
        tags = tags_coll.distinct('tag', {'molecule_id': ObjectId(mol_id['_id'])})
    return exists, tags


def descriptors(tags, presets, conf_option, substructure="") -> dict:
    """Retrieve DFT descriptors from the database

    :param tag: metadata.tag of the db records
    :type tag: str
    :param presets: list of descriptor presets from 'global' (molecule level descriptors), \
    'min_max' (min and max for each atomic descriptor across the molecule), 'substructure' \
    (atomic descriptors for each atom in the substructure)
    :type presets: list
    :param conf_option: conformer averaging option: 'boltzmann' (Boltzmann average), \
    'max' (conformer with highest weight), 'mean' (arithmetic average), 'min' (conformer with smallest weight), \
    'any' (any single conformer), 'std' (std dev. over conformers)
    :type conf_option: str
    :param substructure: substructure SMARTS string
    :type substructure: str
    :return:
    """

    # don't bother with extraction if there are not presets nor conf_option
    if not presets or not conf_option:
        logger.warning(f"One of options 'presets' or 'conf_option' is empty. Not extracting.")
        return {}

    # check that presets are ok
    if not all(p in desc_presets for p in presets):
        logger.warning(f"One of the presets in {presets} is not from allowed list {desc_presets}. Not extracting.")
        return {}

    # check that conf option is ok
    if conf_option not in conf_options:
        logger.warning(f"Conf_option {conf_option} is not one of the allowed options {conf_options}. Not extracting.")
        return {}

    mol_df = db_select_molecules(tags=tags, substructure=substructure)
    descs_df = mol_df.set_index('can')['_ids'].map(lambda l: descriptors_from_list_of_ids(l, conf_option=conf_option))

    data = {}

    if 'global' in presets:
        dg = pd.concat([d['descriptors'] for can, d in descs_df.iteritems()], axis=1)
        dg.columns = descs_df.index
        data['global'] = dg.T

    if 'min_max' in presets:
        dmin = pd.concat([d['atom_descriptors'].min() for can, d in descs_df.iteritems()], axis=1)
        dmax = pd.concat([d['atom_descriptors'].max() for can, d in descs_df.iteritems()], axis=1)
        dmin.columns = descs_df.index
        dmax.columns = descs_df.index
        data['min'] = dmin.T
        data['max'] = dmax.T

    if 'transitions' in presets:
        # select top 3 transitions by oscillation strength
        ts = pd.concat([d['transitions'].sort_values("ES_osc_strength",
                                                     ascending=False).head(3).reset_index(drop=True).unstack()
                        for can, d in descs_df.iteritems()], axis=1)
        ts.index = ts.index.map(lambda i: "_".join(map(str, i)))
        ts.columns = descs_df.index
        data['transitions'] = ts.T

    if 'substructure' in presets and substructure:
        sub = pybel.Smarts(substructure)
        # these matches are numbered from 1, so subtract one from them
        matches = descs_df.index.map(lambda c: sub.findall(pybel.readstring("smi", c))[0])
        matches = matches.map(lambda x: (np.array(x) - 1).tolist())

        # fetch atom labels for this smarts using the first molecule
        sub_labels = pd.Series(descs_df.iloc[0]['labels']).loc[matches[0]].tolist()
        sub_labels = add_numbers_to_repeated_items(sub_labels)

        # create a frame with descriptors large structure in one column, and substructure match
        # indices in the second column
        tmp_df = descs_df.to_frame('descs')
        tmp_df['matches'] = matches

        for i, label in enumerate(sub_labels):
            data[label] = pd.concat([row['descs']['atom_descriptors'].loc[row['matches'][i]]
                                     for c, row in tmp_df.iterrows()], axis=1)
            data[label].columns = descs_df.index
            data[label] = data[label].T

    return data


def descriptors_from_list_of_ids(ids, conf_option='max') -> dict:
    """Get and average descriptors using a list of db ids.

    :param ids: list of db id's that correspond to a given molecule
    :type ids: list
    :param conf_option: conformer averaging option: 'boltzmann' (Boltzmann average), \
     'max' (conformer with highest weight), 'mean' (arithmetic average), 'min' (conformer with smallest weight), \
     'any' (any single conformer), 'std' (std dev. over conformers)
    :type conf_option: str
    :return: dict
    """

    # check that conf option is ok
    if conf_option not in conf_options:
        logger.warning(f"Conf_option {conf_option} is not one of the allowed options {conf_options}. Not extracting.")
        return {}

    # connect to db
    feats_coll = db_connect("qchem_descriptors")

    # fetch db _ids and weights and can
    cursor = feats_coll.find({"_id": {"$in": ids}}, {'weight': 1, 'molecule_id': 1})
    recs = pd.DataFrame(cursor).sort_values('weight', ascending=False)

    # assert that all ids come from the same can, and that weights sum to 1.
    assert len(recs.molecule_id.unique()) == 1
    assert abs(recs.weight.sum() - 1.) < 1e-6

    # set trivial option for the case with only one conformation
    if len(recs) == 1:
        conf_option = 'any'

    # single conf options
    if conf_option in ['min', 'max', 'any']:
        if conf_option == 'max':
            _id = recs['_id'].iloc[0]
        elif conf_option == 'min':
            _id = recs['_id'].iloc[-1]
        else:
            _id = recs['_id'].sample(1).iloc[0]

        # return pandatized record for a chosen id
        return _pandatize_record(feats_coll.find_one({"_id": _id}))

    rec = {}
    if conf_option in ['boltzmann', 'mean', 'std']:
        # fetch db records for these _ids
        cursor = feats_coll.find({"_id": {"$in": ids}})
        recs = [_pandatize_record(record) for record in cursor]
        rec.update({"labels": recs[0]['labels']})

        keys_to_reweight = ['descriptors', 'atom_descriptors', 'modes', 'transitions']

        for key in keys_to_reweight:
            if conf_option == 'boltzmann':
                dfs = pd.concat(r[key] * r['weight'] for r in recs)
                rec[key] = dfs.groupby(dfs.index, sort=False).sum()
            if conf_option in ['mean', 'std']:
                dfs = pd.concat(r[key] for r in recs)
                if conf_option == 'mean':
                    rec[key] = dfs.groupby(dfs.index, sort=False).mean()
                elif conf_option == 'std':
                    rec[key] = dfs.groupby(dfs.index, sort=False).std()

    return rec


def _pandatize_record(record) -> dict:
    """Convert json structures to pandas structures for an individual
    db record of a single conformation.

    :param record: db record of a single conformation
    :return: dict
    """

    del record['descriptors']['stoichiometry']

    record['descriptors'] = pd.Series(record['descriptors']).astype(float)
    record['modes'] = pd.DataFrame(record['modes']).astype(float)
    record['transitions'] = pd.DataFrame(record['transitions']).astype(float)
    record['atom_descriptors'] = pd.DataFrame(record['atom_descriptors']).astype(float)

    record['mode_vectors'] = pd.DataFrame(record['mode_vectors'])
    record['mode_vectors']['atom_idx'] = list(range(len(record['labels']))) * 3 * record['modes'].shape[0]
    record['mode_vectors'] = record['mode_vectors'].set_index(['mode_number', 'axis', 'atom_idx']).unstack(
        ['mode_number', 'axis'])
    record['mode_vectors'] = record['mode_vectors'].droplevel(0, axis=1).astype(float)

    return record
