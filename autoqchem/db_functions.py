import logging

import pandas as pd
import pymongo
from rdkit import Chem

from autoqchem.helper_classes import config

logger = logging.getLogger(__name__)


def db_connect() -> pymongo.collection.Collection:
    """Create a connection to the database and return the table (Collection).

    :return: pymongo.collection.Collection
    """

    db = pymongo.MongoClient(config['mongoDB']['host'],
                             username=config['mongoDB']['user'],
                             password=config['mongoDB']['password'],
                             port=config['mongoDB']['port'])

    return db['autoqchem']['dft_descriptors']


def db_can_summary(tag="", substructure="", db_query={}) -> pd.DataFrame:
    """Get a summary frame of records in the database for a given tag

    :param tag: metadata.tag of the db records
    :type tag: str
    :return: pandas.core.frame.DataFrame
    """

    table = db_connect()

    query = {"$and": [{'metadata.tag': tag}, db_query]}

    # fetch records
    cursor = table.find(query, {'can': 1,
                                'weight': 1,
                                'metadata.gaussian_config': 1,
                                'metadata.max_num_conformers': 1,
                                'metadata.tag': 1})
    records = list(cursor)

    columns = ['can', 'tag', 'DFT_functional', 'DFT_basis_set', 'num_conformers', 'max_num_conformers', 'weights',
               '_ids']
    if not records:
        return pd.DataFrame(columns=columns)

    # reorganize metadata records
    for record in records:
        record['tag'] = record['metadata']['tag']
        record['max_num_conformers'] = record['metadata']['max_num_conformers']
        record['DFT_functional'] = record['metadata']['gaussian_config']['theory']
        record['DFT_basis_set'] = record['metadata']['gaussian_config']['light_basis_set']
        del record['metadata']

    records_df = pd.DataFrame(records)
    records_df = records_df.sort_values(['can', 'weight'], ascending=[True, False])

    # fetch only unique can-tag-config combinations
    agg = records_df.groupby(['can', 'tag', 'DFT_functional', 'DFT_basis_set']).agg({
        'max_num_conformers': ['size', 'first'],
        'weight': lambda x: ["%.6f" % w for w in list(x)],
        '_id': list,
    }).reset_index()
    agg.columns = columns

    # substructure search
    if substructure:
        pattern = Chem.MolFromSmiles(substructure)
        if pattern is not None:
            agg['rdkit_mol'] = agg['can'].map(lambda can: Chem.MolFromSmiles(can))
            agg = agg[agg['rdkit_mol'].map(lambda mol: mol.HasSubstructMatch(pattern))]
            agg = agg.drop('rdkit_mol', axis=1)
        else:
            logger.warning(f"Pattern '{substructure}' could not be initialized by rdkit. "
                           f"Please verify the smiles string.")

    return agg


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


def descriptors_from_can(can, tag, conf_option='boltzmann') -> dict:
    """Get and average descriptors for a given can string and metadata.tag

    :param can: canonical smiles
    :type can: str
    :param tag: metadata.tag of the db records
    :type tag: str
    :param conf_option: conformer averaging options: 'boltzmann' (Boltzmann average), \
     'max' (conformer with highest weight), 'mean' (arithmetic average), 'min' (conformer with smallest weight), \
     'any' (any single conformer), 'std' (std dev. over conformers)
    :return: dict
    """

    table = db_connect()
    items = table.find({"$and": [{'metadata.tag': tag}, {'can': can}]}, {'_id': 1})
    ids = [item['_id'] for item in items]
    return descriptors_from_list_of_ids(ids, conf_option=conf_option)


def descriptors_from_list_of_ids(ids, conf_option='boltzmann') -> dict:
    """Get and average descriptors using a list of db ids.

    :param ids: list of db id's that correspond to a given molecule
    :type ids: list
    :param conf_option: conformer averaging options: 'boltzmann' (Boltzmann average), \
     'max' (conformer with highest weight), 'mean' (arithmetic average), 'min' (conformer with smallest weight), \
     'any' (any single conformer), 'std' (std dev. over conformers)
    :type conf_option: str
    :return: dict
    """

    # validate conf_option
    if conf_option not in ['boltzmann', 'min', 'max', 'mean', 'std', 'any']:
        logger.warning(f"Provided option for reweighting {conf_option} is not one of 'boltzmann',"
                       f" 'min', 'max', 'mean', 'std', 'any'.")
        return

    # connect to db
    table = db_connect()

    # fetch db _ids and weights and can
    cursor = table.find({"_id": {"$in": ids}}, {'weight': 1, 'can': 1})
    recs = pd.DataFrame(cursor).sort_values('weight', ascending=False)

    # assert that all ids come from the same can, and that weights sum to 1.
    assert len(recs.can.unique()) == 1
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
        return _pandatize_record(table.find_one({"_id": _id}))

    if conf_option in ['boltzmann', 'mean', 'std']:
        # fetch db records for these _ids
        cursor = table.find({"_id": {"$in": ids}})
        recs = [_pandatize_record(record) for record in cursor]

        rec = {"can": recs[0]['can'], "labels": recs[0]['labels']}
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
