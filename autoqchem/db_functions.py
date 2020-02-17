import logging
import pybel

import numpy as np
import pandas as pd
import pymongo

from autoqchem.helper_classes import config
from autoqchem.helper_functions import add_numbers_to_repeated_items

logger = logging.getLogger(__name__)

desc_presets = ['global', 'min_max', 'substructure', 'transitions']
desc_presets_long = ['Global', 'Min Max Atomic', 'Substructure Atomic', "Excited State Transitions"]
conf_options = ['boltzmann', 'max', 'min', 'mean', 'std', 'any']
conf_options_long = ['Boltzman Average', 'Lowest Energy Conformer', 'Highest Energy Conformer', 'Arithmetic Average',
                     'Standard Deviation', 'Random']


def db_connect() -> pymongo.collection.Collection:
    """Create a connection to the database and return the table (Collection).

    :return: pymongo.collection.Collection
    """

    db = pymongo.MongoClient(config['mongoDB']['host'],
                             username=config['mongoDB']['user'],
                             password=config['mongoDB']['password'],
                             port=config['mongoDB']['port'])

    return db['autoqchem']['dft_descriptors']


def db_check_exists(can, gaussian_config) -> list:
    """Check if a molecule is already present in the database with the same Gaussian config (theory, basis_sets, etc.)

    :param can: canonical smiles
    :type can: str
    :param gaussian_config: gaussian config dictionary
    :type gaussian_config: dict
    :return: list of tags that are associated with the molecule
    """

    table = db_connect()
    results = table.distinct("metadata", {"can": can, "metadata.gaussian_config": gaussian_config, })
    return [r['tag'] for r in results]


def db_select_molecules(tag="", substructure="", db_query={}) -> pd.DataFrame:
    """Get a summary frame of records in the database for a given tag

    :param tag: metadata.tag of the db records
    :type tag: str
    :param substructure: substructure SMARTS string
    :type substructure: str
    :param db_query: additional MongoDB type query, use only if familiar with MongoDB query language and \
    the database record structure
    :return: pandas.core.frame.DataFrame
    """

    table = db_connect()

    query = {"$and": [{'metadata.tag': tag} if tag else {}, db_query]}

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
        try:
            pattern = pybel.Smarts(substructure)
            agg['pybel_mol'] = agg['can'].map(lambda can: pybel.readstring("smi", can))
            agg = agg[agg['pybel_mol'].map(lambda mol: bool(pattern.findall(mol)))]
            agg = agg.drop('pybel_mol', axis=1)
        except IOError:
            logger.warning(f"Substructure '{substructure}' is an invalid SMARTS pattern.")

    return agg


def descriptors(tag, presets, conf_option, substructure="") -> dict:
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

    mol_df = db_select_molecules(tag=tag, substructure=substructure)
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


def descriptors_from_can(can, tag, conf_option) -> dict:
    """Get and average descriptors for a given can string and metadata.tag

    :param can: canonical smiles
    :type can: str
    :param tag: metadata.tag of the db records
    :type tag: str
    :param conf_option: conformer averaging option: 'boltzmann' (Boltzmann average), \
     'max' (conformer with highest weight), 'mean' (arithmetic average), 'min' (conformer with smallest weight), \
     'any' (any single conformer), 'std' (std dev. over conformers)
    :type conf_option: str
    :return: dict
    """

    table = db_connect()
    items = table.find({"$and": [{'metadata.tag': tag}, {'can': can}]}, {'_id': 1})
    ids = [item['_id'] for item in items]
    return descriptors_from_list_of_ids(ids, conf_option=conf_option)


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

    rec = {}
    if conf_option in ['boltzmann', 'mean', 'std']:
        # fetch db records for these _ids
        cursor = table.find({"_id": {"$in": ids}})
        recs = [_pandatize_record(record) for record in cursor]
        rec.update({"can": recs[0]['can'], "labels": recs[0]['labels']})

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
