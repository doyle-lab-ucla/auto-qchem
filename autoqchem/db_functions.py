import logging

try:
    from openbabel import pybel # ob 3.0.0
except ImportError: # ob 2.4
    import pybel
import re

import numpy as np
import pandas as pd
import pymongo
from bson.objectid import ObjectId
from rdkit import Chem
from rdkit.Chem import rdFMCS

from autoqchem.helper_classes import config
from autoqchem.helper_functions import add_numbers_to_repeated_items

logger = logging.getLogger(__name__)

desc_presets = ['global', 'min_max', 'substructure', 'core', 'labeled', 'transitions']
desc_presets_long = ['Global', 'Min Max Atomic', 'Substructure Atomic', 'Common Core Atomic', 'Labeled Atomic',
                     "Excited State Transitions"]
conf_options = ['boltzmann', 'max', 'min', 'mean', 'std', 'any']
conf_options_long = ['Boltzman Average', 'Lowest Energy Conformer', 'Highest Energy Conformer', 'Arithmetic Average',
                     'Standard Deviation', 'Random']


class InconsistentLabelsException(Exception):
    """Raised when a set of molecules is inconsistently labeled"""
    pass


def db_connect(collection=None) -> pymongo.collection.Collection:
    """Create a connection to the database and return the table (Collection).

    :return: pymongo.collection.Collection
    """

    cli = pymongo.MongoClient(config['mongoDB']['host'],
                              username=config['mongoDB']['user'],
                              password=config['mongoDB']['password'],
                              port=config['mongoDB']['port'])
    if collection is None:
        return cli['autoqchem']
    else:
        return cli['autoqchem'][collection]


def db_upload_molecule(can, tags, metadata, weights, conformations, logs) -> ObjectId:
    """Upload single molecule to DB and all child objects tags, features
    and log files for its conformations"""

    db = db_connect()
    mols_coll = db['molecules']
    tags_coll = db['tags']

    # create molecule record and insert it
    mol_data = {'can': can, 'metadata': metadata}
    ret = mols_coll.insert_one(mol_data)
    mol_id = ret.inserted_id

    # insert tag record
    for tag in tags:
        tags_coll.insert_one({'tag': tag, 'molecule_id': mol_id, 'can': can})

    for weight, conformation, log in zip(weights, conformations, logs):
        db_upload_conformation(mol_id, can, weight, conformation, log, check_mol_exists=False)
    return mol_id


def db_upload_conformation(mol_id, can, weight, conformation, log, check_mol_exists=True):
    """Upload single conformation features and log file to DB, requires a molecule
    to be present"""

    db = db_connect()
    # check if the molecule with a given id exists in the DB
    mols_coll = db["molecules"]
    if check_mol_exists:
        assert mols_coll.find_one({'_id': mol_id}) is not None

    # connect to features and logs collections
    feats_coll = db['qchem_descriptors']
    logs_coll = db['log_files']

    data = {'molecule_id': mol_id, 'weight': weight, 'can': can}

    # update with descriptors
    data.update(conformation)

    # db insertion
    feats_coll.insert_one(data)
    logs_coll.insert_one({'molecule_id': mol_id, 'log': log, 'can': can})


def db_delete_molecule(mol_id):
    """Delete molecule from DB, cascade all child objects: tags, features and log files"""

    db = db_connect()
    if isinstance(mol_id, str):
        mol_id = ObjectId(mol_id)

    print(mol_id)
    db['qchem_descriptors'].delete_many({"molecule_id": mol_id})  # features
    db['log_files'].delete_many({"molecule_id": mol_id})  # log files
    db['tags'].delete_many({"molecule_id": mol_id})  # tags
    db['molecules'].delete_one({"_id": mol_id})  # molecule itself


def db_select_molecules(cls=None, subcls=None, type=None, subtype=None, tags=[], substructure="") -> pd.DataFrame:
    """Get a summary frame of molecules in the database

    :param tags: a list of tags of the db records (if multiple an 'OR' is taken)
    :type tags: list
    :param substructure: substructure SMARTS string
    :type substructure: str
    :return: pandas.core.frame.DataFrame
    """

    db = db_connect()
    tags_coll = db['tags']
    mols_coll = db['molecules']
    feats_coll = db['qchem_descriptors']

    tags_cur = tags_coll.find({'tag': {'$in': tags}} if tags else {})
    tags_df = pd.DataFrame(tags_cur)

    filter = {}
    if cls != "" and cls is not None:
        filter['metadata.class'] = cls
    if subcls != "" and subcls is not None:
        filter['metadata.subclass'] = subcls
    if type != "" and type is not None:
        filter['metadata.type'] = type
    if subtype != "" and subtype is not None:
        filter['metadata.subtype'] = subtype

    filter['_id'] = {'$in': tags_df.molecule_id.tolist()}

    mols_cur = mols_coll.find(filter)
    mols_df = pd.DataFrame(mols_cur)
    if 'name' not in mols_df.columns:
        mols_df['name'] = None

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
    df = pd.concat([grouped['metadata', 'molecule_id', 'name'].first(),
                    grouped['tag'].apply(list)], axis=1).reset_index().drop('metadata_str', axis=1)

    # fetch ids and weights
    feats_cur = feats_coll.find({'molecule_id': {'$in': df.molecule_id.tolist()}},
                                {'_id': 1, 'weight': 1, 'molecule_id': 1})
    feats_df = pd.DataFrame(feats_cur)
    feats_df = feats_df.groupby('molecule_id').agg(list).reset_index()
    feats_df = feats_df.rename(columns={'_id': '_ids', 'weight': 'weights'})

    # merge into df
    df = df.merge(feats_df, on='molecule_id')
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

    db = db_connect()
    mols_coll = db['molecules']
    tags_coll = db['tags']
    mol_id = mols_coll.find_one({"can": can,
                                 "metadata.gaussian_config": gaussian_config,
                                 "metadata.max_num_conformers": max_num_conformers},
                                {"molecule_id": 1})

    exists, tags = False, []
    if mol_id is not None:
        exists = True
        tags = tags_coll.distinct('tag', {'molecule_id': ObjectId(mol_id['_id'])})
    return exists, tags


def descriptors(cls, subcls, type, subtype, tags, presets, conf_option, substructure="") -> dict:
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

    mol_df = db_select_molecules(cls=cls, subcls=subcls, type=type, subtype=subtype,
                                 tags=tags, substructure=substructure)
    # TODO making DB queries inside a loop is very inefficient, this code should be reorganized to use single query
    descs_df = mol_df.set_index('can')['_ids'].map(lambda l: descriptors_from_list_of_ids(l, conf_option=conf_option))

    data = {}

    if 'global' in presets:
        dg = pd.concat([d['descriptors'] for can, d in descs_df.iteritems()], axis=1, sort=True)
        dg.columns = descs_df.index
        data['global'] = dg.T

    if 'min_max' in presets:
        dmin = pd.concat([d['atom_descriptors'].min() for can, d in descs_df.iteritems()], axis=1, sort=True)
        dmax = pd.concat([d['atom_descriptors'].max() for can, d in descs_df.iteritems()], axis=1, sort=True)
        dmin.columns = descs_df.index
        dmax.columns = descs_df.index
        data['min'] = dmin.T
        data['max'] = dmax.T

    if 'transitions' in presets:
        # select top 3 transitions by oscillation strength
        ts = pd.concat([d['transitions'].sort_values("ES_osc_strength",
                                                     ascending=False).head(10).reset_index(drop=True).unstack()
                        for can, d in descs_df.iteritems()], axis=1, sort=True)
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
        sub_labels = [f"atom{i + 1}" for i in range(len(matches[0]))]

        # create a frame with descriptors large structure in one column, and substructure match
        # indices in the second column
        tmp_df = descs_df.to_frame('descs')
        tmp_df['matches'] = matches

        for i, label in enumerate(sub_labels):
            # data[label] = pd.concat([row['descs']['atom_descriptors'].loc[row['matches'][i]]
            #                         for c, row in tmp_df.iterrows()], axis=1)
            to_concat = []
            for c, row in tmp_df.iterrows():
                atom_descs = row['descs']['atom_descriptors']
                atom_descs['labels'] = row['descs']['labels']
                to_concat.append(atom_descs.iloc[row['matches'][i]])
            data[label] = pd.concat(to_concat, axis=1, sort=True)
            data[label].columns = descs_df.index
            data[label] = data[label].T

    if 'core' in presets:
        cans = mol_df['can'].tolist()
        rd_mols = {can: Chem.MolFromSmiles(can) for can in cans}

        # occasionally rdkit cannot create a molecule from can that openbabel can
        # this is typically due to dative bonds, dative
        for can, rd_mol in rd_mols.items():
            if rd_mol is None:
                logger.warning(f"Molecule with can: {can} cannot be constructed directly by rdkit.")
                rd_mols[can] = Chem.MolFromSmarts(can)  # create it from smarts

        # run MCS if there is more than 1 molecule
        if len(rd_mols) > 1:
            core_smarts = rdFMCS.FindMCS(list(rd_mols.values())).smartsString
        else:  # otherwise use the entire smiles as smarts string
            core_smarts = Chem.MolToSmarts(list(rd_mols.values())[0])

        # create an rdkit smarts
        core = Chem.MolFromSmarts(core_smarts)

        # get the first match, if multiple substructure matches exist
        matches = {can: rd_mols[can].GetSubstructMatches(core)[0] for can in cans}
        matches = pd.Series(matches).map(list)

        # create a frame with descriptors large structure in one column, and substructure match
        # indices in the second column
        tmp_df = descs_df.to_frame('descs')
        tmp_df['matches'] = matches

        # fetch atom labels for this smarts using the first molecule
        row = tmp_df.iloc[0]
        row_labels = pd.Series(row['descs']['labels'])
        row_labels = row_labels[~row_labels.str.startswith('H')]  # need to remove hydrogens
        sub_labels = row_labels.iloc[row['matches']].tolist()
        sub_labels = add_numbers_to_repeated_items(sub_labels)

        for i, label in enumerate(sub_labels):
            to_concat = []
            for c, row in tmp_df.iterrows():
                atom_descs = row['descs']['atom_descriptors']
                atom_descs['labels'] = row['descs']['labels']
                atom_descs = atom_descs[~atom_descs['labels'].str.startswith("H")]  # need to remove hydrogens
                to_concat.append(atom_descs.iloc[row['matches'][i]])
            data[label] = pd.concat(to_concat, axis=1, sort=True)
            data[label].columns = descs_df.index
            data[label] = data[label].T

    if 'labeled' in presets:
        # extract the positions of the labeled atoms in the atom lists for each molecule
        labels = descs_df.map(lambda d: [re.sub("\D", "", l) for l in d['labels']])
        labels = labels.map(lambda ls: [(index, l) for index, l in enumerate(ls) if l])
        labels = labels.map(lambda ls: sorted(ls, key=lambda l: l[1]))

        # verify that the atomic labels are consistent across all molecules
        atom_numbers = labels.map(lambda ls: [l[1] for l in ls])
        atom_numbers_dedup = atom_numbers.map(tuple).drop_duplicates()
        if len(atom_numbers_dedup) == 1:
            matches = labels.map(lambda ls: [l[0] for l in ls])

            # create a frame with descriptors large structure in one column, and substructure match
            # indices in the second column
            tmp_df = descs_df.to_frame('descs')
            tmp_df['matches'] = matches

            for i, label in enumerate(atom_numbers_dedup.iloc[0]):
                label = 'A' + label
                data[label] = pd.concat([row['descs']['atom_descriptors'].loc[row['matches'][i]]
                                         for c, row in tmp_df.iterrows()], axis=1, sort=True)
                data[label].columns = descs_df.index
                data[label] = data[label].T
        else:
            logger.warning("Atomic labels are inconsistent. Not all molecules have the same set of labeled atoms")
            raise InconsistentLabelsException
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

    if record['mode_vectors'] is not None:
        record['mode_vectors'] = pd.DataFrame(record['mode_vectors'])
        record['mode_vectors']['atom_idx'] = list(range(len(record['labels']))) * 3 * record['modes'].shape[0]
        record['mode_vectors'] = record['mode_vectors'].set_index(['mode_number', 'axis', 'atom_idx']).unstack(
            ['mode_number', 'axis'])
        record['mode_vectors'] = record['mode_vectors'].droplevel(0, axis=1).astype(float)
    else:
        record['mode_vectors'] = pd.DataFrame()

    return record
