import logging
import re

import pandas as pd
import numpy as np
import pymongo
from bson.objectid import ObjectId
from rdkit import Chem
from rdkit.Chem import rdFMCS

from autoqchem.helper_classes import config, Hartree_in_kcal_per_mol
from autoqchem.helper_functions import add_numbers_to_repeated_items
from autoqchem.rdkit_utils import get_rdkit_mol

logger = logging.getLogger(__name__)

desc_presets = ['global', 'min_max', 'substructure', 'core', 'transitions']
desc_presets_long = ['Global', 'Min Max Atomic', 'Substructure Atomic', 'Common Core Atomic',
                     "Excited State Transitions"]
conf_options = ['boltzmann', 'max', 'min', 'mean', 'std', 'any']
conf_options_long = ['Boltzman Average', 'Lowest Energy Conformer', 'Highest Energy Conformer', 'Arithmetic Average',
                     'Standard Deviation', 'Random']


def db_connect(collection=None) -> pymongo.collection.Collection:
    """Create a connection to the database and return the table (Collection).

    :param collection: database collection name (optional)
    :type collection: str
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


def db_upload_molecule(mol_data, tags, metadata, weights, conformations, logs) -> ObjectId:
    """Upload single molecule to DB and all child objects tags, features and log files for its conformations

    :param mol_data: molecule identity data (inchi, connectivity, etc.)
    :type mol_data: dict
    :param tags: list of tags to assign to this molecule
    :type tags: list
    :param metadata: Gaussian metadata specifying the configuration of the calculation
    :type metadata: dict
    :param weights: conformation weights
    :type weights: list
    :param conformations: list of extracted descriptors for each conformation
    :type conformations: list
    :param logs: list of log files for each conformation
    :type logs: list
    :return: bson.objectid.ObjectId
    """

    db = db_connect()
    mols_coll = db['molecules']
    tags_coll = db['tags']

    # create molecule record and insert it
    mol_data = {'can': mol_data['can'],
                'inchi': mol_data['inchi'],
                'inchikey': mol_data['inchikey'],
                'elements': mol_data['elements'],
                'charges': mol_data['charges'],
                'connectivity_matrix': mol_data['connectivity_matrix'],
                'metadata': metadata}
    ret = mols_coll.insert_one(mol_data)
    mol_id = ret.inserted_id

    # insert tag record
    for tag in tags:
        tags_coll.insert_one({'tag': tag, 'molecule_id': mol_id})

    for weight, conformation, log in zip(weights, conformations, logs):
        db_upload_conformation(mol_id, weight, conformation, log, check_mol_exists=False)
    return mol_id


def db_upload_conformation(mol_id, weight, conformation, log, check_mol_exists=True):
    """Upload single conformation features and log file to DB, requires a molecule to be present

    :param mol_id: molecule id in the DB
    :type mol_id: bson.objectid.ObjectId
    :param weight: conformation weight
    :type weight: float
    :param conformation: extracted descriptors for the conformation
    :type conformation: dict
    :param log: log file for the conformation
    :type log: str
    :param check_mol_exists: check whether the molecule exists in the database, default is TRUE, not recommended to change
    :type check_mol_exists: bool
    """

    db = db_connect()
    # check if the molecule with a given id exists in the DB
    mols_coll = db["molecules"]
    if check_mol_exists:
        assert mols_coll.find_one({'_id': mol_id}) is not None

    # connect to features and logs collections
    feats_coll = db['qchem_descriptors']
    logs_coll = db['log_files']

    data = {'molecule_id': mol_id, 'weight': weight}

    # update with descriptors
    data.update(conformation)

    # db insertion
    feats_coll.insert_one(data)
    try:
        logs_coll.insert_one({'molecule_id': mol_id, 'log': log})
    except pymongo.errors.DocumentTooLarge:
        logger.warning(f"Log file of conformation with weight {weight} too large for DB (limit is 16MB).")


def db_delete_molecule(mol_id):
    """Delete molecule from DB, cascade all child objects: tags, descriptors and log files

    :param mol_id: molecule id in the DB
    :type mol_id: bson.objectid.ObjectId
    """

    db = db_connect()
    if isinstance(mol_id, str):
        mol_id = ObjectId(mol_id)

    print(mol_id)
    db['qchem_descriptors'].delete_many({"molecule_id": mol_id})  # features
    db['log_files'].delete_many({"molecule_id": mol_id})  # log files
    db['tags'].delete_many({"molecule_id": mol_id})  # tags
    db['molecules'].delete_one({"_id": mol_id})  # molecule itself


def db_select_molecules(tags=[], substructure="", smiles="", solvent="ALL",
                        functional="ALL", basis_set="ALL", molecule_ids=[]) -> pd.DataFrame:
    """Get a summary frame of molecules in the database

    :param tags: a list of tags of the db records (if multiple an 'OR' is taken)
    :type tags: list
    :param substructure: substructure SMARTS string
    :type substructure: str
    :param smiles: smiles string
    :type smiles: str
    :param solvent: solvent filter
    :type solvent: str
    :param functional: functional filter
    :type functional: str
    :param basis_set: basis_set filter
    :type basis_set: str
    :param molecule_ids: filter on specific molecule ids in the DB
    :type molecule_ids: list
    :return: pandas.core.frame.DataFrame
    """

    db = db_connect()
    tags_coll = db['tags']
    mols_coll = db['molecules']
    feats_coll = db['qchem_descriptors']

    tags_cur = tags_coll.find({'tag': {'$in': tags}} if tags else {})
    tags_df = pd.DataFrame(tags_cur)

    # filter
    filter = {}
    if molecule_ids:
        filter['_id'] = {'$in': molecule_ids}
    else:
        filter['_id'] = {'$in': tags_df.molecule_id.tolist()}

    if solvent != 'ALL':
        filter['metadata.gaussian_config.solvent'] = re.compile(f"^{re.escape(solvent)}$", re.IGNORECASE)
    if functional != 'ALL':
        filter['metadata.gaussian_config.theory'] = re.compile(f"^{re.escape(functional)}$", re.IGNORECASE)
    if basis_set != 'ALL':
        filter['metadata.gaussian_config.light_basis_set'] = re.compile(f"^{re.escape(basis_set)}$", re.IGNORECASE)
    if smiles != "":
        inchi = Chem.MolToInchi(Chem.MolFromSmiles(smiles))
        filter['inchi'] = inchi

    mols_cur = mols_coll.find(filter)
    mols_df = pd.DataFrame(mols_cur)
    if 'name' not in mols_df.columns:
        mols_df['name'] = None

    if mols_df.empty:
        return mols_df

    if substructure:
        sub = Chem.MolFromSmarts(substructure)
        mols_df['rdmol'] = mols_df['can'].map(Chem.MolFromSmiles)

        # TODO RDKIT molecule creation will fail for dative bonds with metals, this fails silently, not ideal
        mols_df = mols_df.dropna(subset=['rdmol'])

        mols_df = mols_df[mols_df['rdmol'].map(lambda mol: bool(mol.GetSubstructMatches(sub)))]
        mols_df = mols_df.drop('rdmol', axis=1)

        # if substructure filter returned no results
        if mols_df.empty:
            return mols_df

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


def db_check_exists(inchi, gaussian_config, max_num_conformers, conformer_engine) -> tuple:
    """Check if a molecule is already present in the database with the same Gaussian config (function, basis_set, number of conformers, conformer engine)

    :param inchi: inchi of the molecule
    :type inchi: str
    :param gaussian_config: gaussian config
    :type gaussian_config: dict
    :param max_num_conformers: maximum number of conformers to generate
    :type max_num_conformers: int
    :param conformer_engine: conformer engine used for conformer generation (rdkit or openbabel)
    :type conformer_engine: str
    :return: tuple(exists(bool), list of tags that are associated with the molecule if it exists)
    """

    db = db_connect()
    mols_coll = db['molecules']
    tags_coll = db['tags']
    mol_id = mols_coll.find_one({"inchi": inchi,
                                 "metadata.gaussian_config": gaussian_config,
                                 "metadata.max_num_conformers": max_num_conformers,
                                 "metadata.conformer_engine": conformer_engine},
                                {"_id": 1})

    exists, tags = False, []
    if mol_id is not None:
        exists = True
        tags = tags_coll.distinct('tag', {'molecule_id': ObjectId(mol_id['_id'])})
    return exists, tags


def db_get_molecule(inchi, tags=[]):
    """Get an rdkit molecule from DB conformer geometries

    :param inchi: inchi of the molecule
    :type inchi: str
    :param tags: optional list of tags to narrow the search
    :type tags: list
    :return: rdkit.Chem.Mol
    """

    mols_coll = db_connect('molecules')
    tags_coll = db_connect('tags')

    if tags:
        # prefilter molecules that belong to a specific tag
        mol_ids = [record['molecule_id'] for record in tags_coll.find({'tag': {'$in': tags}},
                                                                      {'molecule_id': 1, '_id': 0})]
        m = mols_coll.find_one({'inchi': inchi, '_id': {"$in": mol_ids}})
    else:
        m = mols_coll.find_one({'inchi': inchi})

    if m is None:
        logger.warning(f"Molecule {inchi} not found.")
        return None, None

    return db_get_rdkit_mol(m)


def db_get_rdkit_mol(molecule_record) -> tuple([Chem.Mol, list]):
    """Get an rdkit molecule from DB conformer geometries

    :param molecule_record: molecule record from the DB
    :type molecule_record: dict
    :return: rdkit.Chem.Mol
    """

    feats_coll = db_connect("qchem_descriptors")
    feats = feats_coll.find({'molecule_id': molecule_record['_id']},
                            {'_id': 0, 'descriptors.G': 1, 'labels': 1, 'atom_descriptors.X': 1,
                             'atom_descriptors.Y': 1, 'atom_descriptors.Z': 1})
    feats = list(feats)

    # dummy check
    assert (all(f['labels'] == molecule_record['elements'] for f in feats))

    energies = np.array([f['descriptors']['G'] * Hartree_in_kcal_per_mol for f in feats])
    coords = np.array([np.array([f['atom_descriptors']['X'],
                                 f['atom_descriptors']['Y'], f['atom_descriptors']['Z']]).T for f in feats])

    connectivity = np.reshape(molecule_record['connectivity_matrix'],
                              newshape=(len(molecule_record['elements']), len(molecule_record['elements'])))

    order = energies.argsort()
    rdmol = get_rdkit_mol(molecule_record['elements'], coords[order, :, :],
                          connectivity, molecule_record['charges'])

    return rdmol, energies[order]


def descriptors(tags, presets, conf_option, solvent, functional, basis_set, substructure="") -> dict:
    """Retrieve DFT descriptors from the database

    :param tags: a list of tags of the db records
    :type tags: list
    :param presets: list of descriptor presets from 'global' (molecule level descriptors), \
    'min_max' (min and max for each atomic descriptor across the molecule), 'substructure' \
    (atomic descriptors for each atom in the substructure)
    :type presets: list
    :param conf_option: conformer averaging option: 'boltzmann' (Boltzmann average), \
    'max' (conformer with highest weight), 'mean' (arithmetic average), 'min' (conformer with smallest weight), \
    'any' (any single conformer), 'std' (std dev. over conformers)
    :type conf_option: str
    :param solvent: solvent filter
    :type solvent: str
    :param functional: functional filter
    :type functional: str
    :param basis_set: basis_set filter
    :type basis_set: str
    :param substructure: substructure SMARTS string
    :type substructure: str
    :return: dict
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

    mol_df = db_select_molecules(tags=tags, substructure=substructure, solvent=solvent, functional=functional,
                                 basis_set=basis_set)
    descs_df = descriptors_from_mol_df(mol_df, conf_option)  # this is the heavy query from DB

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
        try:
            ts = pd.concat([d['transitions'].sort_values("ES_osc_strength",
                                                         ascending=False).head(10).reset_index(drop=True).unstack()
                            for can, d in descs_df.iteritems()], axis=1, sort=True)
            ts.index = ts.index.map(lambda i: "_".join(map(str, i)))
            ts.columns = descs_df.index
            data['transitions'] = ts.T
        except KeyError:
            pass

    if 'substructure' in presets or 'core' in presets:
        # make rdmols
        cans = mol_df['can'].tolist()
        rd_mols = {can: Chem.MolFromSmiles(can) for can in cans}

    if 'substructure' in presets and substructure:
        # create an rdkit smarts
        sub = Chem.MolFromSmarts(substructure)

        # get the first match, if multiple substructure matches exist
        matches = {can: rd_mols[can].GetSubstructMatches(sub)[0] for can in cans}
        matches = pd.Series(matches).map(list)

        # create a frame with descriptors large structure in one column, and substructure match
        # indices in the second column
        tmp_df = descs_df.to_frame('descs')
        tmp_df['matches'] = matches

        # fetch atom labels for this smarts using the first molecule
        sub_labels = [f"atom{i + 1}" for i in range(len(matches[0]))]

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

    if 'core' in presets:
        # run MCS if there is more than 1 molecule
        if len(rd_mols) > 1:
            try:
                core_smarts = rdFMCS.FindMCS(list(rd_mols.values())).smartsString
            except ValueError:
                core_smarts = ''
        else:  # otherwise use the entire smiles as smarts string
            core_smarts = Chem.MolToSmarts(list(rd_mols.values())[0])

        if core_smarts:
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
    return data


def descriptors_from_mol_df(mol_df, conf_option='max') -> dict:
    """Get and weight descriptors given a set of molecules and a conformer reweighting option. This function involves a large query from the DB

    :param mol_df: dataframe returned by the autoqchem.db_functions.db_select_molecules function
    :type mol_df: pd.DataFrame
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

    # fetch ids and weights for all conformers for all molecules (lightweight query_
    all_ids = mol_df['_ids'].sum()
    cursor = feats_coll.find({"_id": {"$in": all_ids}}, {'weight': 1, 'molecule_id': 1})
    id_df = pd.DataFrame(cursor).sort_values(['molecule_id', 'weight'], ascending=False)

    # filter ids based on the conf_option (sometimes we don't need to fetch all conformers descriptors form the DB

    # helper function to manage different conf_options
    def filter_ids(group):
        # this function assumes each group is reverse ordered by weight (highest weight comes first)
        assert abs(group.weight.sum() - 1.) < 1e-6

        if conf_option == 'max':
            _ids = [group['_id'].iloc[0]]
        elif conf_option == 'min':
            _ids = [group['_id'].iloc[-1]]
        elif conf_option == 'any':
            _ids = [group['_id'].sample(1).iloc[0]]
        else:
            _ids = group['_id'].tolist()

        return _ids

    filtered_ids = id_df.groupby('molecule_id').apply(filter_ids).sum()

    # query that fetches all decsriptors from the DB (heavyweight query)
    cursor = feats_coll.find({"_id": {"$in": filtered_ids}}, {'molecule_id': 1,
                                                              'descriptors': 1,
                                                              'atom_descriptors': 1,
                                                              'transitions': 1,
                                                              'weight': 1,
                                                              'labels': 1})
    record_df = pd.Series([_pandatize_record(record) for record in cursor]).to_frame('records')

    # merge in can using molecule id
    record_df['molecule_id'] = record_df['records'].map(lambda r: r['molecule_id'])
    record_df = pd.merge(record_df, mol_df[['molecule_id', 'can']], how='left', on='molecule_id')

    # reweight descriptor conformers based on the option
    # helper function to manage different conf_options
    def reweigh_desc(group):

        if conf_option in ['min', 'max', 'any']:
            return group['records'].iloc[0]
        else:
            record = {}
            records = group['records'].tolist()
            record.update({"labels": records[0]['labels']})

            keys_to_reweigh = ['descriptors', 'atom_descriptors', 'transitions']

            for key in keys_to_reweigh:
                if conf_option == 'boltzmann':
                    dfs = pd.concat(r[key] * r['weight'] for r in records)
                    record[key] = dfs.groupby(dfs.index, sort=False).sum()
                if conf_option in ['mean', 'std']:
                    dfs = pd.concat(r[key] for r in records)
                    if conf_option == 'mean':
                        record[key] = dfs.groupby(dfs.index, sort=False).mean()
                    elif conf_option == 'std':
                        record[key] = dfs.groupby(dfs.index, sort=False).std()
            return record

    return record_df.groupby('can').apply(reweigh_desc)


def _pandatize_record(record) -> dict:
    """Convert json structures to pandas structures for an individual
    db record of a single conformation.

    :param record: db record of a single conformation
    :return: dict
    """

    del record['descriptors']['stoichiometry']

    record['descriptors'] = pd.Series(record['descriptors']).astype(float)
    record['atom_descriptors'] = pd.DataFrame(record['atom_descriptors']).astype(float)
    record['transitions'] = pd.DataFrame(record['transitions']).astype(float)

    # block for vibrational modes (currently unused)
    # record['modes'] = pd.DataFrame(record['modes']).astype(float)
    # if record['mode_vectors'] is not None:
    #    record['mode_vectors'] = pd.DataFrame(record['mode_vectors'])
    #    record['mode_vectors']['atom_idx'] = list(range(len(record['labels']))) * 3 * record['modes'].shape[0]
    #    record['mode_vectors'] = record['mode_vectors'].set_index(['mode_number', 'axis', 'atom_idx']).unstack(
    #        ['mode_number', 'axis'])
    #    record['mode_vectors'] = record['mode_vectors'].droplevel(0, axis=1).astype(float)
    # else:
    #    record['mode_vectors'] = pd.DataFrame()

    return record


class InconsistentLabelsException(Exception):
    """Raised when a set of molecules is inconsistently labeled"""
    pass
