import hashlib

from autoqchem.db_functions import *
from autoqchem.molecule import *

app_path = "/home/ubuntu/github/auto-qchem/dash_app"


# app_path = "C:/Users/AndrzejZuranski/Dropbox/DataX_PU/github/auto-qchem/dash_app"


def image(can):
    hash_str = hashlib.md5(can.encode()).hexdigest()
    conv = pybel.ob.OBConversion()

    if not os.path.exists(f"{app_path}/static/{hash_str}.svg"):
        mol = pybel.ob.OBMol()
        conv.SetInFormat("can")
        conv.ReadString(mol, can)
        conv.SetOutFormat("svg")
        conv.WriteFile(mol, f"{app_path}/static/{hash_str}.svg")
    return f"/static/{hash_str}.svg"


def get_table(tag, substructure, solvent, functional, basis_set):
    df = db_select_molecules(tags=[tag] if tag != 'ALL' else [],
                             substructure=substructure, solvent=solvent, functional=functional, basis_set=basis_set)
    if df.empty:
        return df
    df['image'] = df.can.map(image).map(lambda path: f"![]({path})")
    df['detail'] = df['molecule_id'].map(lambda i: f'''|[detail](/detail/{str(i)})|\n|:----:|''')
    df['tags'] = df['tag'].map(lambda t: t.__repr__()[1:-1])

    df_metadata = pd.DataFrame(list(df.metadata))
    for c in ['class', 'subclass', 'type', 'subtype']:
        if c not in df_metadata.columns:
            df_metadata[c] = ''
    df_gaussian_config = pd.DataFrame(list(df_metadata.gaussian_config))
    df = pd.concat([df, df_gaussian_config, df_metadata[['class', 'subclass', 'type', 'subtype',
                                                         'max_num_conformers']]], axis=1)
    df['num_conf/max_conf'] = df.apply(lambda r: f"{r['num_conformers']}/{r['max_num_conformers']}", axis=1)
    df = df.drop(['num_conformers', 'max_num_conformers'], axis=1)
    df = df.loc[:, ~df.columns.duplicated()]  # deduplicate columns

    return df.drop(['molecule_id', 'metadata', '_ids'], axis=1)


def get_tags_dropdown(basis_set, functional, solvent):
    # connect to db
    filter = {}
    mols_coll = db_connect('molecules')
    tags_coll = db_connect('tags')

    if basis_set != 'ALL':
        filter['metadata.gaussian_config.light_basis_set'] = re.compile(f"^{re.escape(basis_set)}$", re.IGNORECASE)
    if functional != 'ALL':
        filter['metadata.gaussian_config.theory'] = re.compile(f"^{re.escape(functional)}$", re.IGNORECASE)
    if solvent != 'ALL':
        filter['metadata.gaussian_config.solvent'] = re.compile(f"^{re.escape(solvent)}$", re.IGNORECASE)

    if not filter:
        mol_ids = mols_coll.distinct('_id')
        available_tags = tags_coll.distinct('tag')
    else:
        mol_ids = mols_coll.distinct('_id', filter)
        available_tags = tags_coll.distinct('tag', {'molecule_id': {"$in": mol_ids}})

    options_tags = [dict(label=f'''All ({len(mol_ids)} molecules)''', value='ALL'),
                    dict(label="-------------------------------", value="", disabled="disabled")] + \
                   [dict(
                       label=f'''{tag} ({len(tags_coll.distinct("molecule_id", {"tag": tag, 'molecule_id': {'$in': mol_ids}}))} molecules)''',
                       value=tag
                   ) for tag in available_tags]

    return options_tags
