from autoqchem.db_functions import *
from autoqchem.molecule import *

app_path = "C:/Users/Andrzej/Software/github/auto-qchem/dash_app"
# app_path = "C:/Users/AndrzejZuranski/Dropbox/DataX_PU/github/auto-qchem/dash_app"

def image(can):
    hash_str = hashlib.md5(can.encode()).hexdigest()

    if not os.path.exists(f"{app_path}/static/{hash_str}.svg"):
        mol = input_to_OBMol(can, input_type="string", input_format="can")
        OBMol_to_file(mol, "svg", f"{app_path}/static/{hash_str}.svg")
    return f"/static/{hash_str}.svg"


def image_3d(can, geom):
    mol = input_to_OBMol(can, input_type="string", input_format="can")
    mol.AddHydrogens()

    # adjust geometry
    for atom in pybel.ob.OBMolAtomIter(mol):
        pos = geom.iloc[atom.GetIdx() - 1]
        atom.SetVector(pos.X, pos.Y, pos.Z)

    pm = pybel.Molecule(mol)
    pybel.ipython_3d = True
    return pm._repr_html_()


def get_table(tag, substructure):
    df = db_select_molecules(tag=tag, substructure=substructure)
    # make a link, as data pass db _id of any conformer (limitation on what to pass in the link)
    # the link action will fetch all conformers and reweight them
    df['descriptors'] = df['_ids'].map(lambda ids: f'''|[descriptors](/descriptors/{str(ids[0])})|\n|:----:|''')
    df['image'] = df.can.map(image).map(lambda path: f"![]({path})")
    return df.drop(['weights', '_ids'], axis=1)
