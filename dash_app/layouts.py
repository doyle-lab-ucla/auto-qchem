import dash_core_components as dcc
import dash_html_components as html
import dash_table as dt
from bson.objectid import ObjectId

from dash_app.functions import *

db = db_connect()


def layout_table(tag, substructure, message=""):
    return html.Div(children=[
        html.Div(id="hidden_div_for_redirect_callback", hidden=True),
        html.H3('Autoqchem DFT descriptors database'),

        html.Datalist(id='collections',
                      children=[html.Option(label=f"{len(db.distinct('can', {'metadata.tag': tag}))} molecules",
                                            value=tag)
                                for tag in list(db.distinct('metadata.tag'))]),

        html.Form(id='form', children=[
            dcc.Input(name="Collection", id="collection", placeholder="Choose molecule collection...",
                      list='collections', style={"width": "300px"}, persistence=True),
            dcc.Input(name="Substructure", id="substructure",
                      placeholder="Filter on substructure, e.g. CCO...", style={"width": "300px"}, persistence=True),
            html.Button('Submit', id='submit_button')]),
        html.P(message) if message else html.Div(),
        dt.DataTable(
            id='table',
            data=get_table(tag, substructure).to_dict('records') if tag is not None else [],
            columns=[dict(name="image", id="image", hideable=True, presentation="markdown"),
                     dict(name="can", id="can", hideable=True),
                     dict(name='DFT functional', id="DFT_functional", hideable=True),
                     dict(name='DFT basis set', id="DFT_basis_set", hideable=True),
                     dict(name="num_conformers", id="num_conformers", hideable=True),
                     dict(name="max_num_conformers", id="max_num_conformers", hideable=True),
                     dict(name="descriptors", id="descriptors", hideable=True, presentation="markdown")
                     ],
            dropdown={'descriptors': {
                'options': [{'label': value, 'value': name} for name, value in zip(list('abc'), ['10', '20', '30'])]}},
            editable=False,
            hidden_columns=['max_num_conformers'],
            page_size=10, page_action="native",
            sort_action="native",
            sort_mode="multi",
            filter_action="native",
        ) if tag is not None else html.Div(),
    ])


def layout_descriptors(id):
    id_info = db.find_one(ObjectId(id), {'can': 1, 'metadata.tag': 1})
    desc = descriptors_from_can(id_info['can'], id_info['metadata']['tag'], conf_option='max')

    d = desc['descriptors'].to_frame().T

    df_atom = desc['atom_descriptors']
    df_atom.insert(0, 'label', desc['labels'])
    df_atom.insert(0, 'atom_idx', range(int(desc['descriptors']['number_of_atoms'])))

    vib = desc['modes']
    trans = desc['transitions']

    return html.Div(children=[

        html.H3(f"Descriptors for {id_info['can']}"),
        html.Img(src=image(id_info['can'])),
        dt.DataTable(data=d.to_dict('records'),
                     columns=[{'name': c, 'id': c, 'hideable': True} for c in d.columns]),
        html.Label("Atom-Level Descriptors:", style={"font-weight": "bold"}),
        dt.DataTable(data=df_atom.to_dict('records'),
                     columns=[{'name': c, 'id': c, 'hideable': True} for c in df_atom.columns],
                     sort_action='native',
                     filter_action='native'),
        html.Label("Excited State Transitions:", style={'font-weight': 'bold'}),
        dt.DataTable(data=trans.to_dict('records'),
                     columns=[{'name': c, 'id': c, 'hideable': True} for c in trans.columns],
                     sort_action='native',
                     filter_action='native'
                     ),
        html.Label("Vibrations:", style={'font-weight': 'bold'}),
        dt.DataTable(data=vib.to_dict('records'),
                     columns=[{'name': c, 'id': c, 'hideable': True} for c in vib.columns],
                     sort_action='native',
                     filter_action='native'
                     ),

    ])
