import dash_core_components as dcc
import dash_html_components as html
import dash_table as dt

from dash_app.functions import *


def layout_table(tags, substructure, message=""):
    tags_coll = db_connect('tags')

    return html.Div(children=[

        html.H3('Auto-QChem DB'),

        html.Datalist(id='conf_options',
                      children=[html.Option(label=desc, value=tag)
                                for desc, tag in zip(conf_options_long, conf_options)]),

        html.Form(id='query-form', style={'display': 'inline-block'}, children=[
            dcc.Dropdown(
                id='tags_dropdown',
                options=[dict(label=f"{tag} ({len(tags_coll.distinct('molecule_id', {'tag': tag}))} molecules)"
                              , value=tag)
                         for tag in list(tags_coll.distinct('tag'))],
                multi=True,
                style={"width": "300px", 'display': 'inline-block', 'verticalAlign': 'top'},
                placeholder="Select tags...",
                persistence=True,
            ),
            dcc.Input(name="tags", id="tags", style={'display': 'none'}),
            dcc.Input(name="substructure", id="substructure",
                      placeholder="Filter on SMARTS substructure...", style={"width": "300px",
                                                                             'verticalAlign': 'top',
                                                                             'display': 'inline-block'},
                      persistence=False,
                      value=substructure
                      ),
            html.Button('Query', id='submit_query-form', style={"width": "150px", 'display': 'inline-block',
                                                                "background": "#F0F8FF"}),
        ], ),

        html.Form(id='export-summary-form', style={'display': 'inline-block'}, method='post', children=[
            dcc.Input(name="export", id='export', style={'display': 'none'}),
            html.Button('Export', id='submit_export-summary-form', style={"width": "150px", "background": "#F0F8FF"})
        ]
                  ),
        html.P(message) if message else html.Div(),

        html.Div(children=[
            html.Details(children=[
                html.Summary("Download descriptors"),
                html.Form(
                    id='export-form', children=[
                        dcc.Dropdown(
                            id='dropdownPresetOptions',
                            options=[dict(label=lab, value=val)
                                     for lab, val in zip(desc_presets_long, desc_presets)],
                            multi=True,
                            style={"width": "300px", 'verticalAlign': 'top', },
                            placeholder="Select descriptor presets...",
                        ),
                        dcc.Input(name="PresetOptions", id="inputPresetOptions", style={'display': 'none'}),
                        dcc.Input(name="ConformerOptions", id="conformerOptions", list='conf_options',
                                  style={"width": "300px", 'verticalAlign': 'top'},
                                  placeholder="Select conformer option..."),
                        html.Br(),
                        html.Button('Download', id='submit_export-form',
                                    style={"width": "150px", "background": "#F0F8FF"})
                    ],
                    method='post',
                )
            ]),
            dt.DataTable(
                id='table',
                data=get_table(tags, substructure).to_dict('records'),
                columns=[dict(name=c, id=c, hideable=False,
                              presentation="markdown" if c in ['image', 'descriptors'] else "input") for c in
                         ['image', 'can', 'name', 'tags', 'theory', 'light_basis_set', 'heavy_basis_set',
                          'generic_basis_set', 'max_light_atomic_number', 'num_conformers',
                          'max_num_conformers', 'descriptors']],
                editable=False,
                page_size=20, page_action="native",
                sort_action="native",
                sort_mode="multi",
                filter_action="native",
            ),
        ]) if tags is not None else html.Div(),
    ])


def layout_descriptors(id):
    feats_coll = db_connect('qchem_descriptors')
    mols_coll = db_connect('molecules')

    can = mols_coll.find_one(ObjectId(id), {'can': 1})['can']
    feat_ids = list(feats_coll.find({'molecule_id': ObjectId(id)}, {'_id': 1}))
    feat_ids = [item['_id'] for item in feat_ids]

    desc = descriptors_from_list_of_ids(feat_ids, conf_option='boltzmann')
    d = desc['descriptors'].to_frame().T

    df_atom = desc['atom_descriptors']
    df_atom.insert(0, 'label', desc['labels'])
    df_atom.insert(0, 'atom_idx', range(df_atom.shape[0]))

    vib = desc['modes']
    trans = desc['transitions']

    return html.Div(children=[

        html.H3(f"Descriptors for {can}"),
        html.Img(src=image(can)),
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
