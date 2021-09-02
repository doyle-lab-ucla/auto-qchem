import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_table as dt

from dash_app.functions import *


def layout_navbar():
    """navigation bar"""

    navbar = dbc.Navbar(
        [
            html.A(
                dbc.Row(
                    [
                        dbc.Col(html.Img(src="/static/logo.png", height="40px")),
                        dbc.Col(dbc.NavbarBrand("Auto-QChem DB", className='ml-2'))
                    ],
                    no_gutters=True,
                    align='center'
                ),
                href="/"),
        ],
        color='primary',
        dark=True
    )
    return navbar


def layout_table(cls, subcls, type, subtype, tags, substructure, message=""):
    """main layout with a table of molecules"""

    tags_coll = db_connect('tags')
    mols_coll = db_connect('molecules')

    query_form = html.Form(id='query-form',
                           children=[
                               dbc.Table(
                                   html.Tbody(
                                       [
                                           html.Tr([
                                               html.Td(
                                                   dcc.Dropdown(id="cls_dropdown",
                                                                options=[dict(label=cls, value=cls)
                                                                         for cls in list(
                                                                        mols_coll.distinct('metadata.class')
                                                                    )],
                                                                placeholder="Select Class",
                                                                multi=False,
                                                                persistence=True,
                                                                )),
                                               html.Td(
                                                   dcc.Dropdown(id="subcls_dropdown",
                                                                options=[dict(label=subcls, value=subcls)
                                                                         for subcls in list(
                                                                        mols_coll.distinct('metadata.subclass')
                                                                    )],
                                                                placeholder="Select SubClass",
                                                                multi=False,
                                                                persistence=True,
                                                                ))
                                           ]),
                                           html.Tr([
                                               html.Td(
                                                   dcc.Dropdown(id="type_dropdown",
                                                                options=[dict(label=type, value=type)
                                                                         for type in list(
                                                                        mols_coll.distinct('metadata.type')
                                                                    )],
                                                                placeholder="Select Type",
                                                                multi=False,
                                                                persistence=True)),
                                               html.Td(
                                                   dcc.Dropdown(id="subtype_dropdown",
                                                                options=[dict(label=subtype, value=subtype)
                                                                         for subtype in list(
                                                                        mols_coll.distinct('metadata.subtype')
                                                                    )],
                                                                placeholder="Select SubType",
                                                                multi=False,
                                                                persistence=True)),
                                           ]),
                                           html.Tr(
                                               html.Td(
                                                   dcc.Dropdown(
                                                       id='tags_dropdown',
                                                       options=[dict(
                                                           label=f'''{tag} ({len(tags_coll.distinct("molecule_id", {"tag": tag}))} molecules)''',
                                                           value=tag)
                                                           for tag in tags_coll.distinct('tag')],
                                                       multi=True,
                                                       placeholder="Select Tags",
                                                       persistence=True,
                                                   ), colSpan=2)
                                           ),
                                           html.Tr(
                                               html.Td(
                                                   dbc.Input(name="substructure", id="substructure",
                                                             placeholder="SMARTS Substructure (optional)",
                                                             style={"width": "100%"},
                                                             persistence=True,
                                                             value=substructure
                                                             ), colSpan=2)
                                           ),
                                           html.Tr(
                                               html.Td(
                                                   dbc.Button('Query', id='submit_query-form',
                                                              color="primary",
                                                              block=True
                                                              ))
                                           )
                                       ],
                                   ),  # html.Tbody
                                   borderless=True,

                               ),
                               dbc.Input(name="cls", id="cls", style={'display': 'none'}),
                               dbc.Input(name="subcls", id="subcls", style={'display': 'none'}),
                               dbc.Input(name="type", id="type", style={'display': 'none'}),
                               dbc.Input(name="subtype", id="subtype", style={'display': 'none'}),
                               dbc.Input(name="tags", id="tags", style={'display': 'none'}, required=True),
                           ])

    export_summary_form = html.Form(id='export-summary-form',
                                    method='post', children=[
            dbc.Table(html.Tbody(html.Tr(html.Td(
                dbc.Button('Export List', id='submit_export-summary-form', color="primary", block=True)))),
                borderless=True),
            dbc.Input(name="export", id='export', style={'display': 'none'}),
        ])

    download_descriptors_form = html.Form(
        id='export-form', children=[
            dbc.Table(html.Tbody([
                html.Tr(html.Td(
                    dcc.Dropdown(
                        id='dropdownPresetOptions',
                        options=[dict(label=lab, value=val)
                                 for lab, val in zip(desc_presets_long, desc_presets)],
                        multi=True,
                        placeholder="Select descriptor presets...",
                    ))),
                html.Tr(html.Td(
                    dcc.Dropdown(
                        id='dropdownConformerOptions',
                        options=[dict(label=desc, value=tag)
                                 for desc, tag in zip(conf_options_long, conf_options)],
                        multi=False,
                        placeholder="Select conformer option..."
                    )
                )),
                html.Tr(html.Td(
                    dbc.Button('Download', id='submit_export-form', color="primary", block=True)
                ))]),
                borderless=True),
            dcc.Input(name="PresetOptions", id="inputPresetOptions", style={'display': 'none'}),
            dcc.Input(name="ConformerOptions", id="inputConformerOptions", style={'display': "none"})
        ],
        method='post',
    )

    queried = any(filter is not None for filter in (cls, subcls, type, subtype, tags))
    if queried:
        cols = ['image', 'can', 'tags', 'theory', 'light_basis_set', 'heavy_basis_set', 'generic_basis_set',
                'max_num_conformers', 'num_conformers', 'descriptors']
        df = get_table(cls, subcls, type, subtype, tags, substructure)[cols]

    content = [
        dbc.Tabs([
                     dbc.Tab(query_form, label="Query")] + ([
                                                                dbc.Tab(download_descriptors_form,
                                                                        label="Download Descriptors"),
                                                                dbc.Tab(export_summary_form,
                                                                        label="Export Molecule List")] if queried else []),
                 style={'margin-top': "5px", "margin-left": "10px"}
                 ),

        # placeholder for any error message
        html.P(message) if message else html.Div(),

        # molecule table
        html.Div(dt.DataTable(
            id="mol_table",
            data=df.to_dict(orient='records'),
            columns=[{'id': c, 'name': c, 'presentation': ('markdown' if c in ['image', 'descriptors'] else 'input')}
                     for c in df.columns],
            page_size=10,
            editable=False,

            style_table={
                "fontFamily": '-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Helvetica Neue",Arial,"Noto Sans",sans-serif,"Apple Color Emoji","Segoe UI Emoji","Segoe UI Symbol","Noto Color Emoji"'},
            style_header={
                'backgroundColor': 'primary',
                'fontWeight': 'bold',
                'padding': '0.75rem'
            },
            style_cell={
                "fontFamily": '-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Helvetica Neue",Arial,"Noto Sans",sans-serif,"Apple Color Emoji","Segoe UI Emoji","Segoe UI Symbol","Noto Color Emoji"',
                'fontWeight': '400',
                'lineHeight': '1.5',
                'color': '#212529',
                'textAlign': 'left',
                'whiteSpace': 'normal',
                'height': 'auto',
                'padding': '0.75rem',
                'border': '1px solid #dee2e6',
                'verticalAlign': 'top',
            },
            style_as_list_view=True,

        ),
            style={"margin": "10px"}) if queried else html.Div()
        # dbc.Table.from_dataframe(df, responsive=True) if queried else html.Div(),
    ]

    content_div = html.Div(children=content, style={"align": "center"})

    return content_div


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

        html.Img(src=image(can)),
        # html.Br(),
        html.Div([
            dbc.Label(f"{can}", style={"font-weight": "bold"}, size='lg'),
            dbc.Table.from_dataframe(d, responsive=True),
            dbc.Label("Atom-Level Descriptors", style={"font-weight": "bold"}, size='lg'),
            dbc.Table.from_dataframe(df_atom, responsive=True),
            dbc.Label("Excited State Transitions", style={'font-weight': 'bold'}, size='lg'),
            dbc.Table.from_dataframe(trans, responsive=True),
            dbc.Label("Vibrations:", style={'font-weight': 'bold'}, size='lg'),
            dbc.Table.from_dataframe(vib, responsive=True),
        ],
            style={'margin': '10px'}
        )
    ])
