import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_table as dt
import dash_bio as dbio

from dash_app.functions import *
from autoqchem.rdkit_utils import extract_from_rdmol

jmolcolors = pd.read_csv("assets/jmolcolors.csv").set_index('atom')['Hex']

def layout_navbar():
    """navigation bar"""

    navbar = dbc.Navbar(
        [
            html.A(
                dbc.Row(
                    [
                        dbc.Col(html.Img(src="/assets/logo.png", height="40px")),
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


def layout_table(tag='ALL', substructure=None, solvent='ALL', functional='ALL', basis_set='ALL', message="",
                 queried=False):
    """main layout with a table of molecules"""

    tags_coll = db_connect('tags')
    mols_coll = db_connect('molecules')

    solvents = list(set(map(str.upper, mols_coll.distinct('metadata.gaussian_config.solvent'))))
    functionals = list(set(map(str.upper, mols_coll.distinct('metadata.gaussian_config.theory'))))
    basis_sets = list(set(map(str.upper, mols_coll.distinct('metadata.gaussian_config.light_basis_set'))))

    options_tags = get_tags_dropdown(basis_set, functional, solvent)

    query_form = dbc.Form(id='query-form', prevent_default_on_submit=False,
                          children=[
                              dbc.Table(
                                  html.Tbody(
                                      [
                                          html.Tr([
                                              html.Td(dbc.Button("Dataset Tag", type='button', style={'width': '100%'}),
                                                      style={'width': '80px'}),
                                              html.Td(
                                                  dbc.Select(
                                                      name='tag',
                                                      id='tags_dropdown',
                                                      options=options_tags,
                                                      required=True,
                                                      valid='ALL',
                                                      persistence="session",
                                                  ))]),
                                          html.Tr([
                                              html.Td(dbc.Button("Solvent", type='Button', style={'width': '100%'}),
                                                      style={'width': '80px'}),
                                              html.Td(
                                                  dbc.Select(
                                                      name='solvent',
                                                      id='solvents_dropdown',
                                                      options=[dict(label=f'ALL', value='ALL'),
                                                               dict(label="-------------------------------", value="",
                                                                    disabled="disabled")] +
                                                              [dict(label=solvent, value=solvent) for solvent in
                                                               solvents],
                                                      required=True,
                                                      value='ALL',
                                                      persistence="session",
                                                  ))]),
                                          html.Tr([
                                              html.Td(dbc.Button("Functional", type='button', style={'width': '100%'}),
                                                      style={'width': '80px'}),
                                              html.Td(
                                                  dbc.Select(
                                                      name='functional',
                                                      id='functionals_dropdown',
                                                      options=[dict(label=f'ALL', value='ALL'),
                                                               dict(label="-------------------------------", value="",
                                                                    disabled="disabled")] +
                                                              [dict(label=functional, value=functional) for functional
                                                               in functionals],
                                                      required=True,
                                                      value='ALL',
                                                      persistence="session",
                                                  ))]),
                                          html.Tr([
                                              html.Td(dbc.Button("Basis Set", type='button', style={'width': '100%'}),
                                                      style={'width': '80px'}),
                                              html.Td(
                                                  dbc.Select(
                                                      name='basis_set',
                                                      id='basis_sets_dropdown',
                                                      options=[dict(label=f'ALL', value='ALL'),
                                                               dict(label="-------------------------------", value="",
                                                                    disabled="disabled")] +
                                                              [dict(label=basis_set, value=basis_set) for basis_set in
                                                               basis_sets],
                                                      required=True,
                                                      value='ALL',
                                                      persistence="session"
                                                  ))]),
                                          html.Tr([
                                              html.Td(
                                                  dbc.Button("Substructure", type='button', style={'width': '100%'}),
                                                  style={'width': '80px'}),
                                              html.Td(
                                                  dbc.Input(name="substructure",
                                                            id="substructure",
                                                            placeholder="SMARTS Substructure",
                                                            style={"width": "100%"},
                                                            persistence="session",
                                                            value=substructure,
                                                            ))]),
                                          html.Tr(
                                              html.Td(
                                                  dbc.Button('Query', id='submit_query-form',
                                                             color="primary",
                                                             block=True,
                                                             type='submit',
                                                             ), colSpan=2)),
                                      ],
                                   ),  # html.Tbody
                                   borderless=True,
                                  responsive=True
                               ),
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
                    dbc.Select(
                        id='dropdownConformerOptions',
                        options=[dict(label=desc, value=val)
                                 for desc, val in zip(conf_options_long, conf_options)],
                        required=True,
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

    if queried:
        cols = ['image', 'can', 'name', 'tags', 'solvent', 'theory', 'light_basis_set', 'heavy_basis_set',
                'num_conf/max_conf', 'detail']
        df = get_table(tag=tag, substructure=substructure, solvent=solvent, functional=functional, basis_set=basis_set)
        if not df.empty:
            df = df[cols]
        else:
            df = pd.DataFrame(columns=cols)

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
        html.Div([dbc.Label(f"Found {df.shape[0]} molecules"), dt.DataTable(
            id="mol_table",
            data=df.to_dict(orient='records'),
            columns=[{'id': c, 'name': c, 'presentation': ('markdown' if c in ['image', 'detail'] else 'input')}
                     for c in df.columns],
            page_size=30,
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

        )],
                 style={"margin": "10px"}) if queried else html.Div()
        # dbc.Table.from_dataframe(df, responsive=True) if queried else html.Div(),
    ]

    content_div = html.Div(children=content, style={"align": "center"})

    return content_div


def layout_descriptors(id):
    mols_df = db_select_molecules(molecule_ids=[ObjectId(id)])
    mol_record = db_connect('molecules').find_one({'_id': ObjectId(id)})

    try:
        rdmol, energies = db_get_rdkit_mol(mol_record)
        # make alternating single/double bonds for aromatics
        Chem.SanitizeMol(rdmol)
        Chem.Kekulize(rdmol, clearAromaticFlags=True)

        # extract geometry info from the molecule
        elements, coords, conn, _ = extract_from_rdmol(rdmol)
        N = len(elements)
        Nconf = len(coords)
        if Nconf <= 10:
            marks = {i: str(i) for i in range(1, Nconf + 1)}
        else:
            marks = {i: str(i) for i in range(1, Nconf + 1, int(Nconf / 10))}

        modelData = {'atoms': [{'name': e, 'element': e, 'positions': c.tolist(), 'serial': i
                                } for i, (e, c) in enumerate(zip(elements, coords[0]))],
                     'bonds': [{'atom1_index': i, 'atom2_index': j, 'bond_order': int(conn[i][j])}
                               for j in range(N) for i in range(N) if (i < j) and (conn[i][j] > 0)]}
        stylesData = {str(i): {"color": f"#{jmolcolors.loc[e]}",
                               "visualization_type": "stick"} for i, e in enumerate(elements)}
        ismol = True
    except AssertionError:
        logger.warning("Cannot create rdmol from DB.")
        ismol = False

    can = mols_df['can'].iloc[0]
    desc = descriptors_from_mol_df(mols_df, conf_option='boltzmann').iloc[0]

    d = desc['descriptors'].to_frame().T
    df_atom = desc['atom_descriptors']
    df_atom.insert(0, 'label', desc['labels'])
    df_atom.insert(0, 'atom_idx', range(df_atom.shape[0]))

    trans = desc['transitions']

    return html.Div(children=[
        dcc.Store(id='store', data={'elements': elements, 'coords': coords, 'conn': conn, 'energies': energies})
        if ismol else dcc.Store(id='store'),

        html.Div([
            dbc.Label(f"{can}", style={"font-weight": "bold"}, size='lg'),
            html.Div(children=[html.P(id='conf_energy', children="", style={"font-weight": "bold", 'align': 'center'}),
                               dbio.Molecule3dViewer(id='molecule3d', modelData=modelData, styles=stylesData),
                               dcc.Slider(id='conf-slider',
                                          min=1, max=Nconf, step=1, value=1,
                                          marks=marks)],
                     style={'width': '500px'}) if ismol else
            html.Div(children=[html.Br(), html.Img(src=image(can))]),
            dbc.Table.from_dataframe(d, responsive=True),
            dbc.Label("Atom-Level Descriptors", style={"font-weight": "bold"}, size='lg'),
            dbc.Table.from_dataframe(df_atom, responsive=True),
            dbc.Label("Excited State Transitions", style={'font-weight': 'bold'}, size='lg'),
            dbc.Table.from_dataframe(trans, responsive=True),
        ],
            style={'margin': '10px'}
        )
    ])
