import os
import time
import re
from urllib.parse import unquote_plus

import dash
import dash_core_components as dcc
import dash_html_components as html
import flask
import pandas as pd
from dash.dependencies import Input, Output, State

from autoqchem.db_functions import Chem, descriptors, InconsistentLabelsException, db_connect
from dash_app.app import app, server
from dash_app.functions import app_path, get_table, get_tags_dropdown
from dash_app.layouts import layout_table, layout_descriptors, layout_navbar

app.layout = html.Div([layout_navbar(),
                       dcc.Location(id='url', refresh=False),
                       html.Div(id='page-content')
                       ])


@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname'),
               Input('url', 'search')])
def display_page(pathname, search):
    if pathname is None:
        return None
    elif pathname == f"/":
        if not search:
            return layout_table(queried=False)
        if search:
            items = [item.split("=") for item in search.split('?')[1].split("&")]
            items_dict = {key: unquote_plus(value) for key, value in items}
            try:
                if items_dict['substructure'] != "":
                    assert Chem.MolFromSmarts(items_dict['substructure']) is not None
                else:
                    pass
                items_dict['queried'] = True
                return layout_table(**items_dict)
            except AssertionError as e:
                return layout_table(None, None, message=f"Substructure '{items_dict['substructure']}'"
                                                        f" is an invalid SMARTS pattern.")
    elif pathname.startswith(f"/detail/"):
        id = pathname.split('/')[-1]
        return layout_descriptors(id)
    else:
        return '404'


@app.callback(Output('tags_dropdown', 'options'),
              [Input('basis_sets_dropdown', 'value'),
               Input('functionals_dropdown', 'value'),
               Input('solvents_dropdown', 'value')],
              prevent_initial_call=False)
def update_tags_dropdown(basis_set, functional, solvent):
    ctx = dash.callback_context
    # determine the source of the callback
    if ctx.triggered:
        return get_tags_dropdown(basis_set, functional, solvent)
    else:
        return dash.no_update


@app.callback(Output('inputPresetOptions', 'value'),
              [Input('dropdownPresetOptions', 'value')])
def pass_value_preset(v):
    if v is not None:
        return ",".join(v)
    else:
        return ""


@app.callback(Output('inputConformerOptions', 'value'),
              [Input('dropdownConformerOptions', 'value')])
def pass_value_conformer(v):
    return v if v else ""


@app.server.route('/', methods=['POST'])
def on_post():
    # fetch query items from url
    try:
        # check if there are url items to make a query
        url_items = [item.split("=") for item in flask.request.url.split('?')[1].split("&")]
    except IndexError:
        return ('', 204)
    items_dict = {key: unquote_plus(value) for key, value in url_items}
    # fetch form items from form
    items_dict.update(flask.request.form)
    # items_dict['tags'] = list(filter(None, items_dict['tags'].split(",")))

    # remove unused files in the dir
    for f in os.listdir(f"{app_path}/static/user_desc/"):
        fpath = f"{app_path}/static/user_desc/{f}"
        try:
            # windows trick, if a file can be renamed, it's free for removal
            os.rename(fpath, fpath)
            os.remove(fpath)
        except OSError:
            # cannot remove file, oh well, let it stay there
            pass

    # get timestamp
    ts = str(time.time()).replace(".", "")

    if 'export' in items_dict:
        path = f"{app_path}/static/user_desc/summary_{ts}.xlsx"
        df = get_table(tag=items_dict['tag'], substructure=items_dict['substructure'],
                       solvent=items_dict['solvent'], functional=items_dict['functional'],
                       basis_set=items_dict['basis_set'])
        data = {'summary': df.drop(['image', 'detail', 'tag'], axis=1)}

    elif ('PresetOptions' in items_dict) and ('ConformerOptions' in items_dict):
        path = f"{app_path}/static/user_desc/descriptors_{ts}.xlsx"
        items_dict['PresetOptions'] = items_dict['PresetOptions'].split(",")
        # extract the descriptors (this can take long)
        try:
            data = descriptors(
                tags=[items_dict['tag']] if items_dict['tag'] != 'ALL' else [],
                presets=items_dict['PresetOptions'],
                conf_option=items_dict['ConformerOptions'],
                solvent=items_dict['solvent'],
                functional=items_dict['functional'],
                basis_set=items_dict['basis_set'],
                substructure=items_dict['substructure']
            )
        except InconsistentLabelsException as e:
            return ('Molecules in the set have inconsistent labels', 200)
    else:
        return ('', 204)

    with pd.ExcelWriter(path) as writer:
        if data:
            for key, df in data.items():
                df.to_excel(writer, sheet_name=key)
        else:
            pd.DataFrame().to_excel(writer, sheet_name="Sheet1")

    # serve file
    return flask.send_file(path, as_attachment=True)


@app.callback(
    Output("molecule3d", "modelData"),
    Input("conf-slider", "value"),
    State("store", "data"),
    prevent_initial_call=True,
)
def conf(value, data):
    elements = data['elements']
    coords = data['coords']
    conn = data['conn']
    N = len(elements)
    modelData = {'atoms': [{'name': e, 'element': e, 'positions': c, 'serial': i
                            } for i, (e, c) in enumerate(zip(elements, coords[value - 1]))],
                 'bonds': [{'atom1_index': i, 'atom2_index': j, 'bond_order': int(conn[i][j])}
                           for j in range(N) for i in range(N) if (i < j) and (conn[i][j] > 0)]}
    return modelData


@app.callback(
    Output("conf_energy", "children"),
    Input("conf-slider", "value"),
    State("store", "data")
)
def conf_energy(value, data):
    energies = data['energies']
    return f"{value}: G = {min(energies):.1f} + {energies[value - 1] - min(energies):.2f} kcal/mol"


if __name__ == '__main__':
    server.run(port=80)
    # app.run_server(debug=True, port=80)
