import os
import time
from urllib.parse import unquote_plus

import dash_core_components as dcc
import dash_html_components as html
import flask
import pandas as pd
from dash.dependencies import Input, Output

from autoqchem.db_functions import pybel, descriptors, InconsistentLabelsException, db_connect
from dash_app.app import app, server
from dash_app.functions import app_path, get_table
from dash_app.layouts import layout_table, layout_descriptors, layout_navbar

app.layout = html.Div([layout_navbar(),
                       dcc.Location(id='url', refresh=False),
                       html.Div(id='page-content')
                       ])


@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname'),
               Input('url', 'search')])
def display_page(pathname, search):
    if pathname == f"/":
        if not search:
            return layout_table(None, None)
        if search:
            items = [item.split("=") for item in search.split('?')[1].split("&")]
            items_dict = {key: unquote_plus(value) for key, value in items}
            items_dict['tags'] = list(filter(None, items_dict['tags'].split(",")))
            try:
                if items_dict['substructure'] != "":
                    pybel.Smarts(items_dict['substructure'])
                else:
                    pass
                return layout_table(items_dict['tags'],
                                    items_dict['substructure'],
                                    )
            except OSError as e:
                return layout_table(None, None, message=f"Substructure '{items_dict['substructure']}'"
                                                        f" is an invalid SMARTS pattern.")
    elif pathname.startswith(f"/descriptors/"):
        id = pathname.split('/')[-1]
        return layout_descriptors(id)
    else:
        return '404'


@app.callback([
    Output('cls_dropdown', 'options'),
    Output('subcls_dropdown', 'options'),
    Output('type_dropdown', 'options'),
    Output('subtype_dropdown', 'options'),
    Output('tags_dropdown', 'options')],
    [Input('cls_dropdown', 'value'),
     Input('subcls_dropdown', 'value'),
     Input('type_dropdown', 'value'),
     Input('subtype_dropdown', 'value')])
def update_type_dropdown(cls, subcls, type, subtype):
    filter = {}
    if cls is not None:
        filter['metadata.class'] = cls
    if subcls is not None:
        filter['metadata.subclass'] = subcls
    if type is not None:
        filter['metadata.type'] = type
    if subtype is not None:
        filter['metadata.subtype'] = subtype

    mols_coll = db_connect('molecules')
    tags_coll = db_connect('tags')
    options_cls = [dict(label=type, value=type) for type in list(mols_coll.distinct('metadata.class', filter))]
    options_subcls = [dict(label=type, value=type) for type in list(mols_coll.distinct('metadata.subclass', filter))]
    options_type = [dict(label=type, value=type) for type in list(mols_coll.distinct('metadata.type', filter))]
    options_subtype = [dict(label=type, value=type) for type in list(mols_coll.distinct('metadata.subtype', filter))]

    if filter:
        mols_ids = mols_coll.distinct('_id', filter)
        available_tags = tags_coll.distinct('tag', {'molecule_id': {"$in": mols_ids}})
        options_tags = [dict(label=f'''{tag} ({len(tags_coll.distinct("molecule_id", {"tag": tag,
                                                                                      'molecule_id': {"$in": mols_ids}}))} molecules)''',
                             value=tag) for tag in available_tags]
    else:
        available_tags = tags_coll.distinct('tag')
        options_tags = [dict(label=f'''{tag} ({len(tags_coll.distinct("molecule_id", {"tag": tag}))} molecules)''',
                             value=tag) for tag in available_tags]

    return options_cls, options_subcls, options_type, options_subtype, options_tags


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


@app.callback(Output('tags', 'value'),
              [Input('tags_dropdown', 'value')])
def pass_value_tags(v):
    return "" if v == "All" else v


@app.callback(Output('cls', 'value'),
              [Input('cls_dropdown', 'value')])
def pass_value_cls(v):
    return v if v else ""


@app.callback(Output('subcls', 'value'),
              [Input('subcls_dropdown', 'value')])
def pass_value_subcls(v):
    return v if v else ""


@app.callback(Output('type', 'value'),
              [Input('type_dropdown', 'value')])
def pass_value_type(v):
    return v if v else ""


@app.callback(Output('subtype', 'value'),
              [Input('subtype_dropdown', 'value')])
def pass_value_subtype(v):
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
    items_dict['tags'] = list(filter(None, items_dict['tags'].split(",")))

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
        df = get_table(
            cls=None,
            subcls=None,
            type=None,
            subtype=None,
            tags=items_dict['tags'],
            substructure=items_dict['substructure'])
        data = {'summary': df.drop(['image', 'descriptors', 'tag'], axis=1)}

    elif ('PresetOptions' in items_dict) and ('ConformerOptions' in items_dict):
        path = f"{app_path}/static/user_desc/descriptors_{ts}.xlsx"
        items_dict['PresetOptions'] = items_dict['PresetOptions'].split(",")
        # extract the descriptors (this can take long)
        try:
            data = descriptors(
                cls=None,
                subcls=None,
                type=None,
                subtype=None,
                tags=items_dict['tags'],
                presets=items_dict['PresetOptions'],
                conf_option=items_dict['ConformerOptions'],
                substructure=items_dict['substructure'])
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


if __name__ == '__main__':
    server.run(port=80)
    # app.run_server(debug=True, port=80)
