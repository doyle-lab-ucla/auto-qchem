import os
import time
from urllib.parse import unquote_plus

import dash_core_components as dcc
import dash_html_components as html
import flask
import pandas as pd
from dash.dependencies import Input, Output

from autoqchem.db_functions import pybel, descriptors, InconsistentLabelsException
from dash_app.app import app
from dash_app.functions import app_path, get_table
from dash_app.layouts import layout_table, layout_descriptors

app.layout = html.Div([
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
                return layout_table(items_dict['tags'], items_dict['substructure'])
            except OSError as e:
                return layout_table(None, None, message=f"Substructure '{items_dict['substructure']}'"
                                                        f" is an invalid SMARTS pattern.")
    elif pathname.startswith(f"/descriptors/"):
        id = pathname.split('/')[-1]
        return layout_descriptors(id)
    else:
        return '404'


@app.callback(Output('inputPresetOptions', 'value'),
              [Input('dropdownPresetOptions', 'value')])
def pass_value1(v):
    if v is not None:
        return ",".join(v)
    else:
        return ""


@app.callback(Output('tags', 'value'),
              [Input('tags_dropdown', 'value')])
def pass_value2(v):
    if v:
        return ",".join(v)
    else:
        return ""


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
        print("I need to export this")
        path = f"{app_path}/static/user_desc/summary_{ts}.xlsx"
        df = get_table(tags=items_dict['tags'], substructure=items_dict['substructure'])
        data = {'summary': df.drop(['image', 'descriptors', 'tag'], axis=1)}

    elif ('PresetOptions' in items_dict) and ('ConformerOptions' in items_dict):
        path = f"{app_path}/static/user_desc/descriptors_{ts}.xlsx"
        items_dict['PresetOptions'] = items_dict['PresetOptions'].split(",")
        # extract the descriptors (this can take long)
        try:
            data = descriptors(items_dict['tags'],
                               items_dict['PresetOptions'],
                               items_dict['ConformerOptions'],
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
