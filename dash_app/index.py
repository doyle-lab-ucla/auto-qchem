from urllib.parse import unquote

import dash_core_components as dcc
import dash_html_components as html
import flask
import pandas as pd
from dash.dependencies import Input, Output

from autoqchem.db_functions import pybel
from dash_app.app import app, server
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
            items_dict = {key: unquote(value) for key, value in items}
            try:
                if items_dict['Substructure'] != "":
                    pybel.Smarts(items_dict['Substructure'])
                return layout_table(items_dict['Collection'], items_dict['Substructure'])
            except OSError as e:
                return layout_table(None, None, message=f"Substructure '{items_dict['Substructure']}'"
                                                        f" is an invalid SMARTS pattern.")
    elif pathname.startswith(f"/descriptors/"):
        id = pathname.split('/')[-1]
        return layout_descriptors(id)
    else:
        return '404'


@app.server.route('/', methods=['POST'])
def on_post():
    print(flask.request.url, flask.request.form)

    url_items = [item.split("=") for item in flask.request.url.split('?')[1].split("&")]
    items_dict = {key: unquote(value) for key, value in url_items}

    path = "static/user_desc/descriptors.xlsx"

    df = pd.DataFrame([1, 2, 3])
    df.to_excel(path)

    return flask.send_file(path, as_attachment=True)


if __name__ == '__main__':
    server.run(port=80)
    # app.run_server(debug=True, port=80)
