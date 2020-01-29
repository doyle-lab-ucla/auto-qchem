import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from autoqchem.db_functions import Chem
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
    print(pathname)

    if search:
        items = [item.split("=") for item in search.split('?')[1].split("&")]
        items_dict = {key: value for key, value in items}
        print(items_dict)
    if pathname == f"/":
        if not search:
            return layout_table(None, None)
        if search:
            if items_dict['Substructure']:
                if Chem.MolFromSmiles(items_dict['Substructure']) is None:
                    return layout_table(None, None, message=f"Substrucre '{items_dict['Substructure']}'"
                                                            f" is an invalid SMILES string.")
            return layout_table(items_dict['Collection'], items_dict['Substructure'])
    elif pathname.startswith(f"/?"):
        print("I'm here")
        items = pathname.split("?")[1].split("&")
        items = [item.split("=") for item in items]
        items_dict = {key: value for key, value in items}
        print(items_dict)
        tag = pathname.split('/')[-1]
        return layout_table(tag)
    elif pathname.startswith(f"/descriptors/"):
        id = pathname.split('/')[-1]
        return layout_descriptors(id)
    else:
        return '404'

if __name__ == '__main__':
    server.run(port=80)
    # app.run_server(debug=True, port=80)
