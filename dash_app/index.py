import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from dash_app.app import app, server
from dash_app.layouts import layout_table, layout_descriptors

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
])


@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == f"/":
        return layout_table(None)
    elif pathname.startswith(f"/collection/"):
        tag = pathname.split('/')[-1]
        return layout_table(tag)
    elif pathname.startswith(f"/descriptors/"):
        id = pathname.split('/')[-1]
        return layout_descriptors(id)
    else:
        return '404'


@app.callback(Output("hidden_div_for_redirect_callback", "children"),
              [Input('field-dropdown', 'value')])
def redirect(value):
    if value is not None:
        return dcc.Location(pathname=f"/collection/{value}", id="")

if __name__ == '__main__':
    server.run(port=80)
    #app.run_server(debug=False, port=80)