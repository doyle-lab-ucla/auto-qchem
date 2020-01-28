import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from dash_app.app import app
from dash_app.layouts import layout_table, layout_descriptors

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
])


@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/':
        return layout_table(None)
    elif pathname.startswith('/collection/'):
        tag = pathname.split('/')[-1]
        return layout_table(tag)
    elif pathname.startswith('/descriptors/'):
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
    app.run_server(debug=False, port=8001, host='0.0.0.0')
