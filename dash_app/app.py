import dash
import dash_bootstrap_components as dbc

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SANDSTONE],
                title='Auto-QChem DB',
                update_title="Loading...",
                meta_tags=[
                    {"name": "og image",
                     "property": "og:image",
                     "content": "assets/acq_thumb.png"}
                ])
server = app.server
app.config.suppress_callback_exceptions = True
