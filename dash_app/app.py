import dash
import dash_bootstrap_components as dbc

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SANDSTONE],
                title='Auto-QChem DB',
                update_title="Loading...",
                meta_tags=[
                    {"name": "image",
                     "property": "og:image",
                     "content": "https://autoqchem.org/assets/acq_thumb.png"},
                    {"name": "description",
                     "property": "og:description",
                     "content": "Auto-QChem is an automated workflow for the generation, storage, and retrieval of Density Functional Theory calculations for organic molecules."},
                    {"name": "author",
                     "content": "Andrzej Żurański"
                     }
                ])
server = app.server
app.config.suppress_callback_exceptions = True

app.scripts.config.serve_locally = False
app.scripts.append_script({
    'external_url': 'https://autoqchem.org/assets/async_src.js'
})
app.scripts.append_script({
    'external_url': 'https://autoqchem.org/assets/gtag.js'
})
