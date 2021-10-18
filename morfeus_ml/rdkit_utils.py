import py3Dmol
from ipywidgets import interact, fixed

from rdkit import Chem


def draw_all(rdmol):
    def draw_one(m, p, confId=-1):
        mb = Chem.MolToMolBlock(m, confId=confId)
        p.removeAllModels()
        p.addModel(mb, 'sdf')
        p.setStyle({'stick': {}})
        p.setBackgroundColor('0xeeeeee')
        p.zoomTo()
        return p.show()

    p = py3Dmol.view(width=400, height=400)
    interact(draw_one, m=fixed(rdmol), p=fixed(p),
             confId=(0, rdmol.GetNumConformers() - 1))
