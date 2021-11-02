from ipywidgets import interact, fixed
from rdkit import Chem
import py3Dmol


def graph_conf(m, confId=0):
    mb = Chem.MolToMolBlock(m, confId=confId)
    p = py3Dmol.view(width=400, height=400)
    p.removeAllModels()
    p.addModel(mb, 'sdf')
    p.setStyle({'stick': {}})
    p.setBackgroundColor('0xeeeeee')
    p.zoomTo()
    return p.show()


def draw(mol: Chem.Mol, energies: list) -> interact:
    """Make a drawing of all conformers in 3d"""

    emin = min(energies)
    energies = energies - emin
    p = py3Dmol.view(width=400, height=400)
    return interact(graph_conf, m=fixed(mol), p=fixed(p),
                    confId=(0, mol.GetNumConformers() - 1))
