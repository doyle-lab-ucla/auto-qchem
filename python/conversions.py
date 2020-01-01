import logging

import pandas as pd
from openbabel import pybel
import rdkit.Chem

from python.enumerators import *

logger = logging.getLogger(__name__)
conv = pybel.ob.OBConversion()

def input_to_OBMol(input, input_type, input_format):
    """create OBMol from input (type and format as supported in the enumerators)"""

    mol = pybel.ob.OBMol()
    conv.SetInFormat(input_format.value)

    if input_type == input_types.file:
        conv.ReadFile(mol, input)
    elif input_type == input_types.string:
        conv.ReadString(mol, input)

    return mol

def OBMol_to_string(mol, format):
    """from OBMol to a string format"""

    conv.SetOutFormat(format)
    return conv.WriteString(mol).strip()


def OBMol_to_geom_df(mol):
    """extract geometry dataframe from OBMol"""

    pt = rdkit.Chem.GetPeriodicTable()
    array = [[pt.GetElementSymbol(a.GetAtomicNum()),
              a.GetIsotope(), a.x(), a.y(), a.z()] for a in pybel.ob.OBMolAtomIter(mol)]
    df = pd.DataFrame(array, columns=['Atom', 'Isotope'] + list('XYZ'))
    df['Atom'] = df['Atom'] + df['Isotope'].astype(str).where(df['Isotope'] > 0, '')
    df = df.drop('Isotope', axis=1)
    return df



