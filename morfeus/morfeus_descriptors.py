import os
import argparse
import pickle

import pandas as pd
import numpy as np

from rdkit import Chem

from morfeus.conformer import ConformerEnsemble
from morfeus import Sterimol, BuriedVolume, XTB


def compute(smiles, n_confs=None, solvent=None):
    """Calculates a conformer ensemble using rdkit, then optimizes
    each conformer geometry using GFN2-xTB and calculates their properties"""

    # create conformer ensemble
    ce = ConformerEnsemble.from_rdkit(smiles, n_confs=n_confs,
                                      n_threads=os.cpu_count() - 1)

    # optimize conformer geometries
    ce.optimize_qc_engine(program="xtb",
                          model={'method': 'GFN2-xTB', "solvent": solvent},
                          procedure="berny")
    ce.prune_rmsd()

    # compute energies of the single point calculations
    ce.sp_qc_engine(program='xtb', model={'method': 'GFN2-xTB', "solvent": solvent})

    # sort on energy and generate an rdkit molecule
    ce.sort()
    ce.generate_mol()

    # compute xtb for the conformers
    for conf in ce:
        conf.xtb = XTB(ce.elements, conf.coordinates, solvent=solvent, version='2')

    return ce


def get_descriptors(conf_ensemble, substructure=None, substructure_labels=None,
                    sterimol_pairs=[]):
    """Extract descriptorrs from the conformational ensemble"""

    def get_substructure_match(conf_ensemble, substructure):
        """Helper function to find substructure match"""

        sub = Chem.MolFromSmarts(substructure)
        matches = conf_ensemble.mol.GetSubstructMatches(sub)
        if len(matches) > 1:
            raise ValueError(f"Substructer {substructure} is matched more than once in the molecule.")
        elif len(matches) == 0:
            raise ValueError(f"Substructer {substructure} not found in the molecule.")
        else:
            match = [m + 1 for m in matches[0]]  # add one to conform with morfeus 1-indexing vs rdkit 0-indexing

        print(match)
        return match

    def make_substructure_labels(match, substructure_labels):
        """Helper function to assign substructre atom labels"""

        if substructure_labels is None:
            labels = [f"atom{i}" for i in range(len(match))]
        else:
            if len(substructure_labels) != len(match):
                raise ValueError(f"Length of labels ({len(substructure_labels)}) is different\
                than the lenght of the substructure {substructure} match ({len(match)})")
            else:
                labels = substructure_labels

        return labels

        # prep for the substructure

    if substructure is not None:
        match = get_substructure_match(conf_ensemble, substructure)
        labels = make_substructure_labels(match, substructure_labels)
    else:
        match = None

    # loop over the conformers and get properties
    for conf in conf_ensemble.conformers:

        # XTB
        x = conf.xtb

        # molecular features
        conf.properties['IP'] = x.get_ip(corrected=True)
        conf.properties['EA'] = x.get_ea()
        conf.properties['HOMO'] = x.get_homo()
        conf.properties['LUMO'] = x.get_lumo()

        dip = x.get_dipole()
        conf.properties['dipole'] = np.sqrt(dip.dot(dip))

        conf.properties['electro'] = x.get_global_descriptor("electrophilicity", corrected=True)
        conf.properties['nucleo'] = x.get_global_descriptor("nucleophilicity", corrected=True)

        if match is not None:
            # atomic features
            charges = x.get_charges()
            electro = x.get_fukui('electrophilicity')
            nucleo = x.get_fukui('nucleophilicity')
            for idx, label in zip(match, labels):
                # charge, electrophilicity, nucleophilicity
                conf.properties[f"{label}_charge"] = charges[idx]
                conf.properties[f"{label}_electro"] = electro[idx]
                conf.properties[f"{label}_nucleo"] = nucleo[idx]

                # buried volumes
                conf.properties[f"{label}_VBur"] = BuriedVolume(conf_ensemble.elements, conf.coordinates,
                                                                idx).fraction_buried_volume

            # Sterimols
            for pair in sterimol_pairs:
                i1, i2 = pair[0], pair[1]
                match1, match2 = match[pair[0]], match[pair[1]]
                label = f"{labels[i1]}{labels[i2]}"
                s = Sterimol(conf_ensemble.elements, conf.coordinates, match1, match2)

                conf.properties[f"{label}_L"] = s.L_value
                conf.properties[f"{label}_B1"] = s.B_1_value
                conf.properties[f"{label}_B5"] = s.B_5_value
                conf.properties[f"{label}_length"] = s.bond_length

    props = {}
    for key in conf_ensemble.get_properties().keys():
        props[f"{key}_Boltz"] = conf_ensemble.boltzmann_statistic(key)
        props[f"{key}_Emin"] = conf_ensemble.get_properties()[key][0]

    return pd.Series(props)


if __name__ == "__main__":
    """executable script"""

    parser = argparse.ArgumentParser(description='Compute Conformational Ensemble and its Features')
    parser.add_argument("smiles", type=str, help="SMILES string of the molecule")
    parser.add_argument("name", type=str, help="Name of te molecule")
    parser.add_argument("output_path", type=str, help="Storage output path")
    parser.add_argument("--n_confs", type=int, help="Optional number of conformers to initially generate with RDKit. "
                                                    "If not specified a default is used based on the number of"
                                                    " rotatable bonds (50 if n_rot <=7, 200 if n_rot <=12, "
                                                    "300 for larger n_rot)",
                        default=None)
    parser.add_argument("--solvent", type=str, help="XTB supported solvents and their names can be found at "
                                                    "https://xtb-docs.readthedocs.io/en/latest/gbsa.html",
                        default=None)
    parser.add_argument("--descriptors", action='store_true',
                        help="Compute molecular descriptors")
    parser.add_argument("--substructure", type=str,
                        help="Substructure atoms get their individual descriptors")
    parser.add_argument("--substructure_labels", type=str,
                        nargs='+', help='List of labels for substructre query atoms', default=None)


    def sterimol_pair(s):
        try:
            i1, i2 = map(int, s.split(","))
            return i1, i2
        except:
            raise argparse.ArgumentTypeError("Sterimol pairs must be index1, index2")


    parser.add_argument("--sterimol_pairs", type=sterimol_pair, nargs="+",
                        help='Pair of atoms from the substructure match to extract the sterimol parameters (0-indexed)')

    args = parser.parse_args()
    print(args.__dict__)

    # compute the conformational ensemble
    m = compute(args.smiles, n_confs=args.n_confs, solvent=args.solvent)
    print(m.get_energies())

    # save output
    with open(f"{args.output_path}/{args.name}.pkl", "wb") as f:
        pickle.dump(m, f)

    if args.descriptors:
        descs = get_descriptors(m, substructure=args.substructure,
                                substructure_labels=args.substructure_labels,
                                sterimol_pairs=args.sterimol_pairs)

        with open(f"{args.output_path}/{args.name}_descriptors.csv", "wb") as f:
            descs.to_frame(args.smiles).T.to_csv(f, index_label="smiles")
