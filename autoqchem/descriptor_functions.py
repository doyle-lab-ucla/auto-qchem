import logging

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from autoqchem.molecule import GetVdwRad

logger = logging.getLogger(__name__)


def occupied_volume(geometry_df, atom_idx, r, mesh_density=30) -> float:
    """Compute occupied volume fraction within a sphere of radius 'r' for an atom at position 'atom_idx'. Each atom \
    radius is taken to be its Van der Waals radius.

    :param geometry_df: geometry dataframe, must contain 'X', 'Y', 'Z' and 'AN' (atomic number) columns
    :type geometry_df: pd.DataFrame
    :param atom_idx: index of the atom to use as 'central' atom
    :type atom_idx: int
    :param r: occupied volume radius in Angstroms
    :type r: float
    :param mesh_density: density of the mesh for numerical integration (MAX=100)
    :type mesh_density: int
    :return: float, occupied volume fraction
    """

    # make sure mesh_density is not outrageous
    max_mesh_density = 100
    if mesh_density > max_mesh_density:
        mesh_density = max_mesh_density
        logger.warning(f"Mesh density {mesh_density} is larger than allowed "
                       f"max of {max_mesh_density}. Using {max_mesh_density} instead.")

    # fetch Van der Waals radii for atoms, r
    atom_r = geometry_df['AN'].map(GetVdwRad)

    # isolate coordinates
    coords = geometry_df[list('XYZ')]

    # create cubic mesh, then make it spherical, then move it into central atom
    ticks = np.linspace(-r, r, mesh_density)
    x, y, z = np.meshgrid(ticks, ticks, ticks)
    mesh = np.vstack((x.ravel(), y.ravel(), z.ravel())).T
    mesh = mesh[cdist(mesh, np.array([[0., 0., 0.]]), metric='sqeuclidean').ravel() < r ** 2]
    mesh = mesh + coords.iloc[atom_idx].values

    # filter atoms that are certainly not in the mesh, d > R + r
    atom_distances = cdist(coords.iloc[[atom_idx]], coords)[0]
    mesh_overlap_indices = (atom_distances - atom_r) < r

    # compute distance of every atom to every point in the mesh (this is the longest operation)
    distances_sq = pd.DataFrame(cdist(coords[mesh_overlap_indices], mesh, metric='sqeuclidean'),
                                index=atom_r[mesh_overlap_indices].index)
    # mesh cells are occupied if their distances are less then Van der Waals radius
    # the below comparison requires matching indexes in the distances_sq matrix and atom_r series
    occupancy = distances_sq.lt(atom_r[mesh_overlap_indices] ** 2, axis=0)

    # mesh cells are occupied if they are occupied by at least 1 atom
    occupied = occupancy.any()

    return occupied.sum() / mesh.shape[0]
