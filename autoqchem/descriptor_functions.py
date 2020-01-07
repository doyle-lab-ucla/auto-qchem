import logging

import numpy as np
import pandas as pd
from openbabel import pybel
from scipy.spatial.distance import cdist

logger = logging.getLogger(__name__)


def occupied_volume(geometry_df, atom_idx, r, mesh_density=30):
    """compute occupied volume"""

    # make sure mesh_density is not outrageous
    MAX_MESH_DENSITY = 100
    if mesh_density > MAX_MESH_DENSITY:
        mesh_density = MAX_MESH_DENSITY
        logger.warning(f"Mesh density {mesh_density} is larger than allowed "
                       f"max of {MAX_MESH_DENSITY}. Using {MAX_MESH_DENSITY} instead.")

    # compute Van der Waals radii for atoms, r
    atom_r = geometry_df['AN'].map(pybel.ob.GetVdwRad)

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
    distances_sq = pd.DataFrame(cdist(mesh, coords[mesh_overlap_indices], metric='sqeuclidean'))

    # mesh cells are occupied if their distances are less then Van der Waals radius
    occupancy = distances_sq < atom_r[mesh_overlap_indices] ** 2

    # mesh cells are occupied if they are occupied by at least 1 atom
    occupied = occupancy.any(axis=1)

    return occupied.sum() / mesh.shape[0]
