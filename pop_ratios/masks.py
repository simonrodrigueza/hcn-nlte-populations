"""
Mapping of ProDiMo emission-region masks onto the rows of a pop_*.out file.

The pop files use a different spatial grid (120x121 cells) than the model grid
returned by prodimopy, so ``getLineOriginMask`` cannot be applied to them directly.
Instead, a k-d tree built on the pop-file cell centres is queried with the
coordinates of the cells inside the mask, keeping only the matching rows.

In addition, our LTE model was run with an older ProDiMo version than the NLTE
model, so their spatial grids are not in phase (even with ``refine_zz`` switched
off). ``ltemask_to_nonltemask`` rebuilds the emission-region mask of an LTE line
directly on the NLTE grid. For two models from the same ProDiMo version this
step may be unnecessary.
"""

import numpy as np
from scipy.spatial import cKDTree

from .io import COL_RA, COL_RB, COL_ZA, COL_ZB

AU_TO_M = 1.495978707e11


def create_tree(pop_data):
    """
    Build a cKDTree on the (r, z) centre of every cell of a pop file, computed
    from its borders ra, rb, za, zb (in metres).
    """
    r_center = 0.5 * (pop_data[:, COL_RA] + pop_data[:, COL_RB])
    z_center = 0.5 * (pop_data[:, COL_ZA] + pop_data[:, COL_ZB])
    return cKDTree(np.column_stack((r_center, z_center)))


def apply_mask_to_pop(pop_data, tree, model, mask, tol_m=None):
    """
    Return the rows of ``pop_data`` that correspond to the cells inside ``mask``.

    ``mask`` follows the prodimopy convention (0 = inside the emission region)
    on the grid of ``model``. Model coordinates (au) are converted to metres and
    the nearest pop-file cell is looked up in ``tree``. If ``tol_m`` is given,
    a warning is printed for matches farther than that distance.
    """
    radius = model.x[mask == 0] * AU_TO_M
    heights = model.z[mask == 0] * AU_TO_M
    obj_coords = np.column_stack((radius, heights))

    distances, indexes = tree.query(obj_coords, k=1, workers=-1)

    if tol_m is not None and np.any(distances > tol_m):
        bad_n = np.sum(distances > tol_m)
        print(f"{bad_n} cells without a close match (max. distance = {distances.max():.3e} m)")

    return pop_data[indexes, :]


def ltemask_to_nonltemask(model_lte, line_lte, model_nonlte):
    """
    Emission-region mask of ``line_lte`` (a LineEstimate of ``model_lte``),
    expressed on the spatial grid of ``model_nonlte``.

    The first part reproduces the prodimopy ``getLineOriginMask`` algorithm
    (radial range r15-r85 from the cumulative column flux, vertical range
    z15-z85 per column). That method only returns the final mask, so it is
    rebuilt here; the last part interpolates z15/z85 onto the NLTE radii and
    evaluates the mask on the NLTE grid.
    """
    # --- same as prodimopy getLineOriginMask, on the LTE grid ---
    r_lte = model_lte.x[:, 0]
    rinfo_index = {info.ix: info for info in line_lte.rInfo}

    fcolumn_lte = np.array([rinfo_index[ix].Fcolumn if ix in rinfo_index else 0.0
                            for ix in range(model_lte.nx)])
    z15_lte = np.array([rinfo_index[ix].z15 if ix in rinfo_index else 0.0
                        for ix in range(model_lte.nx)])
    z85_lte = np.array([rinfo_index[ix].z85 if ix in rinfo_index else 0.0
                        for ix in range(model_lte.nx)])

    cum_flux = np.cumsum(fcolumn_lte)
    if cum_flux[-1] <= 0:
        raise ValueError(f"LTE line {line_lte.wl:.8f} um does not have a positive flux")
    cum_fraction = cum_flux / cum_flux[-1]
    unique_fraction, unique_index = np.unique(cum_fraction, return_index=True)
    r_to_interpolate = r_lte[unique_index]
    r15 = np.interp(0.15, unique_fraction, r_to_interpolate)
    r85 = np.interp(0.85, unique_fraction, r_to_interpolate)

    # --- interpolate onto the NLTE grid and build the mask there ---
    r_nonlte = model_nonlte.x[:, 0]
    z15_nonlte = np.interp(r_nonlte, r_lte, z15_lte)
    z85_nonlte = np.interp(r_nonlte, r_lte, z85_lte)
    z_max_nonlte = np.maximum(z15_nonlte, z85_nonlte)
    z_min_nonlte = np.minimum(z15_nonlte, z85_nonlte)

    inside_radius = (r_nonlte >= r15) & (r_nonlte <= r85)
    inside_height = (
        (model_nonlte.z <= z_max_nonlte[:, np.newaxis])
        & (model_nonlte.z >= z_min_nonlte[:, np.newaxis])
    )
    included = inside_radius[:, np.newaxis] & inside_height

    return np.where(included, 0, 1).astype(np.uint8)
