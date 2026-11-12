"""Reading of ProDiMo ``pop_*.out`` level-population files."""

import numpy as np
import pandas as pd

# Column layout of the pop_*.out files (see the "columns=" header line):
#   id, ra, rb, za, zb, rd, tk, nm, vr, vz, va, db, td, lp(1), lp(2), ...
# Coordinates are cell borders in metres, nm is the molecule number density and
# lp(j) is the fractional population of level j (1-based).
COL_RA, COL_RB, COL_ZA, COL_ZB = 1, 2, 3, 4
COL_NMOL = 7
COL_FIRST_LEVEL = 13


def level_column(j):
    """Column index of the fractional population of level ``j`` (1-based, as in LineEstimates)."""
    return COL_FIRST_LEVEL + (j - 1)


def load_pop_file(path):
    """
    Load a pop_*.out file (pop_LTE_<species>.out or pop_PRO_<species>.out) into a numpy array.

    The data block starts after the line containing only '@'. Files are large
    (~1 GB for HCN), so this can take a while.
    Returns a 2D float64 array with one row per spatial cell.
    """
    print('reading', path)
    initial_line = None
    with open(path, 'r') as f:
        for i, line in enumerate(f):
            if line.strip() == '@':
                initial_line = i + 1
                break
    if initial_line is None:
        raise ValueError(f"'@' separator not found in {path}; is this a ProDiMo pop file?")

    df = pd.read_csv(path, sep=r'\s+', skiprows=initial_line,
                     header=None, engine='c', dtype=np.float64)
    data = df.to_numpy()
    print('pop file loaded:', path, 'shape:', data.shape)
    return data
