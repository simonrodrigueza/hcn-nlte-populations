"""
Mean NLTE/LTE level-population ratios over line emission regions.

For a transition u -> l, the mean excitation inside a region is the ratio of the
total number of molecules in the upper and lower levels, summing over cells

    N_j = sum_cells  x_j * n_mol * V_cell ,   V_cell = pi (rb^2 - ra^2) (zb - za)

so the result is weighted by volume and molecule density. The quantity returned
for every line is (N_u / N_l)_NLTE / (N_u / N_l)_LTE, where NLTE comes from
pop_PRO_<species>.out and LTE from pop_LTE_<species>.out of the same NLTE model.

Two variants are provided:
  * mean_population_emission_region: each line is evaluated in its own emission
    region, so different lines may probe different parts of the disk
    (e.g. inner vs outer disk).
  * mean_population_rep: all lines are evaluated in one fixed region, the
    emission region of a single representative line of an independent LTE
    model (by default its brightest line, see select_brightest_line).
"""

import numpy as np

from .io import COL_NMOL, COL_RA, COL_RB, COL_ZA, COL_ZB, level_column
from .masks import apply_mask_to_pop, ltemask_to_nonltemask

DEFAULT_SPECIES = 'HCN_H'
DEFAULT_WLRANGE = (13.7, 14.1)   # micron, HCN Q branch


def _cell_volumes(rows):
    """Volume of each (annular) cell, from its borders."""
    return np.pi * (rows[:, COL_RB]**2 - rows[:, COL_RA]**2) * (rows[:, COL_ZB] - rows[:, COL_ZA])


def _level_ratio(rows, vol, line):
    """Weighted N_u / N_l for ``line`` over ``rows`` (0 if the lower level is empty)."""
    total_u = np.sum(rows[:, level_column(line.jup)] * rows[:, COL_NMOL] * vol)
    total_l = np.sum(rows[:, level_column(line.jlow)] * rows[:, COL_NMOL] * vol)
    return total_u / total_l if total_l != 0 else 0


def _nlte_lte_ratio(rows_pro, rows_lte, vol, line):
    """(N_u/N_l)_NLTE / (N_u/N_l)_LTE, NaN if the LTE ratio vanishes."""
    mean_pop_pro = _level_ratio(rows_pro, vol, line)
    mean_pop_lte = _level_ratio(rows_lte, vol, line)
    return mean_pop_pro / mean_pop_lte if mean_pop_lte != 0 else np.nan


def mean_population_emission_region(pop_lte_data, pop_pro_data, tree, model_nonlte,
                                    species=DEFAULT_SPECIES, wlrange=DEFAULT_WLRANGE):
    """
    NLTE/LTE population ratio of every line of ``species`` in ``wlrange`` of the
    NLTE model, each one averaged over its own emission region.
    Returns an array ordered as model_nonlte.selectLineEstimates(species, wlrange).
    """
    lines_nonlte = model_nonlte.selectLineEstimates(species, wlrange=list(wlrange))
    mean_pop = []

    for line in lines_nonlte:
        mask_nonlte = model_nonlte.getLineOriginMask(line)
        rows_pop_pro = apply_mask_to_pop(pop_pro_data, tree, model_nonlte, mask_nonlte)
        rows_pop_lte = apply_mask_to_pop(pop_lte_data, tree, model_nonlte, mask_nonlte)
        vol = _cell_volumes(rows_pop_lte)
        mean_pop.append(_nlte_lte_ratio(rows_pop_pro, rows_pop_lte, vol, line))

    return np.array(mean_pop)


def select_brightest_line(model, species=DEFAULT_SPECIES, wlrange=DEFAULT_WLRANGE):
    """
    The LineEstimate with the highest flux among the lines of ``species`` in
    ``wlrange`` of ``model``. Meant to be used as the representative line
    (``line_lte``) of mean_population_rep.
    """
    lines = model.selectLineEstimates(species, wlrange=list(wlrange))
    return max(lines, key=lambda line: line.flux)


def mean_population_rep(pop_lte_data, pop_pro_data, tree, model_lte, model_nonlte, line_lte,
                        species=DEFAULT_SPECIES, wlrange=DEFAULT_WLRANGE):
    """
    NLTE/LTE population ratio of every line of ``species`` in ``wlrange`` of the
    NLTE model, all averaged over one fixed region: the emission region
    (r15-r85, z15-z85) of ``line_lte``, a LineEstimate of the LTE model chosen
    by the caller (e.g. with select_brightest_line).
    Returns an array ordered as model_nonlte.selectLineEstimates(species, wlrange).
    """
    lines_nonlte = model_nonlte.selectLineEstimates(species, wlrange=list(wlrange))

    mask_nonlte = ltemask_to_nonltemask(model_lte, line_lte, model_nonlte)
    rows_pop_pro = apply_mask_to_pop(pop_pro_data, tree, model_nonlte, mask_nonlte)
    rows_pop_lte = apply_mask_to_pop(pop_lte_data, tree, model_nonlte, mask_nonlte)
    vol = _cell_volumes(rows_pop_lte)

    return np.array([_nlte_lte_ratio(rows_pop_pro, rows_pop_lte, vol, line)
                     for line in lines_nonlte])
