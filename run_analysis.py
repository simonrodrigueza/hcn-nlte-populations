"""
Compute and plot NLTE/LTE mean level-population ratios of the HCN 14 um lines
in their emission regions, for a pair of ProDiMo models (one NLTE, one LTE).

Example:
    python run_analysis.py --nlte-model path/to/M_R70_100 --lte-model path/to/M_R70_10000

Outputs (in --outdir):
    mean_pop_emission_region.pdf/.png  ratio of every NLTE line in its own emission region
    mean_pop_rep.pdf/.png              ratio of every NLTE line in the emission region of
                                       the brightest LTE line
    mean_pop_ratios.csv                the numbers behind both plots
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import prodimopy.read as pread

from pop_ratios import (create_tree, load_pop_file, mean_population_emission_region,
                        mean_population_rep, select_brightest_line)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--nlte-model', required=True,
                   help='ProDiMo NLTE model directory (must contain pop_LTE_<species>.out '
                        'and pop_PRO_<species>.out)')
    p.add_argument('--lte-model', required=True,
                   help='ProDiMo LTE model directory (provides the representative line)')
    p.add_argument('--species', default='HCN_H', help='line species (default: HCN_H)')
    p.add_argument('--wl-min', type=float, default=13.7, help='min. wavelength in um (default: 13.7)')
    p.add_argument('--wl-max', type=float, default=14.1, help='max. wavelength in um (default: 14.1)')
    p.add_argument('--outdir', default='output', help='output directory (default: output)')
    p.add_argument('--no-show', action='store_true', help='do not open the plot windows')
    return p.parse_args()


def plot_ratios(wl, ratios, common, title, ylabel, legend_loc, highlight_wl=None):
    """Scatter of ratio vs wavelength; lines shared with the LTE model in red."""
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.set_title(title)
    ax.set_xlabel(r'$\lambda$ [$\mu$m]')
    ax.set_ylabel(ylabel)

    ax.scatter(wl[~common], ratios[~common], color='black', alpha=0.7, label='NLTE exclusive')
    ax.scatter(wl[common], ratios[common], color='crimson', label='Common LTE/NLTE')

    if highlight_wl is not None:
        match = np.where(wl == highlight_wl)[0]
        if match.size > 0:
            ax.scatter(wl[match], ratios[match], color='gold', edgecolor='black', marker='*',
                       s=200, zorder=5, label='Representative line')

    ax.axhline(y=1, color='gray', linestyle='--', linewidth=1)
    ax.legend(loc=legend_loc)
    ax.set_yscale('log')
    fig.tight_layout()
    return fig


def main():
    args = parse_args()
    species = args.species
    wlrange = (args.wl_min, args.wl_max)
    os.makedirs(args.outdir, exist_ok=True)

    # --- models and population files ---
    nonlte_dir = os.path.abspath(args.nlte_model)
    lte_dir = os.path.abspath(args.lte_model)
    model_nonlte = pread.read_prodimo(nonlte_dir, name='nonlte')
    model_lte = pread.read_prodimo(lte_dir, name='lte')

    pop_lte_data = load_pop_file(os.path.join(nonlte_dir, f'pop_LTE_{species}.out'))
    pop_pro_data = load_pop_file(os.path.join(nonlte_dir, f'pop_PRO_{species}.out'))
    # Both pop files share the same cell ordering, so one tree serves for both.
    tree = create_tree(pop_lte_data)

    # --- lines and LTE/NLTE overlap ---
    # Lines are matched by wavelength. In the HCN Q branch there are no lines sharing
    # the same wavelength; check this before using another band.
    lines_nonlte = model_nonlte.selectLineEstimates(species, wlrange=list(wlrange))
    lines_lte = model_lte.selectLineEstimates(species, wlrange=list(wlrange))
    wl_nonlte = np.array([line.wl for line in lines_nonlte])
    wl_lte_set = {line.wl for line in lines_lte}
    common = np.array([wl in wl_lte_set for wl in wl_nonlte])

    # --- 1) every line in its own emission region ---
    ratio_own = mean_population_emission_region(pop_lte_data, pop_pro_data, tree, model_nonlte,
                                                species=species, wlrange=wlrange)
    fig1 = plot_ratios(wl_nonlte, ratio_own, common,
                       'Mean population ratio for each region of emission (NLTE lines)',
                       r'$n_{NLTE} / n_{LTE}$', legend_loc='center left')
    print('First plot created')

    # --- 2) every line in the emission region of the brightest LTE line ---
    rep_line = select_brightest_line(model_lte, species=species, wlrange=wlrange)
    ratio_rep = mean_population_rep(pop_lte_data, pop_pro_data, tree, model_lte, model_nonlte,
                                    rep_line, species=species, wlrange=wlrange)
    fig2 = plot_ratios(wl_nonlte, ratio_rep, common,
                       f'Mean population ratio for brightest LTE line region ({rep_line.wl:.2f} um)',
                       r'$n_{NLTE} / n_{LTE}$', legend_loc='upper right',
                       highlight_wl=rep_line.wl)
    print('Second plot created')

    # --- save ---
    for fig, name in ((fig1, 'mean_pop_emission_region'), (fig2, 'mean_pop_rep')):
        for ext in ('pdf', 'png'):
            fig.savefig(os.path.join(args.outdir, f'{name}.{ext}'), dpi=200)

    pd.DataFrame({
        'wl_um': wl_nonlte,
        'jup': [line.jup for line in lines_nonlte],
        'jlow': [line.jlow for line in lines_nonlte],
        'in_lte_model': common,
        'ratio_own_region': ratio_own,
        'ratio_rep_region': ratio_rep,
    }).to_csv(os.path.join(args.outdir, 'mean_pop_ratios.csv'), index=False)
    print(f'Results saved in {os.path.abspath(args.outdir)}')

    if not args.no_show:
        plt.show()


if __name__ == '__main__':
    main()
