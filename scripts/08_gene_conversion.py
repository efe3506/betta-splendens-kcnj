#!/usr/bin/env python3
"""Sliding-window test for gene conversion in the kcnj10 subfamily."""
import itertools
import os
import subprocess

from common import read_fasta, write_fasta

SRC = 'results/03_phylogeny/kcnj_labeled.faa'
OUT = 'results/04_structure'
FOCAL = 'Bs_LOC114844899'
TAXA = [FOCAL, 'Bs_kcnj10a', 'Dr_kcnj10a', 'Dr_kcnj10b', 'Hs_KCNJ10']
WINDOW, STEP, MIN_SITES = 60, 10, 10


def identity(a, b):
    pairs = [(x, y) for x, y in zip(a, b) if x != '-' and y != '-']
    if len(pairs) < MIN_SITES:
        return None
    return 100.0 * sum(x == y for x, y in pairs) / len(pairs)


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs('figures', exist_ok=True)
    seqs = dict(read_fasta(SRC))
    write_fasta(f'{OUT}/kcnj10_subset.faa', [(t, seqs[t]) for t in TAXA])
    with open(f'{OUT}/kcnj10_subset_aln.faa', 'w') as fh:
        subprocess.run(['mafft', '--localpair', '--maxiterate', '1000', '--quiet',
                        f'{OUT}/kcnj10_subset.faa'], stdout=fh, check=True)
    aln = dict(read_fasta(f'{OUT}/kcnj10_subset_aln.faa'))

    with open(f'{OUT}/kcnj10_pairwise_identity.tsv', 'w') as fh:
        fh.write('seq1\tseq2\tidentity_pct\n')
        for a, b in itertools.combinations(TAXA, 2):
            fh.write(f'{a}\t{b}\t{identity(aln[a], aln[b]):.1f}\n')

    rows = []
    length = len(aln[FOCAL])
    for start in range(0, length - WINDOW + 1, STEP):
        w = slice(start, start + WINDOW)
        vals = (identity(aln[FOCAL][w], aln['Bs_kcnj10a'][w]),
                identity(aln[FOCAL][w], aln['Dr_kcnj10b'][w]),
                identity(aln[FOCAL][w], aln['Dr_kcnj10a'][w]),
                identity(aln['Bs_kcnj10a'][w], aln['Dr_kcnj10a'][w]))
        if None not in vals:
            rows.append((start + 1,) + vals)

    with open(f'{OUT}/gene_conversion_windows.tsv', 'w') as fh:
        fh.write('window_start\tvs_Bs_kcnj10a\tvs_Dr_kcnj10b\tvs_Dr_kcnj10a\t'
                 'control_Bs_kcnj10a_vs_Dr_kcnj10a\n')
        for r in rows:
            fh.write(f'{r[0]}\t' + '\t'.join(f'{v:.1f}' for v in r[1:]) + '\n')

    closer_b = sum(1 for r in rows if r[2] > r[1])
    print(f'{len(rows)} windows; {closer_b} in which {FOCAL} is closer to Dr_kcnj10b '
          f'than to Bs_kcnj10a')

    import figstyle
    import matplotlib.pyplot as plt
    x = [r[0] for r in rows]
    delta = [r[1] - r[2] for r in rows]
    fig, ax = plt.subplots(2, 1, figsize=(9, 6), sharex=True,
                           gridspec_kw={'height_ratios': [2, 1]})
    ax[0].plot(x, [r[1] for r in rows], lw=2, label='vs Bs_kcnj10a')
    ax[0].plot(x, [r[2] for r in rows], lw=2, label='vs Dr_kcnj10b')
    ax[0].plot(x, [r[3] for r in rows], lw=1, ls='--', alpha=.7, label='vs Dr_kcnj10a')
    ax[0].plot(x, [r[4] for r in rows], lw=1, ls=':', color='gray',
               label='control: Bs_kcnj10a vs Dr_kcnj10a')
    ax[0].set_ylabel('% identity')
    figstyle.title(ax[0], f'Sliding-window identity for {FOCAL} '
                          f'(window {WINDOW} aa, step {STEP} aa)')
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=.3)
    ax[1].axhline(0, color='k', lw=.8)
    ax[1].fill_between(x, delta, 0, where=[d > 0 for d in delta], alpha=.6,
                       label='closer to Bs_kcnj10a')
    ax[1].fill_between(x, delta, 0, where=[d < 0 for d in delta], alpha=.6,
                       label='closer to Dr_kcnj10b')
    ax[1].set_xlabel('alignment position (aa)')
    ax[1].set_ylabel('difference (points)')
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=.3)
    fig.tight_layout()
    figstyle.save(fig, 'figures/figS2_gene_conversion.png')


if __name__ == '__main__':
    main()
