#!/usr/bin/env python3
"""Targeted analysis of AlphaFold DB models (no structure prediction is run here)."""
import json
import os
import subprocess
import urllib.request

import numpy as np

from common import read_fasta

OUT = 'results/06_protein_model'
TARGETS = [('A0A6P7NPA5', 'kcnj4b'), ('A0A6P7L678', 'kcnj2a'),
           ('A0A6P7N863', 'kcnj16b'), ('A0A6P7L4L5', 'kcnj16a'),
           ('A0A6P7L185', 'kcnj10b'), ('A0A6P7LKZ5', 'kcnj10a'),
           ('A0A6P7PE36', 'kcnj15')]
PAIRS = [('kcnj16b', 'kcnj16a'), ('kcnj10b', 'kcnj10a'), ('kcnj4b', 'kcnj2a')]
FILTER = 'TIGYG'
CONFIDENT = 70.0
AA = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E',
      'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
      'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}


def fetch_model(accession, gene):
    path = f'{OUT}/{gene}_{accession}.pdb'
    if not os.path.exists(path):
        meta = json.load(urllib.request.urlopen(
            f'https://alphafold.ebi.ac.uk/api/prediction/{accession}'))[0]
        urllib.request.urlretrieve(meta['pdbUrl'], path)
        print(f'  downloaded {gene} ({accession}, model v{meta["latestVersion"]})')
    return path


def read_ca(path):
    seq, plddt, xyz = [], [], []
    for line in open(path):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            seq.append(AA.get(line[17:20], 'X'))
            plddt.append(float(line[60:66]))
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return ''.join(seq), np.array(plddt), np.array(xyz)


def align(name_a, seq_a, name_b, seq_b):
    tmp = f'{OUT}/_pair.fa'
    with open(tmp, 'w') as fh:
        fh.write(f'>{name_a}\n{seq_a}\n>{name_b}\n{seq_b}\n')
    out = subprocess.run(['mafft', '--localpair', '--maxiterate', '1000', '--quiet', tmp],
                         capture_output=True, text=True, check=True).stdout
    with open(tmp, 'w') as fh:
        fh.write(out)
    aln = dict(read_fasta(tmp))
    os.remove(tmp)
    return aln[name_a].upper(), aln[name_b].upper()


def rmsd_after_superposition(P, Q):
    P, Q = P - P.mean(0), Q - Q.mean(0)
    V, _, W = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(V @ W))
    R = V @ np.diag([1, 1, d]) @ W
    return float(np.sqrt(((P @ R - Q) ** 2).sum() / len(P)))


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs('figures', exist_ok=True)
    models = {gene: read_ca(fetch_model(acc, gene)) for acc, gene in TARGETS}

    with open(f'{OUT}/plddt_summary.tsv', 'w') as fh:
        fh.write('protein\tlength\tmean_plddt\tpct_below_50\tfilter_position\tfilter_plddt\n')
        for _, gene in TARGETS:
            seq, plddt, _ = models[gene]
            i = seq.find(FILTER)
            row = (gene, len(seq), plddt.mean(), 100 * (plddt < 50).mean(),
                   i + 1, plddt[i:i + 5].mean())
            fh.write('{}\t{}\t{:.1f}\t{:.1f}\t{}\t{:.1f}\n'.format(*row))
            print('  {:<8} {:>4} aa  mean pLDDT {:.1f}  <50: {:.1f}%  filter at {} '
                  '(pLDDT {:.1f})'.format(*row))

    with open(f'{OUT}/structural_comparison.tsv', 'w') as fh:
        fh.write('target\treference\trmsd_A\tn_residues\tseq_identity_pct\n')
        for target, ref in PAIRS:
            (s1, p1, x1), (s2, p2, x2) = models[target], models[ref]
            a1, a2 = align(target, s1, ref, s2)
            i = j = 0
            P, Q, same, cols = [], [], 0, 0
            for c1, c2 in zip(a1, a2):
                if c1 != '-' and c2 != '-':
                    cols += 1
                    same += c1 == c2
                    if p1[i] >= CONFIDENT and p2[j] >= CONFIDENT:
                        P.append(x1[i])
                        Q.append(x2[j])
                if c1 != '-':
                    i += 1
                if c2 != '-':
                    j += 1
            rmsd = rmsd_after_superposition(np.array(P), np.array(Q))
            fh.write(f'{target}\t{ref}\t{rmsd:.2f}\t{len(P)}\t{100 * same / cols:.1f}\n')
            print(f'  {target} vs {ref}: RMSD {rmsd:.2f} A over {len(P)} residues, '
                  f'identity {100 * same / cols:.1f}%')

    import figstyle
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(TARGETS), 1, figsize=(11, 1.5 * len(TARGETS)), sharex=True)
    for ax, (_, gene) in zip(axes, TARGETS):
        seq, plddt, _ = models[gene]
        x = np.arange(1, len(seq) + 1)
        ax.fill_between(x, plddt, 0, where=plddt >= CONFIDENT, color='#4C72B0', alpha=.8)
        ax.fill_between(x, plddt, 0, where=plddt < CONFIDENT, color='#C44E52', alpha=.8)
        ax.axhline(CONFIDENT, color='k', lw=.6, ls='--')
        i = seq.find(FILTER)
        ax.axvspan(i + 1, i + 5, color='gold', alpha=.7)
        ax.set_ylim(0, 100)
        ax.set_yticks([50, 100])
        ax.set_ylabel(f'$\\it{{{gene}}}$', rotation=0, ha='right', va='center', fontsize=10)
        ax.tick_params(labelsize=7)
    axes[-1].set_xlabel('residue number')
    figstyle.title(axes[0], 'AlphaFold pLDDT profiles', fontsize=11)
    fig.tight_layout()
    figstyle.save(fig, 'figures/figS3_plddt_profiles.png')


if __name__ == '__main__':
    main()
