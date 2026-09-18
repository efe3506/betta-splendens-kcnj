#!/usr/bin/env python3
"""Give the 70 proteins short, readable labels for alignment and tree building."""
import collections
import os
import re

from common import read_fasta

SRC = 'results/02_cross_assembly/all_queries.faa'
OUT = 'results/03_phylogeny'


def make_label(header):
    qid = header.split()[0]
    if '|XP_' in qid:
        return 'Bs_' + qid.split('|')[0]
    m = re.search(r'GN=(\S+)', header)
    gene = m.group(1) if m else qid.replace('|', '_')
    if 'OS=Danio rerio' in header:
        return f'Dr_{gene}'
    if 'OS=Homo sapiens' in header:
        return f'Hs_{gene}'
    return f'Xx_{gene}'


def main():
    os.makedirs(OUT, exist_ok=True)
    seen, per_species = collections.Counter(), collections.Counter()
    with open(f'{OUT}/kcnj_labeled.faa', 'w') as fa, open(f'{OUT}/label_map.tsv', 'w') as mp:
        mp.write('label\toriginal_header\n')
        for header, seq in read_fasta(SRC):
            label = make_label(header)
            seen[label] += 1
            if seen[label] > 1:
                label = f'{label}_{seen[label]}'
            per_species[label[:2]] += 1
            fa.write(f'>{label}\n{seq}\n')
            mp.write(f'{label}\t{header}\n')
    print(', '.join(f'{k}: {v}' for k, v in sorted(per_species.items())))


if __name__ == '__main__':
    main()
