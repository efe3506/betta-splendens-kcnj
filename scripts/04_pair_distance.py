#!/usr/bin/env python3
"""Distance and relative orientation of two genes in each assembly."""
import csv
import sys

BASE = 'results/02_cross_assembly'
ASSEMBLIES = ['REF', 'NCU', 'LOEWE', 'BGI']


def main(gene_a, gene_b):
    out = open(f'{BASE}/pair_{gene_a}_{gene_b}.tsv', 'w')
    out.write('assembly\tcontig\tdistance_bp\torientation\n')
    for assembly in ASSEMBLIES:
        loci = {}
        for r in csv.DictReader(open(f'{BASE}/loci_{assembly}.tsv'), delimiter='\t'):
            if r['gene'] in (gene_a, gene_b):
                loci[r['gene']] = (r['contig'], int(r['start']), int(r['end']), r['strand'])
        if len(loci) < 2:
            print(f'{assembly}: one of the genes was not found')
            continue
        a, b = loci[gene_a], loci[gene_b]
        if a[0] != b[0]:
            print(f'{assembly}: genes are on different contigs')
            continue
        first, second = sorted((a, b), key=lambda x: x[1])
        distance = second[1] - first[2]
        if first[3] == second[3]:
            orientation = 'tandem'
        elif first[3] == '+':
            orientation = 'convergent'
        else:
            orientation = 'divergent'
        print(f'{assembly:<6} {a[0]:<20} {distance:>8,} bp  {orientation}')
        out.write(f'{assembly}\t{a[0]}\t{distance}\t{orientation}\n')
    out.close()


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit('usage: 04_pair_distance.py <geneA> <geneB>')
    main(sys.argv[1], sys.argv[2])
