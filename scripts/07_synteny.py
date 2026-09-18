#!/usr/bin/env python3
"""Microsynteny: compare the gene neighbourhoods of kcnj loci in B. splendens and D. rerio."""
import collections
import os
import re

OUT = 'results/04_structure'
GFFS = {'Betta': 'data/annotation/REF_genomic.gff', 'Danio': 'data/annotation/DANIO_genomic.gff'}
FLANK = 10
FOCAL = {
    'Betta': ['kcnj10a', 'LOC114844899', 'kcnj16a', 'LOC114860957',
              'kcnj2a', 'LOC114861067', 'kcnj15', 'kcnj6'],
    'Danio': ['kcnj10a', 'kcnj10b', 'kcnj16a', 'kcnj2a', 'kcnj2b', 'kcnj15', 'kcnj6'],
}


def read_genes(path):
    genes = collections.defaultdict(list)
    for line in open(path):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'gene':
            continue
        m = re.search(r'gene=([^;]+)', f[8])
        if m:
            genes[f[0]].append((int(f[3]), int(f[4]), f[6], m.group(1)))
    for contig in genes:
        genes[contig].sort()
    return genes


def window(genes, focal):
    for contig, rows in genes.items():
        for i, (start, end, strand, name) in enumerate(rows):
            if name == focal:
                lo, hi = max(0, i - FLANK), min(len(rows), i + FLANK + 1)
                return contig, (start, end), rows[lo:hi], i - lo
    return None, None, [], 0


def normalise(name):
    n = name.lower()
    if re.search(r'\d[ab]$', n):
        n = n[:-1]
    return re.sub(r'\.\d+$', '', n)


def main():
    os.makedirs(OUT, exist_ok=True)
    neighbours = {}
    with open(f'{OUT}/synteny_neighbors.tsv', 'w') as fh:
        fh.write('species\tfocal\tcontig\tfocal_start\tfocal_end\tneighbor\t'
                 'n_start\tn_end\tstrand\toffset\n')
        for species, path in GFFS.items():
            genes = read_genes(path)
            print(f'{species}: {sum(len(v) for v in genes.values())} genes')
            for focal in FOCAL[species]:
                contig, coord, rows, centre = window(genes, focal)
                if contig is None:
                    print(f'  not found: {focal}')
                    continue
                neighbours[(species, focal)] = {normalise(r[3]) for r in rows} - {normalise(focal)}
                for k, (s, e, strand, name) in enumerate(rows):
                    fh.write(f'{species}\t{focal}\t{contig}\t{coord[0]}\t{coord[1]}\t'
                             f'{name}\t{s}\t{e}\t{strand}\t{k - centre}\n')

    with open(f'{OUT}/synteny_scores.tsv', 'w') as fh:
        fh.write('betta_gene\tdanio_gene\tshared_n\tshared_names\n')
        for (sp1, bg), bset in sorted(neighbours.items()):
            if sp1 != 'Betta':
                continue
            scores = []
            for (sp2, dg), dset in neighbours.items():
                if sp2 == 'Danio':
                    shared = sorted(bset & dset)
                    scores.append((len(shared), dg, shared))
            scores.sort(reverse=True)
            for n, dg, shared in scores[:3]:
                fh.write(f'{bg}\t{dg}\t{n}\t{",".join(shared)}\n')
                if n:
                    print(f'  {bg:<14} {dg:<9} {n}  {", ".join(shared)}')


if __name__ == '__main__':
    main()
