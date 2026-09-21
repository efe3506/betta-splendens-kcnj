#!/usr/bin/env python3
"""Assign protein-to-genome alignments to annotated genes.

usage: 18_hits_to_genes.py <miniprot.tsv> <out_prefix> TAG=annotation.gff [TAG=annotation.gff ...]

The input has the columns genome, sequence, start, end, strand, attributes (miniprot mRNA
lines). Genomes without an annotation are summarised only.
"""
import collections
import re
import sys

from common import display_name

PAD = 5000


def protein_coding_genes(gff):
    genes = collections.defaultdict(list)
    with open(gff) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'gene' or 'gene_biotype=protein_coding' not in f[8]:
                continue
            name = re.search(r'(?:^|;)Name=([^;]+)', f[8])
            genes[f[0]].append((int(f[3]), int(f[4]), name.group(1) if name else '?'))
    return genes


def main():
    hits_path, prefix = sys.argv[1], sys.argv[2]
    annotations = {tag: protein_coding_genes(path)
                   for tag, path in (arg.split('=', 1) for arg in sys.argv[3:])}

    rows = []
    with open(hits_path) as fh:
        for line in fh:
            genome, seqid, start, end, strand, attr = line.rstrip('\n').split('\t')
            target = re.search(r'Target=(\S+) (\d+) (\d+)', attr)
            rows.append({'genome': genome, 'seqid': seqid, 'start': int(start), 'end': int(end),
                         'query': target.group(1), 'span': f'{target.group(2)}-{target.group(3)}',
                         'identity': 100 * float(re.search(r'Identity=([\d.]+)', attr).group(1)),
                         'frameshifts': (re.search(r'Frameshift=(\d+)', attr) or [0, '0'])[1]})

    unassigned = 0
    with open(f'{prefix}_hits.tsv', 'w') as fh:
        fh.write('genome\tquery\tseqid\tstart\tend\tquery_span\tidentity_pct\tframeshifts\tgene\n')
        for r in sorted(rows, key=lambda r: (r['genome'], r['query'], -r['identity'])):
            if r['genome'] not in annotations:
                continue
            # a Kir gene if one overlaps, otherwise whatever overlaps
            over = [display_name(name) for s, e, name in annotations[r['genome']][r['seqid']]
                    if s - PAD <= r['end'] and r['start'] <= e + PAD]
            kir = [g for g in over if re.match(r'kcnj\d', g)]
            gene = ','.join(kir or over) or 'NONE'
            unassigned += not kir
            fh.write(f"{r['genome']}\t{r['query']}\t{r['seqid']}\t{r['start']}\t{r['end']}\t"
                     f"{r['span']}\t{r['identity']:.1f}\t{r['frameshifts']}\t{gene}\n")

    summary = collections.defaultdict(list)
    for r in rows:
        summary[(r['query'], r['genome'])].append(r['identity'])
    with open(f'{prefix}_summary.tsv', 'w') as fh:
        fh.write('query\tgenome\tn_alignments\tbest_identity_pct\n')
        for (query, genome), values in sorted(summary.items()):
            fh.write(f'{query}\t{genome}\t{len(values)}\t{max(values):.1f}\n')
    print(f'{len(rows)} alignments; outside an annotated kcnj gene: {unassigned}')


if __name__ == '__main__':
    main()
