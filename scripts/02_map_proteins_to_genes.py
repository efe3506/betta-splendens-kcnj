#!/usr/bin/env python3
"""Reduce hmmsearch protein hits to genes using the RefSeq annotation."""
import collections
import re
from urllib.parse import unquote

from common import read_fasta, write_fasta

OUT = 'results/01_identification'
DOMTBL = f'{OUT}/REF_IRK_domtbl.txt'
GFF = 'data/annotation/REF_genomic.gff'
PROTEINS = 'data/annotation/REF_protein.faa'


def attr(field, text):
    m = re.search(rf'{field}=([^;]+)', text)
    return unquote(m.group(1)) if m else 'NA'


def main():
    hits, per_profile = set(), collections.defaultdict(set)
    for line in open(DOMTBL):
        if line.startswith('#'):
            continue
        f = line.split()
        hits.add(f[0])
        per_profile[f[3]].add(f[0])
    print(f'{len(hits)} proteins with a Kir domain hit')
    for profile, prots in sorted(per_profile.items()):
        print(f'  {profile:<6} {len(prots)} proteins')

    prot2gene, genes = {}, {}
    for line in open(GFF):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        if f[2] == 'gene':
            gid = re.search(r'GeneID:(\d+)', f[8])
            genes[attr('gene', f[8])] = (f[0], f[3], f[4], f[6],
                                         gid.group(1) if gid else 'NA',
                                         attr('description', f[8]))
        elif f[2] == 'CDS':
            pid = attr('protein_id', f[8])
            if pid in hits:
                prot2gene[pid] = (attr('gene', f[8]), f[0], attr('product', f[8]))

    with open(f'{OUT}/REF_IRK_proteins.txt', 'w') as fh:
        fh.write('\n'.join(sorted(hits)) + '\n')
    with open(f'{OUT}/REF_IRK_gene_map.tsv', 'w') as fh:
        fh.write('protein\tgene\tseqid\tproduct\n')
        for p, (g, c, d) in sorted(prot2gene.items(), key=lambda x: (x[1][0], x[0])):
            fh.write(f'{p}\t{g}\t{c}\t{d}\n')

    gene_set = sorted({g for g, _, _ in prot2gene.values()})
    with open(f'{OUT}/REF_kcnj_genes.tsv', 'w') as fh:
        fh.write('gene\tseqid\tstart\tend\tstrand\tgene_id\tdescription\n')
        for g in gene_set:
            fh.write(g + '\t' + '\t'.join(genes[g]) + '\n')
    print(f'{len(gene_set)} genes')

    records = read_fasta(PROTEINS)
    annotated = {h.split()[0] for h, _ in records if 'inward rectifier' in h.lower()}
    print(f'{len(annotated)} proteins described as "inward rectifier" in RefSeq; '
          f'{len(annotated - hits)} of them missed by the HMM search; '
          f'{len(hits - annotated)} HMM hits carry a different description')

    best = {}
    for header, seq in records:
        pid = header.split()[0]
        if pid not in prot2gene:
            continue
        gene = prot2gene[pid][0]
        if gene not in best or len(seq) > len(best[gene][1]):
            best[gene] = (f'{gene}|{pid}', seq)
    write_fasta(f'{OUT}/REF_kcnj_canonical.faa', [best[g] for g in sorted(best)])
    print(f'{len(best)} representative proteins written')


if __name__ == '__main__':
    main()
