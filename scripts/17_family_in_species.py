#!/usr/bin/env python3
"""Kir genes of five further species, found with the same domain search as in B. splendens."""
import collections
import os
import re
import subprocess

from common import read_fasta, write_fasta

SPECIES = {'Lo': 'Lepisosteus oculatus', 'At': 'Anabas testudineus', 'Ol': 'Oryzias latipes',
           'Am': 'Astyanax mexicanus', 'Dr': 'Danio rerio'}
IN_TREE = ['Lo', 'At', 'Ol', 'Am']      # zebrafish is already in the tree (UniProt sequences)
HMM = 'data/queries/IRK_domains.hmm'
BASE = 'results/03_phylogeny/kcnj_labeled.faa'
OUT = 'results/07_taxa'
MIN_LEN, MAX_LEN = 150, 600


def protein_to_gene(gff):
    mapping = {}
    with open(gff) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'CDS':
                continue
            pid = re.search(r'protein_id=([^;]+)', f[8])
            gene = re.search(r'(?:^|;)gene=([^;]+)', f[8])
            if pid and gene:
                mapping[pid.group(1)] = gene.group(1)
    return mapping


def main():
    os.makedirs(OUT, exist_ok=True)
    added, counts = [], []
    for tag, species in SPECIES.items():
        proteome = f'data/species/{tag}/protein.faa'
        seqs = {h.split()[0]: s for h, s in read_fasta(proteome)}
        genes = protein_to_gene(f'data/species/{tag}/genomic.gff')
        table = f'{OUT}/{tag}_IRK_tbl.txt'
        subprocess.run(['hmmsearch', '--cut_ga', '--noali', '--tblout', table, HMM, proteome],
                       stdout=subprocess.DEVNULL, check=True)
        with open(table) as fh:
            hits = {line.split()[0] for line in fh if not line.startswith('#')}

        by_gene = collections.defaultdict(list)
        for pid in hits:
            by_gene[genes.get(pid, pid)].append(pid)

        kept = 0
        with open(f'{OUT}/{tag}_kir_genes.tsv', 'w') as fh:
            fh.write('gene\tprotein\tlength_aa\tn_isoforms\tin_length_range\n')
            for gene, pids in sorted(by_gene.items()):
                in_range = [p for p in pids if MIN_LEN <= len(seqs[p]) <= MAX_LEN]
                rep = max(in_range or pids, key=lambda p: len(seqs[p]))
                fh.write(f'{gene}\t{rep}\t{len(seqs[rep])}\t{len(pids)}\t'
                         f'{"yes" if in_range else "no"}\n')
                if in_range:
                    kept += 1
                    if tag in IN_TREE:
                        added.append((f'{tag}_{gene}', seqs[rep]))
                else:
                    print(f'  {tag} {gene}: outside {MIN_LEN}-{MAX_LEN} aa, not used in the tree')
        counts.append((tag, species, len(hits), len(by_gene), kept))
        print(f'{tag} ({species}): {len(hits)} proteins, {len(by_gene)} genes')

    with open(f'{OUT}/gene_counts.tsv', 'w') as fh:
        fh.write('tag\tspecies\tproteins\tgenes\tgenes_in_length_range\n')
        for row in counts:
            fh.write('\t'.join(map(str, row)) + '\n')

    base = read_fasta(BASE)
    write_fasta(f'{OUT}/kir_all.faa', base + added)
    print(f'{len(base)} + {len(added)} sequences written to {OUT}/kir_all.faa')

    with open(f'{OUT}/Dr_kir_genes.tsv') as fh:
        hmm = {line.split('\t')[0] for line in list(fh)[1:]}
    uniprot = {h[3:] for h, _ in base if h.startswith('Dr_')}
    with open(f'{OUT}/Dr_hmm_vs_uniprot.tsv', 'w') as fh:
        fh.write('gene\tin_domain_search\tin_uniprot_set\n')
        for gene in sorted(hmm | uniprot):
            fh.write(f'{gene}\t{"yes" if gene in hmm else "no"}\t'
                     f'{"yes" if gene in uniprot else "no"}\n')
    print(f'zebrafish: {len(hmm)} genes by domain search, {len(uniprot)} in the UniProt set, '
          f'{len(hmm ^ uniprot)} differences')


if __name__ == '__main__':
    main()
