#!/usr/bin/env python3
"""Kir4/Kir5 subfamily sequence sets, with and without two additional percomorph species."""
import os
import subprocess

from common import read_fasta, write_fasta

SRC = 'results/03_phylogeny/kcnj_labeled.faa'
OUT = 'results/04_structure'
SUBFAMILY = ['Bs_kcnj10a', 'Bs_LOC114844899', 'Bs_kcnj15', 'Bs_kcnj16a', 'Bs_LOC114860957',
             'Dr_kcnj10a', 'Dr_kcnj10b', 'Dr_kcnj15', 'Dr_kcnj16a',
             'Hs_KCNJ10', 'Hs_KCNJ15', 'Hs_KCNJ16']
QUERIES = ['Bs_kcnj10a', 'Bs_LOC114844899']
PROTEOMES = {'At': 'data/annotation/ANABAS_protein.faa',
             'Ol': 'data/annotation/MEDAKA_protein.faa'}
TOP_N, MIN_IDENTITY = 6, 55.0


def main():
    os.makedirs(OUT, exist_ok=True)
    seqs = dict(read_fasta(SRC))
    subset = [(name, seqs[name]) for name in SUBFAMILY]
    write_fasta(f'{OUT}/kir45.faa', subset)

    query_path = f'{OUT}/kcnj10_blast_query.faa'
    write_fasta(query_path, [(q, seqs[q]) for q in QUERIES])

    added, hit_rows = [], []
    for tag, proteome in PROTEOMES.items():
        table = subprocess.run(
            ['blastp', '-query', query_path, '-subject', proteome, '-evalue', '1e-20',
             '-max_hsps', '1', '-outfmt', '6 qseqid sseqid pident length evalue bitscore'],
            capture_output=True, text=True, check=True).stdout
        per_query = {}
        for line in table.strip().split('\n'):
            f = line.split('\t')
            if len(f) == 6 and float(f[2]) >= MIN_IDENTITY:
                per_query.setdefault(f[0], []).append(f)
        chosen = []
        for query in QUERIES:
            hits = sorted(per_query.get(query, []), key=lambda h: -float(h[5]))[:TOP_N]
            for h in hits:
                hit_rows.append([tag] + h)
                if h[1] not in chosen:
                    chosen.append(h[1])
        records = {hd.split()[0]: (hd, sq) for hd, sq in read_fasta(proteome)}
        for accession in chosen:
            header, seq = records[accession]
            added.append((f'{tag}_{accession.split(".")[0]}', seq))
            print(f'  {tag}_{accession:<18} {len(seq):>4} aa  {header[len(accession) + 1:][:70]}')

    write_fasta(f'{OUT}/kir45_expanded.faa', subset + added)
    with open(f'{OUT}/kcnj10_blast_hits.tsv', 'w') as fh:
        fh.write('species\tquery\tsubject\tidentity\tlength\tevalue\tbitscore\n')
        for r in hit_rows:
            fh.write('\t'.join(r) + '\n')
    print(f'{len(subset)} + {len(added)} sequences written')


if __name__ == '__main__':
    main()
