#!/usr/bin/env python3
"""Build the Kir query set from a raw UniProt download."""
import re
import sys

from common import read_fasta, write_fasta

MIN_LEN, MAX_LEN = 150, 550


def gene_symbol(header):
    m = re.search(r'GN=(\S+)', header)
    return m.group(1) if m else header


def main(raw, out):
    records = read_fasta(raw)
    by_gene = {}
    for header, seq in records:
        by_gene.setdefault(gene_symbol(header), []).append((header, seq))

    kept, reinstated, dropped = {}, [], []
    for gene, recs in by_gene.items():
        ok = [r for r in recs if MIN_LEN <= len(r[1]) <= MAX_LEN]
        if ok:
            kept[gene] = max(ok, key=lambda r: len(r[1]))
        elif all(len(r[1]) > MAX_LEN for r in recs):
            kept[gene] = min(recs, key=lambda r: len(r[1]))
            reinstated.append(gene)
        else:
            dropped.append(gene)

    write_fasta(out, [kept[g] for g in sorted(kept)])
    print(f'{len(records)} records, {len(by_gene)} gene symbols in {raw}')
    print(f'{len(kept)} sequences written to {out}')
    if reinstated:
        print(f'reinstated (all records > {MAX_LEN} aa): {", ".join(sorted(reinstated))}')
    if dropped:
        print(f'dropped (only fragments < {MIN_LEN} aa): {", ".join(sorted(dropped))}')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit('usage: 01_filter_queries.py <raw_uniprot.faa> <out.faa>')
    main(sys.argv[1], sys.argv[2])
