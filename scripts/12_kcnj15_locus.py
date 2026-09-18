#!/usr/bin/env python3
"""Compare the kcnj15 locus across the four assemblies."""
import os
import subprocess
import sys

from Bio.Seq import Seq

from common import read_fasta

OUT = 'results/05_kcnj15_locus'
ASSEMBLIES = ['REF', 'NCU', 'LOEWE', 'BGI']
ANCHOR, ANCHOR_CONTIG = 'NCU', 'CM045483.1'
UTR_VARIANT, LEAD_VARIANT = 9596738, 9590677
CDS_START, CDS_END = 9595158, 9596181
PAD = 200


def genome(name):
    return f'data/genomes/{name}.fna'


def fetch(name, region):
    out = subprocess.run(['samtools', 'faidx', genome(name), region],
                         capture_output=True, text=True, check=True).stdout
    return ''.join(out.split('\n')[1:]).upper()


def revcomp(seq):
    return seq.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def first_base(seq):
    return len(seq) - len(seq.lstrip('-'))


def last_base(seq):
    return len(seq.rstrip('-')) - 1


def main(mode):
    os.makedirs(OUT, exist_ok=True)
    if mode == 'utr':
        centre, flank = UTR_VARIANT, 1000
    elif mode == 'lead':
        centre, flank = LEAD_VARIANT, 1000
    elif mode == 'cds':
        centre = (CDS_START + CDS_END) // 2
        flank = (CDS_END - CDS_START) // 2 + 400
    else:
        sys.exit('usage: 12_kcnj15_locus.py utr|lead|cds')

    start = centre - flank
    anchor_seq = fetch(ANCHOR, f'{ANCHOR_CONTIG}:{start}-{centre + flank}')
    anchor_path = f'{OUT}/anchor_{mode}.fa'
    with open(anchor_path, 'w') as fh:
        fh.write(f'>{ANCHOR}\n{anchor_seq}\n')

    regions = {ANCHOR: anchor_seq}
    for name in ASSEMBLIES:
        if name == ANCHOR:
            continue
        paf = subprocess.run(['minimap2', '-x', 'asm10', '-t', '4', genome(name), anchor_path],
                             capture_output=True, text=True, check=True).stdout
        hits = [line.split('\t') for line in paf.strip().split('\n') if line]
        if not hits:
            print(f'{name}: no alignment')
            continue
        best = max(hits, key=lambda f: int(f[9]))
        contig, t0, t1, strand = best[5], int(best[7]), int(best[8]), best[4]
        seq = fetch(name, f'{contig}:{max(1, t0 - PAD)}-{t1 + PAD}')
        regions[name] = revcomp(seq) if strand == '-' else seq
        print(f'{name}: {contig}:{t0}-{t1} ({strand})')

    fasta = f'{OUT}/kcnj15_{mode}_regions.fa'
    with open(fasta, 'w') as fh:
        for name in ASSEMBLIES:
            if name in regions:
                fh.write(f'>{name}\n{regions[name]}\n')
    aln_path = f'{OUT}/kcnj15_{mode}_aln.fa'
    with open(aln_path, 'w') as fh:
        subprocess.run(['mafft', '--auto', '--quiet', fasta], stdout=fh, check=True)
    aln = {h: s.upper() for h, s in read_fasta(aln_path)}
    names = [n for n in ASSEMBLIES if n in aln]

    column = [i for i, ch in enumerate(aln[ANCHOR]) if ch != '-']
    target = column[centre - start]

    print(f'\nanchor position {ANCHOR_CONTIG}:{centre} = alignment column {target + 1}')
    with open(f'{OUT}/variant_genotypes_{mode}.tsv', 'w') as fh:
        fh.write('assembly\tallele\tcontext\n')
        for n in names:
            context = f'{aln[n][target - 10:target]}[{aln[n][target]}]{aln[n][target + 1:target + 11]}'
            print(f'  {n:<6} {aln[n][target]}   {context}')
            fh.write(f'{n}\t{aln[n][target]}\t{context}\n')

    first = max(first_base(aln[n]) for n in names)
    last = min(last_base(aln[n]) for n in names)
    with open(f'{OUT}/all_differences_{mode}.tsv', 'w') as fh:
        fh.write('aln_column\t' + '_'.join(names) + '\ttype\n')
        for i in range(first, last + 1):
            chars = [aln[n][i] for n in names]
            if len(set(chars)) > 1:
                kind = 'indel' if '-' in chars else 'substitution'
                fh.write(f'{i + 1}\t{"".join(chars)}\t{kind}\n')
                print(f'  column {i + 1:>5}  {"".join(chars)}  {kind}'
                      f'{"   <- reported variant" if i == target else ""}')
                if kind == 'indel':
                    ref = next(n for n in names if aln[n][i] != '-')
                    print(f'      context in {ref}: {aln[ref][i - 12:i]}[{aln[ref][i]}]'
                          f'{aln[ref][i + 1:i + 13]}')

    if mode == 'cds':
        c0, c1 = column[CDS_START - start], column[CDS_END - start]
        print('\ntranslation of the coding region (first stop codon ends the product):')
        with open(f'{OUT}/cds_translation.tsv', 'w') as fh:
            fh.write('assembly\tcds_bp\tproduct_aa\n')
            for n in names:
                cds = aln[n][c0:c1 + 1].replace('-', '')
                protein = str(Seq(cds[:len(cds) // 3 * 3]).translate())
                product = protein.split('*')[0]
                print(f'  {n:<6} {len(cds)} bp -> {len(product)} aa')
                fh.write(f'{n}\t{len(cds)}\t{len(product)}\n')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '')
