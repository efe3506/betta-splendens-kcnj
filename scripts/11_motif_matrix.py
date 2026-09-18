#!/usr/bin/env python3
"""Motif presence/absence matrix from the MEME text output."""
import re

from common import display_name

MEME_TXT = 'results/04_structure/meme/meme.txt'
OUT = 'results/04_structure/motif_presence.tsv'


def main():
    text = open(MEME_TXT).read()
    n_motifs = len(re.findall(r'^MOTIF ', text, re.M))
    block = text.split('Combined block diagrams')[1]

    diagrams, current = {}, None
    for line in block.split('\n'):
        m = re.match(r'^(\S+\|\S+)\s+[\d.e+-]+\s+(.*)$', line)
        if m:
            current = m.group(1)
            diagrams[current] = m.group(2)
        elif current and line.startswith('    '):
            diagrams[current] += line.strip()

    rows = []
    for seq_id, diagram in diagrams.items():
        gene = display_name(seq_id.split('|')[0])
        rows.append((gene, {int(x) for x in re.findall(r'\[(\d+)\(', diagram)}))
    rows.sort()

    with open(OUT, 'w') as fh:
        fh.write('gene\t' + '\t'.join(f'M{i}' for i in range(1, n_motifs + 1)) + '\n')
        for gene, motifs in rows:
            fh.write(gene + '\t' + '\t'.join('1' if i in motifs else '0'
                                             for i in range(1, n_motifs + 1)) + '\n')

    print(f'{len(rows)} proteins, {n_motifs} motifs')
    for i in range(1, n_motifs + 1):
        missing = [g for g, motifs in rows if i not in motifs]
        if missing and len(missing) <= len(rows) // 2:
            print(f'  M{i} absent from: {", ".join(missing)}')
        elif missing:
            present = [g for g, motifs in rows if i in motifs]
            print(f'  M{i} present only in: {", ".join(present)}')


if __name__ == '__main__':
    main()
