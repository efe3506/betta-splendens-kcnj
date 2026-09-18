"""Shared constants and helpers used by the analysis scripts."""
import re

PROPOSED_NAMES = {
    'LOC114861067': 'kcnj2b',
    'LOC114863156': 'kcnj4',
    'LOC114860957': 'kcnj16b',
    'LOC114846053': 'kcnj19b',
    'LOC114860740': 'kcnj19a',
    'LOC114844899': 'kcnj10c',
}

SUBFAMILY_BY_NUMBER = {
    1: 'Kir1',
    2: 'Kir2', 4: 'Kir2', 12: 'Kir2', 14: 'Kir2', 17: 'Kir2', 18: 'Kir2',
    3: 'Kir3', 5: 'Kir3', 6: 'Kir3', 9: 'Kir3', 19: 'Kir3', 20: 'Kir3', 21: 'Kir3',
    10: 'Kir4', 15: 'Kir4',
    16: 'Kir5',
    8: 'Kir6', 11: 'Kir6',
    13: 'Kir7',
}
SUBFAMILY_ORDER = ['Kir1', 'Kir2', 'Kir3', 'Kir4', 'Kir5', 'Kir6', 'Kir7']
SUBFAMILY_COLORS = dict(zip(SUBFAMILY_ORDER, [
    '#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B3', '#937860', '#DA8BC3']))


def display_name(gene):
    return PROPOSED_NAMES.get(gene, gene)


def subfamily(gene):
    m = re.match(r'kcnj(\d+)', gene.lower())
    return SUBFAMILY_BY_NUMBER.get(int(m.group(1))) if m else None


def read_fasta(path):
    records, header, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(chunks)))
                header, chunks = line[1:], []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        records.append((header, ''.join(chunks)))
    return records


def write_fasta(path, records):
    with open(path, 'w') as fh:
        for header, seq in records:
            fh.write(f'>{header}\n{seq}\n')
