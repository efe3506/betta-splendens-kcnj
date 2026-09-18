#!/usr/bin/env python3
"""Collapse miniprot alignments into loci and compare loci across assemblies."""
import collections
import re

from common import read_fasta

BASE = 'results/02_cross_assembly'
ASSEMBLIES = ['REF', 'NCU', 'LOEWE', 'BGI']
MIN_IDENTITY = 0.40
MIN_OVERLAP = 0.50


def query_labels():
    label, length = {}, {}
    for header, seq in read_fasta(f'{BASE}/all_queries.faa'):
        qid = header.split()[0]
        length[qid] = len(seq)
        if '|XP_' in qid:
            label[qid] = qid.split('|')[0]
        else:
            m = re.search(r'GN=(\S+)', header)
            label[qid] = m.group(1) if m else qid
    return label, length


def read_alignments(assembly):
    rows = []
    for line in open(f'{BASE}/{assembly}_miniprot.gff'):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'mRNA':
            continue
        target = re.search(r'Target=(\S+) (\d+) (\d+)', f[8])
        identity = re.search(r'Identity=([0-9.]+)', f[8])
        identity = float(identity.group(1)) if identity else 0.0
        if identity < MIN_IDENTITY or not target:
            continue
        fs = re.search(r'Frameshift=(\d+)', f[8])
        stop = re.search(r'StopCodon=(\d+)', f[8])
        rows.append({
            'contig': f[0], 'start': int(f[3]), 'end': int(f[4]), 'strand': f[6],
            'score': float(f[5]), 'identity': identity, 'query': target.group(1),
            'qstart': int(target.group(2)), 'qend': int(target.group(3)),
            'frameshift': int(fs.group(1)) if fs else 0,
            'stop': int(stop.group(1)) if stop else 0,
        })
    return rows


def overlaps_seed(locus, r):
    if r['contig'] != locus['contig'] or r['strand'] != locus['strand']:
        return False
    overlap = min(locus['seed_end'], r['end']) - max(locus['seed_start'], r['start'])
    shorter = min(locus['seed_end'] - locus['seed_start'], r['end'] - r['start'])
    return shorter > 0 and overlap > MIN_OVERLAP * shorter


def build_loci(rows):
    loci = []
    for r in sorted(rows, key=lambda x: -x['score']):
        placed = False
        for locus in loci:
            if overlaps_seed(locus, r):
                locus['hits'].append(r)
                placed = True
                break
        if not placed:
            loci.append({'contig': r['contig'], 'strand': r['strand'],
                         'seed_start': r['start'], 'seed_end': r['end'], 'hits': [r]})
    return sorted(loci, key=lambda l: (l['contig'], l['seed_start']))


def main():
    label, qlen = query_labels()
    counts = collections.defaultdict(dict)
    quality = collections.defaultdict(dict)

    for assembly in ASSEMBLIES:
        loci = build_loci(read_alignments(assembly))
        with open(f'{BASE}/loci_{assembly}.tsv', 'w') as fh:
            fh.write('contig\tstart\tend\tstrand\tgene\tidentity\tqcov\tframeshift\tstop\t'
                     'span_bp\tscore\tn_hits\tother_labels\n')
            for locus in loci:
                best = locus['hits'][0]
                gene = label.get(best['query'], best['query'])
                qcov = (best['qend'] - best['qstart'] + 1) / qlen[best['query']]
                others = sorted({label.get(h['query'], h['query']) for h in locus['hits']} - {gene})
                fh.write(f"{best['contig']}\t{best['start']}\t{best['end']}\t{best['strand']}\t"
                         f"{gene}\t{best['identity']:.3f}\t{qcov:.3f}\t{best['frameshift']}\t"
                         f"{best['stop']}\t{best['end'] - best['start']}\t{best['score']:.0f}\t"
                         f"{len(locus['hits'])}\t{','.join(others)}\n")
                counts[gene][assembly] = counts[gene].get(assembly, 0) + 1
                prev = quality[gene].get(assembly)
                if prev is None or best['identity'] > prev[0]:
                    quality[gene][assembly] = (best['identity'], qcov,
                                               best['frameshift'], best['stop'])
        print(f'{assembly}: {len(loci)} loci')

    with open(f'{BASE}/presence_matrix.tsv', 'w') as fh:
        fh.write('gene\t' + '\t'.join(ASSEMBLIES) + '\n')
        for gene in sorted(counts):
            fh.write(gene + '\t' + '\t'.join(str(counts[gene].get(a, 0)) for a in ASSEMBLIES) + '\n')

    with open(f'{BASE}/identity_matrix.tsv', 'w') as fh:
        fh.write('gene\t' + '\t'.join(f'{a}_ident\t{a}_qcov' for a in ASSEMBLIES) + '\n')
        for gene in sorted(quality):
            cells = []
            for a in ASSEMBLIES:
                v = quality[gene].get(a)
                cells.append(f'{v[0] * 100:.1f}\t{v[1] * 100:.1f}' if v else 'NA\tNA')
            fh.write(gene + '\t' + '\t'.join(cells) + '\n')

    flagged = 0
    with open(f'{BASE}/frameshift_flags.tsv', 'w') as fh:
        fh.write('gene\tassembly\tidentity\tframeshifts\tstop_codons\n')
        for gene in sorted(quality):
            for a in ASSEMBLIES:
                v = quality[gene].get(a)
                if v and (v[2] or v[3]):
                    fh.write(f'{gene}\t{a}\t{v[0] * 100:.1f}\t{v[2]}\t{v[3]}\n')
                    flagged += 1
    print(f'{len(counts)} genes; {flagged} gene/assembly pairs flagged for frameshift or stop')


if __name__ == '__main__':
    main()
