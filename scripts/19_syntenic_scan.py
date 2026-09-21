#!/usr/bin/env python3
"""Search the syntenic interval of a missing gene with tblastn.

A protein-to-genome aligner finds intact genes; a degraded remnant is found only by searching
the interval in which the gene is expected. Flanking genes follow the percomorph arrangement:

  kcnj20 : zbtb32 | kcnj20 | igflr1         kcnj4a : kdelr3 | kcnj4 | sox10
  kcnj4b : rnaseh2a ... kcnj4b | smdt1b     kcnj16b: kcnj2b | kcnj16b ... otop2, ush1ga

usage: 19_syntenic_scan.py [betta|otophysi]
"""
import os
import subprocess
import sys

OUTDIR = 'results/09_syntenic_scan'
QUERIES = 'results/07_taxa/kir_all.faa'
EVALUE = '1e-3'
PAD = 20000          # B. splendens: flank added to the interval between the two genes
FLANK = 200000       # otophysans: searched on either side of each flanking gene

# B. splendens: coordinates in REF (end of the left gene, start of the right gene)
BETTA_INTERVALS = {
    'kcnj20_zbtb32-igflr1': ('NC_040896.2', 15037184, 15042124, ['At_kcnj20', 'Ol_kcnj20']),
    'kcnj4a_kdelr3-sox10': ('NC_040888.2', 1426570, 1435269, ['At_kcnj4', 'Ol_kcnj4']),
}
BETTA_GENOMES = {g: f'data/genomes/{g}.fna' for g in ('REF', 'NCU', 'LOEWE', 'BGI')}

# otophysans: positions of the flanking genes in the RefSeq annotations
B_QUERIES = {
    'kcnj4b': ['Bs_LOC114863156', 'At_LOC113162111', 'Ol_LOC101162982', 'Dr_kcnj4', 'Am_kcnj4'],
    'kcnj16b': ['Bs_LOC114860957', 'At_LOC113160616', 'Ol_LOC101171333', 'Dr_kcnj16a',
                'Am_kcnj16'],
}
OTOPHYSI = {
    'Dr': ('data/species/Dr/genome.fna', {
        'kcnj4b': [('smdt1b', 'NC_141021.1', 56704533, 56705546),
                   ('rnaseh2a', 'NC_141021.1', 54207408, 54215171)],
        'kcnj16b': [('kcnj2b', 'NC_141023.1', 12038559, 12048504),
                    ('otop2-ush1ga', 'NC_141023.1', 61727078, 61795283)],
    }),
    'Am': ('data/species/Am/genome.fna', {
        'kcnj4b': [('smdt1b', 'NC_064428.1', 34674404, 34677059),
                   ('rnaseh2a', 'NC_064428.1', 32569148, 32581904)],
        'kcnj16b': [('kcnj2b', 'NC_064426.1', 11253079, 11256655),
                    ('otop2-ush1ga', 'NC_064426.1', 43125276, 43167322)],
    }),
}


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, text=True, capture_output=True, **kw).stdout


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith('>'):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: ''.join(v) for k, v in seqs.items()}


def write_queries(names, path):
    allq = read_fasta(QUERIES)
    missing = [n for n in names if n not in allq]
    if missing:
        raise SystemExit(f'query not found: {missing}')
    with open(path, 'w') as fh:
        for n in names:
            fh.write(f'>{n}\n{allq[n]}\n')


def extract(genome, seqid, start, end, path):
    start = max(1, start)
    fa = run(['samtools', 'faidx', genome, f'{seqid}:{start}-{end}'])
    open(path, 'w').write(fa)
    seq = ''.join(fa.split('\n')[1:])
    return len(seq), seq.upper().count('N')


def tblastn(query, window, offset):
    out = run(['tblastn', '-query', query, '-subject', window, '-evalue', EVALUE, '-seg', 'no',
               '-outfmt', '6 qseqid qstart qend sstart send pident length evalue bitscore'])
    rows = []
    for line in out.strip().split('\n'):
        if not line:
            continue
        q, qs, qe, ss, se, pid, ln, ev, bit = line.split('\t')
        a, b = sorted((int(ss), int(se)))
        strand = '+' if int(ss) < int(se) else '-'
        rows.append([q, qs, qe, a + offset - 1, b + offset - 1, strand, pid, ln, ev, bit])
    return rows


def locate(ref_window, genome):
    """Region of another assembly that corresponds to a REF window (minimap2 asm10)."""
    paf = run(['minimap2', '-x', 'asm10', '-t', '8', genome, ref_window])
    best = {}
    for line in paf.strip().split('\n'):
        if not line:
            continue
        f = line.split('\t')
        tname, ts, te, nmatch = f[5], int(f[7]), int(f[8]), int(f[9])
        b = best.setdefault(tname, [ts, te, 0])
        b[0], b[1], b[2] = min(b[0], ts), max(b[1], te), b[2] + nmatch
    if not best:
        return None
    tname, (ts, te, _) = max(best.items(), key=lambda kv: kv[1][2])
    return tname, ts + 1, te


def scan_betta(hits, windows):
    for name, (seqid, left, right, qnames) in BETTA_INTERVALS.items():
        qpath = f'{OUTDIR}/q_{name}.faa'
        write_queries(qnames, qpath)
        ref_win = f'{OUTDIR}/REF_{name}.fna'
        s, e = left - PAD, right + PAD
        n, nn = extract(BETTA_GENOMES['REF'], seqid, s, e, ref_win)
        for g, genome in BETTA_GENOMES.items():
            if g == 'REF':
                loc, win = (seqid, s, e), ref_win
            else:
                loc = locate(ref_win, genome)
                if loc is None:
                    windows.append([g, name, '-', '-', '-', 0, 0, 'no minimap2 alignment'])
                    continue
                win = f'{OUTDIR}/{g}_{name}.fna'
                n, nn = extract(genome, loc[0], loc[1], loc[2], win)
            note = f'interval {seqid}:{left}-{right} +-{PAD}' if g == 'REF' else 'minimap2 asm10'
            windows.append([g, name, loc[0], loc[1], loc[2], n, nn, note])
            for r in tblastn(qpath, win, loc[1]):
                hits.append([g, name, loc[0]] + r)
            print(f'  {g:6} {name:24} {loc[0]}:{loc[1]}-{loc[2]}  ({n} bp, N={nn})')


def scan_otophysi(hits, windows):
    for tag, (genome, genes) in OTOPHYSI.items():
        if not os.path.exists(genome):
            print(f'skipped {tag}: {genome} not found')
            continue
        for gene, anchors in genes.items():
            qpath = f'{OUTDIR}/q_{gene}.faa'
            write_queries(B_QUERIES[gene], qpath)
            for anchor, seqid, a, b in anchors:
                name = f'{gene}_near_{anchor}'
                s, e = max(1, a - FLANK), b + FLANK
                win = f'{OUTDIR}/{tag}_{name}.fna'
                n, nn = extract(genome, seqid, s, e, win)
                windows.append([tag, name, seqid, s, e, n, nn, f'{anchor} +-{FLANK}'])
                for r in tblastn(qpath, win, s):
                    hits.append([tag, name, seqid] + r)
                print(f'  {tag:6} {name:28} {seqid}:{s}-{e}  ({n} bp, N={nn})')


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    hits, windows = [], []
    if what in ('all', 'betta'):
        scan_betta(hits, windows)
    if what in ('all', 'otophysi'):
        scan_otophysi(hits, windows)

    suffix = '' if what == 'all' else f'_{what}'
    with open(f'{OUTDIR}/windows{suffix}.tsv', 'w') as fh:
        fh.write('genome\twindow\tseqid\tstart\tend\tlength_bp\tn_count\tnote\n')
        for w in windows:
            fh.write('\t'.join(map(str, w)) + '\n')
    with open(f'{OUTDIR}/hits{suffix}.tsv', 'w') as fh:
        fh.write('genome\twindow\tseqid\tquery\tqstart\tqend\tstart\tend\tstrand\tpident\t'
                 'aln_len\tevalue\tbitscore\n')
        for h in sorted(hits, key=lambda h: (h[0], h[1], h[6])):
            fh.write('\t'.join(map(str, h)) + '\n')
    print(f'{len(hits)} alignments written to {OUTDIR}/hits{suffix}.tsv')


if __name__ == '__main__':
    main()
