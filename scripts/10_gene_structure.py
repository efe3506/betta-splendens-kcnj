#!/usr/bin/env python3
"""Gene structure table, chromosome map (Fig. 1) and exon-intron diagram (Fig. 3)."""
import collections
import os
import re

from common import (SUBFAMILY_COLORS, SUBFAMILY_ORDER, display_name, read_fasta, subfamily)

GFF = 'data/annotation/REF_genomic.gff'
CANONICAL = 'results/01_identification/REF_kcnj_canonical.faa'
OUT = 'results/04_structure'
TANDEM_MB = 0.3


def parse_gff(wanted):
    gene_rec, rna_of_protein = {}, {}
    exon_spans = collections.defaultdict(list)
    cds_len = collections.Counter()
    seq_len, chrom_name = {}, {}
    for line in open(GFF):
        if line.startswith('##sequence-region'):
            p = line.split()
            seq_len[p[1]] = int(p[3])
            continue
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        if f[2] == 'region' and 'chromosome=' in f[8]:
            chrom_name.setdefault(f[0], re.search(r'chromosome=([^;]+)', f[8]).group(1))
            continue
        m = re.search(r'gene=([^;]+)', f[8])
        if not m or m.group(1) not in wanted:
            continue
        parent = re.search(r'Parent=([^;]+)', f[8])
        if f[2] == 'gene':
            gene_rec[m.group(1)] = (f[0], int(f[3]), int(f[4]), f[6])
        elif f[2] == 'exon' and parent:
            exon_spans[parent.group(1)].append((int(f[3]), int(f[4])))
        elif f[2] == 'CDS' and parent:
            cds_len[parent.group(1)] += int(f[4]) - int(f[3]) + 1
            pid = re.search(r'protein_id=([^;]+)', f[8])
            if pid:
                rna_of_protein[pid.group(1)] = parent.group(1)
    return gene_rec, rna_of_protein, exon_spans, cds_len, seq_len, chrom_name


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs('figures', exist_ok=True)
    canonical = dict(h.split('|') for h, _ in read_fasta(CANONICAL))
    gene_rec, rna_of_protein, exon_spans, cds_len, seq_len, chrom_name = parse_gff(set(canonical))

    rows = []
    for gene, protein in canonical.items():
        rna = rna_of_protein[protein]
        contig, start, end, strand = gene_rec[gene]
        name = display_name(gene)
        rows.append({
            'gene': name, 'refseq_id': gene if gene != name else '-', 'protein': protein,
            'subfamily': subfamily(name), 'chr': chrom_name.get(contig, '?'), 'contig': contig,
            'start': start, 'end': end, 'strand': strand, 'span': end - start + 1,
            'exons': len(exon_spans[rna]), 'cds': cds_len[rna],
            'aa': cds_len[rna] // 3 - 1, 'spans': sorted(exon_spans[rna]),
        })
    rows.sort(key=lambda r: (SUBFAMILY_ORDER.index(r['subfamily']), r['gene']))

    with open(f'{OUT}/gene_structure.tsv', 'w') as fh:
        fh.write('gene\tsubfamily\trefseq_id\tprotein\tchr\tchrom_acc\tstart\tend\tstrand\t'
                 'span_bp\texons\tcds_bp\tprotein_aa\n')
        for r in rows:
            fh.write(f"{r['gene']}\t{r['subfamily']}\t{r['refseq_id']}\t{r['protein']}\t"
                     f"{r['chr']}\t{r['contig']}\t{r['start']}\t{r['end']}\t{r['strand']}\t"
                     f"{r['span']}\t{r['exons']}\t{r['cds']}\t{r['aa']}\n")
    print(f'{len(rows)} genes written to {OUT}/gene_structure.tsv')

    import figstyle
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    chroms = sorted((c for c in chrom_name if chrom_name[c].isdigit()),
                    key=lambda c: int(chrom_name[c]))
    fig, ax = plt.subplots(figsize=(14, 6.5))
    for i, contig in enumerate(chroms):
        ax.plot([i, i], [0, seq_len[contig] / 1e6], color='#d9d9d9', lw=7,
                solid_capstyle='round', zorder=1)
        prev_y = prev_x = last_label_y = None
        for r in sorted((r for r in rows if r['contig'] == contig), key=lambda r: r['start']):
            y = r['start'] / 1e6
            tandem = prev_y is not None and abs(y - prev_y) < TANDEM_MB
            x = i + (0.13 if tandem else 0.0)
            ax.plot(x, y, 'o', color=SUBFAMILY_COLORS[r['subfamily']], ms=7, zorder=3,
                    markeredgecolor='white', markeredgewidth=.6)
            if tandem:
                ax.plot([prev_x, x], [prev_y, y], color='#333333', lw=1.4, zorder=2)
            label_y = y
            if last_label_y is not None and abs(label_y - last_label_y) < 0.9:
                label_y = last_label_y + 0.9
            if abs(label_y - y) > 0.05:
                ax.plot([x, x + 0.07], [y, label_y], color='#999999', lw=.6, zorder=2)
            ax.annotate(f"$\\it{{{r['gene']}}}$", (x, label_y), xytext=(8, 0),
                        textcoords='offset points', fontsize=8.5, va='center', zorder=4)
            prev_y, prev_x, last_label_y = y, x, label_y
    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels([chrom_name[c] for c in chroms], fontsize=8)
    ax.set_xlim(-0.7, len(chroms) - 0.3)
    ax.invert_yaxis()
    ax.set_xlabel('chromosome')
    ax.set_ylabel('position (Mb)')
    figstyle.title(ax, 'Chromosomal distribution of $\\it{kcnj}$ genes in '
                       '$\\it{Betta\\ splendens}$ (fBetSpl5.4)', fontsize=12)
    ax.legend(handles=[Line2D([], [], marker='o', ls='', color=SUBFAMILY_COLORS[s], label=s)
                       for s in SUBFAMILY_ORDER], fontsize=9, title='subfamily')
    ax.grid(axis='y', alpha=.3)
    fig.tight_layout()
    fig.savefig('figures/fig1_chromosome_map.png')

    fig, ax = plt.subplots(figsize=(10, 9))
    for k, r in enumerate(rows):
        y = len(rows) - k
        g0, g1 = r['spans'][0][0], max(e for _, e in r['spans'])
        scale = 1.0 / (g1 - g0)
        ax.plot([0, 1], [y, y], color='#999999', lw=1, zorder=1)
        for s, e in r['spans']:
            ax.add_patch(plt.Rectangle(((s - g0) * scale, y - .3), (e - s) * scale, .6,
                                       color=SUBFAMILY_COLORS[r['subfamily']], zorder=2))
        ax.text(-0.02, y, f"$\\it{{{r['gene']}}}$", ha='right', va='center', fontsize=10)
        ax.text(1.02, y, f"{r['exons']} exons · {r['span'] / 1000:.1f} kb",
                ha='left', va='center', fontsize=8, color='#555555')
    ax.set_xlim(-0.18, 1.3)
    ax.set_ylim(0, len(rows) + 1)
    ax.axis('off')
    figstyle.title(ax, 'Exon–intron structure of $\\it{kcnj}$ genes', fontsize=12)
    fig.tight_layout()
    fig.savefig('figures/fig3_gene_structure.png')


if __name__ == '__main__':
    main()
