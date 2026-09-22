#!/usr/bin/env python3
"""Fig. S1: motif architecture of the 23 Kir proteins (MEME), proteins drawn to scale.

One site per motif and protein is drawn; where MEME's combined block diagram lists a second,
weaker site for the same motif (zoops model), the site with the better p-value is kept.
"""
import re


MEME_TXT = 'results/04_structure/meme/meme.txt'
FASTA = 'results/01_identification/REF_kcnj_canonical.faa'
OUT = 'figures/figS1_motif_architecture.png'

from common import PROPOSED_NAMES as RENAME, SUBFAMILY_ORDER as ORDER, subfamily

MOTIF_COLORS = {1: '#4C72B0', 2: '#C44E52', 3: '#55A868', 4: '#DD8452', 5: '#8172B3',
                6: '#937860', 7: '#64B5CD', 8: '#DA8BC3', 9: '#CCB974', 10: '#3a3a3a'}
FILTER = r'[TS]IG[YF]G'


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith('>'):
            name = line[1:].split()[0]
            seqs[name] = ''
        else:
            seqs[name] += line.strip()
    return seqs


def main():
    import figstyle  # noqa: F401
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    text = open(MEME_TXT).read()
    widths = {int(n): int(w) for n, w in re.findall(r'MEME-(\d+)\s+width =\s+(\d+)', text)}
    block = text.split('Combined block diagrams')[1].split('SEQUENCE NAME')[1]

    diagrams, current = {}, None
    for line in block.split('\n'):
        m = re.match(r'^(\S+\|\S+)\s+[\d.e+-]+\s+(.*)$', line)
        if m:
            current = m.group(1)
            diagrams[current] = m.group(2).rstrip('\\')
        elif current and line.startswith('    '):
            diagrams[current] += line.strip().rstrip('\\')

    seqs = read_fasta(FASTA)
    full_id = {k[:len(d)]: k for k in seqs for d in diagrams if k.startswith(d)}

    proteins = []
    for seq_id, diagram in diagrams.items():
        seq = seqs[full_id.get(seq_id, seq_id)]
        gene = RENAME.get(seq_id.split('|')[0], seq_id.split('|')[0])
        pos, best = 0, {}
        for token in diagram.split('_'):
            m = re.match(r'[\[<](\d+)\(([\d.e+-]+)\)', token)
            if m:
                motif, pval = int(m.group(1)), float(m.group(2))
                # zoops: one site per motif; keep the best-scoring one
                if motif not in best or pval < best[motif][0]:
                    best[motif] = (pval, pos)
                pos += widths[motif]
            elif token.isdigit():
                pos += int(token)
        sites = [(motif, p, widths[motif]) for motif, (_, p) in best.items()]
        hit = re.search(FILTER, seq)
        proteins.append({'gene': gene, 'sub': subfamily(gene), 'length': len(seq),
                         'sites': sites, 'filter': hit.start() if hit else -1})
    proteins.sort(key=lambda p: (ORDER.index(p['sub']), p['gene']))

    fig, ax = plt.subplots(figsize=(11, 9.5))
    y = 0
    previous = None
    ticks = []
    for p in proteins:
        if previous and p['sub'] != previous:
            y += 0.6
        previous = p['sub']
        ax.plot([0, p['length']], [y, y], color='#9a9a9a', lw=1.4, zorder=1)
        for motif, start, width in p['sites']:
            ax.add_patch(plt.Rectangle((start, y - .32), width, .64, zorder=2,
                                       facecolor=MOTIF_COLORS[motif], edgecolor='white',
                                       linewidth=.5))
            ax.text(start + width / 2, y, str(motif), ha='center', va='center', fontsize=6.5,
                    color='white', fontweight='bold', zorder=3)
        if p['filter'] >= 0:
            ax.plot(p['filter'] + 2.5, y - .50, marker='^', ms=4.5, color='k', zorder=4)
        ax.text(-12, y, p['gene'], ha='right', va='center', fontsize=9.5, fontstyle='italic')
        ax.text(p['length'] + 8, y, f"{p['length']} aa", ha='left', va='center', fontsize=7.5,
                color='#666666')
        ticks.append((p['sub'], y))
        y += 1

    for sub in ORDER:
        ys = [t for s, t in ticks if s == sub]
        ax.plot([-95, -95], [min(ys) - .35, max(ys) + .35], color='#555555', lw=2)
        ax.text(-103, (min(ys) + max(ys)) / 2, sub, ha='right', va='center', fontsize=10,
                fontweight='bold', color='#333333')

    ax.set_ylim(y - .3, -1)
    ax.set_xlim(-150, 660)
    ax.set_yticks([])
    ax.set_xticks(range(0, 601, 100))
    ax.set_xlabel('position (amino acids)')
    for side in ('left', 'right', 'top'):
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_bounds(0, 600)

    handles = [Patch(facecolor=MOTIF_COLORS[m], label=f'M{m}') for m in sorted(MOTIF_COLORS)]
    handles.append(plt.Line2D([], [], marker='^', ls='', color='k', ms=6,
                              label='selectivity filter'))
    ax.legend(handles=handles, ncol=11, fontsize=8.5, frameon=False, loc='upper center',
              bbox_to_anchor=(0.45, -0.07), handlelength=1.2, columnspacing=1.0,
              handletextpad=.4)
    fig.savefig(OUT)
    print(OUT)


if __name__ == '__main__':
    main()
