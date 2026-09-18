#!/usr/bin/env python3
"""Fig. 4: the kcnj16-kcnj2 gene pair in both duplicated blocks of two species."""
import collections
import csv
import os

from common import PROPOSED_NAMES

SRC = 'results/04_structure/synteny_neighbors.tsv'
BLOCKS = [
    ('Betta splendens', 'chromosome 19  —  $\\it{a}$ block', 'Betta', ['kcnj16a', 'kcnj2a']),
    ('Danio rerio', '$\\it{kcnj16a}$/$\\it{kcnj2a}$ locus  —  $\\it{a}$ block', 'Danio',
     ['kcnj16a', 'kcnj2a']),
    ('Betta splendens', 'chromosome 8  —  $\\it{b}$ block', 'Betta',
     ['LOC114860957', 'LOC114861067']),
    ('Danio rerio', '$\\it{kcnj2b}$ locus  —  $\\it{b}$ block', 'Danio', ['kcnj2b']),
]
FOCAL_COLORS = {'kcnj16a': '#C44E52', 'kcnj16b': '#C44E52',
                'kcnj2a': '#DD8452', 'kcnj2b': '#DD8452'}
WINDOW_KB = 200
ROW = 1.3


def main():
    os.makedirs('figures', exist_ok=True)
    table = collections.defaultdict(list)
    for r in csv.DictReader(open(SRC), delimiter='\t'):
        table[(r['species'], r['focal'])].append(r)

    tracks = []
    for species, region, key, focals in BLOCKS:
        genes, contig = {}, None
        for focal in focals:
            for r in table[(key, focal)]:
                contig = r['contig']
                genes[(r['neighbor'], r['n_start'])] = (
                    PROPOSED_NAMES.get(r['neighbor'], r['neighbor']),
                    int(r['n_start']), int(r['n_end']), r['strand'])
        tracks.append((species, region, contig, sorted(genes.values(), key=lambda g: g[1])))

    import figstyle
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyArrow

    fig, ax = plt.subplots(figsize=(15, 2.7 * len(tracks) + 1.6))
    scale = 0.9 / (WINDOW_KB * 500)
    placed = []
    for i, (species, region, contig, genes) in enumerate(tracks):
        y = -i * ROW
        focal = [g for g in genes if g[0] in FOCAL_COLORS]
        centre = (min(g[1] for g in focal) + max(g[2] for g in focal)) / 2
        k16 = next((g for g in focal if g[0].startswith('kcnj16')), None)
        k2 = next((g for g in focal if g[0].startswith('kcnj2')), None)
        mirrored = bool(k16 and k2 and k16[1] > k2[1])

        ax.plot([-0.95, 0.95], [y, y], color='#d8d8d8', lw=2.0, solid_capstyle='round', zorder=0)
        positions, neighbour_labels = {}, []
        for name, start, end, strand in genes:
            x0, x1 = (start - centre) * scale, (end - centre) * scale
            if mirrored:
                x0, x1, strand = -x1, -x0, '-' if strand == '+' else '+'
            if x1 < -0.93 or x0 > 0.93:
                continue
            is_focal = name in FOCAL_COLORS
            width = max(x1 - x0, 0.016)
            ax.add_patch(FancyArrow(
                x0 if strand == '+' else x1, y, width if strand == '+' else -width, 0,
                width=.17 if is_focal else .085, head_width=.26 if is_focal else .13,
                head_length=min(width * .45, .03), length_includes_head=True, zorder=3,
                facecolor=FOCAL_COLORS.get(name, '#c2c2c2'),
                linewidth=1.0 if is_focal else 0, edgecolor='black' if is_focal else 'none'))
            cx = (x0 + x1) / 2
            if is_focal:
                level = sum(1 for px in positions.values() if abs(px - cx) < 0.18)
                label_y = y + .42 + level * .30
                ax.plot([cx, cx], [y + .13, label_y - .06], color='#999999', lw=.7, zorder=2)
                ax.text(cx, label_y, name, ha='center', fontsize=12.5, fontstyle='italic',
                        fontweight='bold', zorder=4)
                positions[name] = cx
            elif width > 0.030 and not name.startswith(('LOC', 'si:', 'zgc')):
                if any(abs(cx - px) < 0.055 for px, _ in neighbour_labels):
                    continue
                row = 1 if any(abs(cx - px) < 0.13 and r == 0 for px, r in neighbour_labels) else 0
                ax.plot([cx, cx], [y - .13, y - (.30 + row * .22)], color='#cccccc', lw=.6)
                ax.text(cx, y - (.34 + row * .22), name, ha='center', fontsize=8,
                        color='#8a8a8a', fontstyle='italic', zorder=4)
                neighbour_labels.append((cx, row))
        placed.append((y, positions))

        ax.text(-1.03, y + .12, species, ha='right', va='center', fontsize=12, fontstyle='italic')
        note = '  [shown reverse-complemented]' if mirrored else ''
        ax.text(-1.03, y - .22, f'{region}  ({contig}){note}', ha='right', va='center',
                fontsize=8.5, color='#666666')
        if any(n.startswith('kcnj2') for n in positions) and \
                not any(n.startswith('kcnj16') for n in positions):
            ax.text(-1.03, y - .46, 'no $\\it{kcnj16}$ copy in this block', ha='right',
                    va='center', fontsize=10, color='#C44E52')

    for (y1, p1), (y2, p2) in zip(placed, placed[1:]):
        for name in set(p1) & set(p2):
            ax.plot([p1[name], p2[name]], [y1 - .20, y2 + .20], color='#4C72B0', lw=1.2,
                    ls='--', alpha=.55, zorder=1)

    bar_y = -ROW * (len(tracks) - 1) - .95
    ax.plot([0.55, 0.55 + 50_000 * scale], [bar_y] * 2, color='k', lw=2.4)
    ax.text(0.55 + 25_000 * scale, bar_y - .22, '50 kb', ha='center', fontsize=10)
    ax.legend(handles=[
        Line2D([], [], color='#C44E52', lw=8, label='$\\it{kcnj16}$ genes'),
        Line2D([], [], color='#DD8452', lw=8, label='$\\it{kcnj2}$ genes'),
        Line2D([], [], color='#c2c2c2', lw=8, label='flanking genes'),
        Line2D([], [], color='#4C72B0', lw=1.5, ls='--', label='orthologous pair')],
        fontsize=10, loc='lower left', frameon=False, ncol=4, bbox_to_anchor=(-0.02, -0.02))
    figstyle.title(ax, 'Conservation of the $\\it{kcnj16}$–$\\it{kcnj2}$ gene pair',
                   fontsize=13, pad=14)
    ax.set_xlim(-1.45, 1.02)
    ax.set_ylim(-ROW * len(tracks) - .35, 1.35)
    ax.axis('off')
    fig.tight_layout()
    fig.savefig('figures/fig4_synteny.png')
    print('figures/fig4_synteny.png')


if __name__ == '__main__':
    main()
