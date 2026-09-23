#!/usr/bin/env python3
"""Fig. 6 and Fig. S6: the 169-sequence tree (six fish species and human).

Support values belong to branches. Bio.Phylo does not move node labels when a tree is
re-rooted, so they are read from the bipartitions of the unrooted tree before rooting.
"""
import os

from Bio import Phylo

from common import PROPOSED_NAMES

TREE = 'results/07_taxa/kir_all.treefile'
OUT_FOCAL = 'figures/fig3_focal_clades.png'
OUT_FULL = 'figures/figS6_extended_tree.png'

RENAME = {f'Bs_{loc}': f'Bs_{name}' for loc, name in PROPOSED_NAMES.items()}
# trimAl cuts sequence names at ':'; full gene symbols for display
DISPLAY = {'Am_si': 'Am_si:ch211-113j13.2', 'Am_zgc': 'Am_zgc:162160'}
COLORS = {'Bs': '#C44E52', 'At': '#DD8452', 'Ol': '#937860',
          'Dr': '#4C72B0', 'Am': '#64B5CD', 'Lo': '#55A868', 'Hs': '#6e6e6e'}
SPECIES = {'Bs': 'Betta splendens', 'At': 'Anabas testudineus', 'Ol': 'Oryzias latipes',
           'Dr': 'Danio rerio', 'Am': 'Astyanax mexicanus', 'Lo': 'Lepisosteus oculatus',
           'Hs': 'Homo sapiens'}
ORDER = ['Bs', 'At', 'Ol', 'Dr', 'Am', 'Lo', 'Hs']

FOCAL = [
    ('kcnj16', ['Bs_kcnj16a', 'Bs_LOC114860957', 'Hs_KCNJ16']),
    ('kcnj4', ['Bs_LOC114863156', 'Dr_kcnj4', 'Hs_KCNJ4']),
    ('kcnj10', ['Bs_kcnj10a', 'Bs_LOC114844899', 'Dr_kcnj10b', 'Hs_KCNJ10']),
]


def split_supports(tree):
    """Bipartition -> support label, taken from the unrooted tree."""
    tips = frozenset(t.name for t in tree.get_terminals())
    table = {}
    for cl in tree.get_nonterminals():
        s = frozenset(t.name for t in cl.get_terminals())
        lab = cl.name if cl.name else cl.confidence
        if lab is not None and 1 < len(s) < len(tips) - 1:
            table[s] = table[tips - s] = str(lab)
    return table


def support(clade, table):
    """(SH-aLRT, UFBoot) of the branch above this clade, looked up by its bipartition."""
    lab = table.get(frozenset(t.name for t in clade.get_terminals()))
    if lab is None:
        return None, None
    try:
        if '/' in lab:
            a, b = lab.split('/')
            return float(a), float(b)
        return None, float(lab)
    except ValueError:
        return None, None


def layout(root):
    leaves = root.get_terminals()
    ypos = {id(l): i for i, l in enumerate(leaves)}

    def y_of(cl):
        if not cl.is_terminal():
            ypos[id(cl)] = sum(y_of(c) for c in cl) / len(cl)
        return ypos[id(cl)]
    y_of(root)

    xpos = {}

    def x_of(cl, x0):
        xpos[id(cl)] = x0
        for c in cl:
            x_of(c, x0 + (c.branch_length or 0.0))
    x_of(root, 0.0)
    return leaves, xpos, ypos


def draw_tree(ax, fig, root, numbers, fontsize, table):
    from matplotlib.transforms import offset_copy
    leaves, xpos, ypos = layout(root)
    star_tr = offset_copy(ax.transData, fig=fig, x=-7, y=0, units='points')

    def draw(cl):
        if cl.is_terminal():
            return
        x, y = xpos[id(cl)], ypos[id(cl)]
        ys = [ypos[id(c)] for c in cl]
        ax.plot([x, x], [min(ys), max(ys)], color='#3a3a3a', lw=1.1, zorder=1)
        for c in cl:
            ax.plot([x, xpos[id(c)]], [ypos[id(c)]] * 2, color='#3a3a3a', lw=1.1, zorder=1)
            draw(c)
        sh, uf = support(cl, table)
        if uf is None or cl is root:
            return
        if numbers:
            txt = f'{sh:g}/{uf:g}' if sh is not None else f'{uf:g}'
            ax.annotate(txt, (x, y), xytext=(3, 4), textcoords='offset points',
                        ha='left', va='bottom', fontsize=fontsize - 2, color='#444444',
                        bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=.85),
                        zorder=5)
        elif uf >= 95:
            ax.plot(x, y, 'o', ms=3.4, color='#1a1a1a', zorder=4)
        elif uf >= 70:
            ax.plot(x, y, 'o', ms=3.4, mfc='white', mec='#1a1a1a', mew=.9, zorder=4)
    draw(root)

    xmax = max(xpos[id(l)] for l in leaves)
    label_x = xmax * 1.03
    for l in leaves:
        x, y = xpos[id(l)], ypos[id(l)]
        name = RENAME.get(l.name, DISPLAY.get(l.name, l.name))
        tag, gene = name.split('_', 1)
        col = COLORS.get(tag, 'k')
        bold = 'bold' if tag == 'Bs' else 'normal'
        ax.plot([x, label_x - xmax * .005], [y, y], color='#d0d0d0', lw=.5, ls=':', zorder=0)
        if l.name in RENAME:
            ax.plot(label_x, y, marker='*', ms=fontsize - 1, color=col, transform=star_tr,
                    clip_on=False, zorder=4)
        ax.annotate(tag, (label_x, y), fontsize=fontsize, va='center', color=col,
                    fontweight=bold)
        ax.annotate(gene, (label_x, y), xytext=(fontsize * 2.2, 0), textcoords='offset points',
                    fontsize=fontsize, va='center', color=col, fontweight=bold,
                    fontstyle='italic')
    ax.set_ylim(len(leaves) - .4, -.8)
    ax.set_xlim(-xmax * .03, xmax * 1.55)
    ax.axis('off')
    return xmax


def legend(ax, numbers, **kw):
    from matplotlib.lines import Line2D

    def italic(text):
        return '$\\it{' + text.replace(' ', '\\ ') + '}$'
    handles = [Line2D([], [], color=COLORS[k], lw=4, label=f'{italic(SPECIES[k])} ({k})')
               for k in ORDER]
    if not numbers:
        handles += [Line2D([], [], marker='o', ls='', ms=6, color='#1a1a1a',
                           label='UFBoot ≥ 95'),
                    Line2D([], [], marker='o', ls='', ms=6, mfc='white', mec='#1a1a1a',
                           label='UFBoot 70–94')]
    handles.append(Line2D([], [], marker='*', ls='', ms=10, color=COLORS['Bs'],
                          label='named in this study'))
    ax.legend(handles=handles, frameon=False, **kw)


def scale_bar(ax, y, length):
    ax.plot([0, length], [y, y], color='k', lw=1.8)
    ax.annotate(f'{length} substitutions/site', (length / 2, y), xytext=(0, 4),
                textcoords='offset points', ha='center', va='bottom', fontsize=8)


def main():
    os.makedirs('figures', exist_ok=True)
    import figstyle
    import matplotlib.pyplot as plt

    t = Phylo.read(TREE, 'newick')
    table = split_supports(t)
    t.root_at_midpoint()
    t.ladderize()

    # --- focal clades ---
    clades = [(title, t.common_ancestor(names)) for title, names in FOCAL]
    sizes = [len(c.get_terminals()) for _, c in clades]
    fig, axes = plt.subplots(len(clades), 1, figsize=(8.2, 0.34 * sum(sizes) + 2.2),
                             gridspec_kw={'height_ratios': sizes})
    for ax, (title, clade), letter in zip(axes, clades, 'abc'):
        draw_tree(ax, fig, clade, numbers=True, fontsize=9, table=table)
        ax.text(-0.02, 1.0, letter, transform=ax.transAxes, fontsize=13, fontweight='bold',
                va='top', ha='right')
        scale_bar(ax, len(clade.get_terminals()) - .7, 0.2)
        figstyle.title(ax, title, fontsize=11, loc='left')
    legend(axes[-1], numbers=True, fontsize=8.5, loc='upper center',
           bbox_to_anchor=(0.5, -0.04), ncol=4, handletextpad=.5, columnspacing=1.4)
    fig.tight_layout()
    figstyle.save(fig, OUT_FOCAL)
    print(f'{OUT_FOCAL} written ({sizes} taxa per panel)')

    # --- whole tree ---
    # the two halves of the midpoint-rooted tree side by side, same leaf spacing
    halves = list(t.root.clades)
    sizes = [len(c.get_terminals()) for c in halves]
    fig, axes = plt.subplots(1, len(halves), figsize=(7.2 * len(halves), 0.135 * max(sizes) + 2.4))
    for ax, clade, letter in zip(axes, halves, 'ab'):
        draw_tree(ax, fig, clade, numbers=False, fontsize=6.5, table=table)
        ax.set_ylim(max(sizes) - .4, -.8)
        ax.text(0.0, 1.0, letter, transform=ax.transAxes, fontsize=13, fontweight='bold',
                va='top', ha='left')
        scale_bar(ax, len(clade.get_terminals()) + .6, 0.5)
    small = axes[sizes.index(min(sizes))]
    legend(small, numbers=False, fontsize=8.5, loc='lower left', bbox_to_anchor=(0.0, 0.02),
           ncol=2, handletextpad=.5, columnspacing=1.4)
    fig.tight_layout()
    figstyle.save(fig, OUT_FULL)
    print(f'{OUT_FULL} written ({sum(sizes)} taxa in two panels)')


if __name__ == '__main__':
    main()
