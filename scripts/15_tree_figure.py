#!/usr/bin/env python3
"""Fig. 2: maximum-likelihood tree of the Kir (kcnj) family."""
import os

from Bio import Phylo

from common import PROPOSED_NAMES, subfamily

TREE = 'results/03_phylogeny/kcnj.treefile'
COLORS = {'Bs': '#C44E52', 'Dr': '#4C72B0', 'Hs': '#6e6e6e'}
SPECIES = {'Bs': 'Betta splendens', 'Dr': 'Danio rerio', 'Hs': 'Homo sapiens'}


def split_supports(tree):
    """Bipartition -> support label, taken from the unrooted tree.

    Support values belong to branches. Bio.Phylo does not move node labels when a tree is
    re-rooted, so they are looked up by bipartition instead of being read from the nodes.
    """
    tips = frozenset(t.name for t in tree.get_terminals())
    table = {}
    for clade in tree.get_nonterminals():
        members = frozenset(t.name for t in clade.get_terminals())
        label = clade.name if clade.name else clade.confidence
        if label is not None and 1 < len(members) < len(tips) - 1:
            table[members] = table[tips - members] = str(label)
    return table


def ufboot(clade, table):
    label = table.get(frozenset(t.name for t in clade.get_terminals()), '')
    try:
        return float(label.split('/')[1])
    except (IndexError, ValueError):
        return None


def main():
    os.makedirs('figures', exist_ok=True)
    tree = Phylo.read(TREE, 'newick')
    table = split_supports(tree)
    tree.root_at_midpoint()
    tree.ladderize()
    leaves = tree.get_terminals()

    y = {id(leaf): i for i, leaf in enumerate(leaves)}

    def place(clade):
        if not clade.is_terminal():
            y[id(clade)] = sum(place(c) for c in clade) / len(clade)
        return y[id(clade)]
    place(tree.root)

    x = {}

    def depth(clade, x0=0.0):
        x[id(clade)] = x0 + (clade.branch_length or 0.0)
        for child in clade:
            depth(child, x[id(clade)])
    depth(tree.root)

    import figstyle
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.transforms import offset_copy

    fig, ax = plt.subplots(figsize=(10.5, 16))
    star_position = offset_copy(ax.transData, fig=fig, x=-7, y=0, units='points')

    def draw(clade):
        if clade.is_terminal():
            return
        cx, cy = x[id(clade)], y[id(clade)]
        ys = [y[id(c)] for c in clade]
        ax.plot([cx, cx], [min(ys), max(ys)], color='#3a3a3a', lw=1.1, zorder=1)
        for child in clade:
            ax.plot([cx, x[id(child)]], [y[id(child)]] * 2, color='#3a3a3a', lw=1.1, zorder=1)
            draw(child)
        support = ufboot(clade, table)
        if support is not None and clade is not tree.root:
            if support >= 95:
                ax.plot(cx, cy, 'o', ms=4.2, color='#1a1a1a', zorder=4)
            elif support >= 70:
                ax.plot(cx, cy, 'o', ms=4.2, mfc='white', mec='#1a1a1a', mew=1.0, zorder=4)
    draw(tree.root)

    x_max = max(x.values())
    label_x = x_max * 1.03
    groups = {}
    for leaf in leaves:
        tag, gene = leaf.name.split('_', 1)
        named_here = gene in PROPOSED_NAMES
        gene = PROPOSED_NAMES.get(gene, gene)
        colour = COLORS.get(tag, 'k')
        weight = 'bold' if tag == 'Bs' else 'normal'
        ly = y[id(leaf)]
        ax.plot([x[id(leaf)], label_x - x_max * 0.005], [ly, ly], color='#d0d0d0',
                lw=.5, ls=':', zorder=0)
        if named_here:
            ax.plot(label_x, ly, marker='*', ms=8, color=colour, transform=star_position,
                    clip_on=False, zorder=4)
        ax.annotate(tag, (label_x, ly), xytext=(0, 0), textcoords='offset points',
                    fontsize=8.6, va='center', color=colour, fontweight=weight)
        ax.annotate(gene, (label_x, ly), xytext=(19, 0), textcoords='offset points',
                    fontsize=8.6, va='center', color=colour, fontweight=weight,
                    fontstyle='italic')
        groups.setdefault(subfamily(gene), []).append(ly)

    bracket_x = x_max * 1.40
    for name, ys in groups.items():
        ys = sorted(ys)
        if name is None or ys[-1] - ys[0] + 1 != len(ys):
            continue
        ax.plot([bracket_x] * 2, [ys[0] - .3, ys[-1] + .3], color='#555555', lw=2.2)
        for end in (ys[0] - .3, ys[-1] + .3):
            ax.plot([bracket_x, bracket_x - x_max * .012], [end, end], color='#555555', lw=2.2)
        ax.text(bracket_x + x_max * .02, (ys[0] + ys[-1]) / 2, name, fontsize=11,
                va='center', fontweight='bold', color='#333333')

    ax.plot([0, 0.5], [-1.0, -1.0], color='k', lw=2)
    ax.text(0.25, -1.55, '0.5 substitutions/site', ha='center', fontsize=9)

    handles = [Line2D([], [], color=COLORS[k], lw=4,
                      label='$\\it{' + SPECIES[k].replace(' ', '\\ ') + '}$ (' + k + ')')
               for k in ('Bs', 'Dr', 'Hs')]
    handles += [Line2D([], [], marker='o', ls='', ms=7, color='#1a1a1a', label='UFBoot ≥ 95'),
                Line2D([], [], marker='o', ls='', ms=7, mfc='white', mec='#1a1a1a',
                       label='UFBoot 70–94'),
                Line2D([], [], marker='*', ls='', ms=11, color=COLORS['Bs'],
                       label='named in this study')]
    ax.legend(handles=handles, fontsize=10, loc='upper center', bbox_to_anchor=(0.45, 0.0),
              ncol=3, frameon=False, handletextpad=.6, columnspacing=2.0)
    figstyle.title(ax, 'Maximum-likelihood phylogeny of the Kir (kcnj) gene family',
                   fontsize=13, pad=16)
    ax.set_ylim(-1.8, len(leaves) + 0.5)
    ax.set_xlim(-x_max * .03, x_max * 1.62)
    ax.invert_yaxis()
    ax.axis('off')
    fig.tight_layout()
    fig.savefig('figures/fig2_phylogeny.png')
    print(f'figures/fig2_phylogeny.png ({len(leaves)} taxa)')


if __name__ == '__main__':
    main()
