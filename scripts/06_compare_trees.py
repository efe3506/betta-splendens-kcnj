#!/usr/bin/env python3
"""Robustness check: is every B. splendens sequence sister to the same sequences in two trees?."""
import sys

from Bio import Phylo


def sisters(path):
    tree = Phylo.read(path, 'newick')
    tree.root_at_midpoint()
    parent = {}
    for clade in tree.find_clades(order='level'):
        for child in clade:
            parent[child] = clade
    out = {}
    for leaf in tree.get_terminals():
        if not leaf.name.startswith('Bs_') or leaf not in parent:
            continue
        node = parent[leaf]
        group = []
        for child in node:
            if child is not leaf:
                group += [t.name for t in child.get_terminals()]
        group.sort()
        support = node.name or (str(node.confidence) if node.confidence is not None else '-')
        out[leaf.name] = (','.join(group), support)
    return out


def main(tree_a, tree_b, out_path):
    a, b = sisters(tree_a), sisters(tree_b)
    same = 0
    with open(out_path, 'w') as fh:
        fh.write('sequence\tsister_tree1\tsupport_tree1\tsister_tree2\tsupport_tree2\tagree\n')
        for name in sorted(a):
            s1, p1 = a[name]
            s2, p2 = b.get(name, ('-', '-'))
            agree = s1 == s2
            same += agree
            fh.write(f'{name}\t{s1}\t{p1}\t{s2}\t{p2}\t{"yes" if agree else "NO"}\n')
    print(f'{same} of {len(a)} B. splendens sequences have the same sister group in both trees')


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit('usage: 06_compare_trees.py <tree1> <tree2> <out.tsv>')
    main(*sys.argv[1:])
