#!/usr/bin/env python3
"""Codon alignments and labelled trees for the tests of selection (HyPhy RELAX and aBSREL).

For each family the teleost and gar members of its clade in the 169-sequence tree are used.
Codons are threaded onto a MAFFT protein alignment. Branches of the "b" clade are labelled
{test}, those of the "a" clade {reference}; the gar branch is left unlabelled. A second
alignment without any gapped column is written for the sensitivity analysis.
"""
import os
import re
import subprocess

from Bio import Phylo
from Bio.Seq import Seq

TREE = 'results/07_taxa/kir_all.treefile'
OUTDIR = 'results/10_selection'
FAMILIES = {
    # family: (tips that define the clade, tips of the "b" clade)
    'kcnj16': (['Bs_kcnj16a', 'Bs_LOC114860957', 'Lo_kcnj16a'],
               ['Bs_LOC114860957', 'At_LOC113160616', 'Ol_LOC101171333']),
    'kcnj4': (['Bs_LOC114863156', 'Dr_kcnj4', 'Lo_LOC107078893'],
              ['Bs_LOC114863156', 'At_LOC113162111', 'Ol_LOC101162982']),
}


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith('>'):
            name = line[1:].rstrip('\n')
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: ''.join(v) for k, v in seqs.items()}


def label_to_protein():
    """Tree label -> RefSeq protein accession."""
    m = {}
    for head in read_fasta('results/01_identification/REF_kcnj_canonical.faa'):
        gene, pid = head.split()[0].split('|')
        m[f'Bs_{gene}'] = pid
    for tag in ('At', 'Ol', 'Am', 'Lo', 'Dr'):
        path = f'results/07_taxa/{tag}_kir_genes.tsv'
        for line in list(open(path))[1:]:
            f = line.rstrip('\n').split('\t')
            m[f'{tag}_{f[0]}'] = f[1]
    return m


def cds_by_protein(tag):
    """cds_from_genomic.fna carries [protein_id=...] in the header."""
    out = {}
    for head, seq in read_fasta(f'data/species/{tag}/cds.fna').items():
        pid = re.search(r'\[protein_id=([^\]]+)\]', head)
        if pid:
            out[pid.group(1)] = seq.upper()
    return out


def newick(clade, test, ref):
    if clade.is_terminal():
        tag = '{test}' if clade.name in test else '{reference}' if clade.name in ref else ''
        return clade.name + tag
    tips = {t.name for t in clade.get_terminals()}
    tag = '{test}' if tips <= test else '{reference}' if tips <= ref else ''
    return '(' + ','.join(newick(c, test, ref) for c in clade) + ')' + tag


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    l2p = label_to_protein()
    cds_cache = {}
    tree = Phylo.read(TREE, 'newick')

    for fam, (anchors, b_tips) in FAMILIES.items():
        clade = tree.common_ancestor(anchors)
        members = [t.name for t in clade.get_terminals() if not t.name.startswith('Hs_')]
        test = set(b_tips)
        ref = {m for m in members if m not in test and not m.startswith('Lo_')}

        prot, codon, rows = {}, {}, []
        for lab in members:
            tag, pid = lab[:2], l2p.get(lab)
            if tag not in cds_cache:
                cds_cache[tag] = cds_by_protein(tag)
            cds = cds_cache[tag].get(pid)
            if cds is None:
                raise SystemExit(f'{lab}: no CDS for {pid}')
            cds = cds[:len(cds) - len(cds) % 3]
            aa = str(Seq(cds).translate())
            if aa.endswith('*'):
                aa, cds = aa[:-1], cds[:-3]
            if '*' in aa:
                raise SystemExit(f'{lab}: internal stop codon')
            prot[lab], codon[lab] = aa, cds
            role = 'test' if lab in test else 'reference' if lab in ref else 'outgroup'
            rows.append((lab, pid, len(aa), role))

        pfa = f'{OUTDIR}/{fam}_protein.faa'
        with open(pfa, 'w') as fh:
            for lab, aa in prot.items():
                fh.write(f'>{lab}\n{aa}\n')
        aln = subprocess.run(['mafft', '--localpair', '--maxiterate', '1000', '--quiet', pfa],
                             check=True, text=True, capture_output=True).stdout
        paln = {}
        for block in aln.strip().split('>')[1:]:
            head, *seq = block.split('\n')
            paln[head.strip()] = ''.join(seq)

        caln = {}
        for lab, a in paln.items():
            it = iter(range(0, len(codon[lab]), 3))
            caln[lab] = ['---' if ch == '-' else codon[lab][next(it):][:3] for ch in a]
        with open(f'{OUTDIR}/{fam}_codon.fna', 'w') as fh:
            for lab, cols in caln.items():
                fh.write(f'>{lab}\n{"".join(cols)}\n')
        # sensitivity analysis: drop every column that contains a gap
        keep = [j for j in range(len(next(iter(caln.values()))))
                if all(cols[j] != '---' for cols in caln.values())]
        with open(f'{OUTDIR}/{fam}_codon_nogap.fna', 'w') as fh:
            for lab, cols in caln.items():
                fh.write(f'>{lab}\n{"".join(cols[j] for j in keep)}\n')

        sub = Phylo.read(TREE, 'newick').common_ancestor(anchors)
        for t in [t for t in sub.get_terminals() if t.name not in members]:
            sub.prune(t)
        with open(f'{OUTDIR}/{fam}_tagged.nwk', 'w') as fh:
            fh.write(newick(sub, test, ref) + ';\n')
        with open(f'{OUTDIR}/{fam}_members.tsv', 'w') as fh:
            fh.write('label\tprotein\tlength_aa\trole\n')
            for r in rows:
                fh.write('\t'.join(map(str, r)) + '\n')

        ncol = len(next(iter(paln.values())))
        print(f'{fam}: {len(members)} sequences, {ncol} codon columns ({len(keep)} without gaps); '
              f'test={len(test)} reference={len(ref)}')


if __name__ == '__main__':
    main()
