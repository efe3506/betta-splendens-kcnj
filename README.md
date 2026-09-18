# The Kir (*kcnj*) gene family of *Betta splendens*

Code and reference results for a genome-wide survey of the inwardly rectifying potassium
channel (Kir, *kcnj*) gene family in the Siamese fighting fish, *Betta splendens*, using
four public genome assemblies.

The whole analysis is one bash script, `run_pipeline.sh`, which calls a set of small Python
scripts in `scripts/`. It needs no input other than an internet connection: every sequence
is downloaded from public databases.

## What the analysis does

| Step | Question | Main tools |
|---|---|---|
| 1 | Download the four assemblies and the RefSeq annotation | NCBI Datasets |
| 2 | How complete is each assembly? | seqkit, compleasm |
| 3 | Build the query set (human + zebrafish Kir proteins) and the Pfam profiles | UniProt, InterPro |
| 4 | Which annotated *B. splendens* proteins are Kir channels? | HMMER, BLAST |
| 5 | Is every gene present and intact in all four assemblies? | miniprot |
| 6 | Whose orthologue is each gene? | MAFFT, trimAl, IQ-TREE |
| 7 | Are the orthology calls supported by gene neighbourhood? Is there gene conversion? Where do the unnamed Kir4/Kir5 genes belong? | custom scripts, IQ-TREE |
| 8 | Exon structure and conserved motifs | MEME |
| 9 | A closer look at the *kcnj15* locus in the four assemblies | minimap2, MAFFT |
| 10 | Do predicted structures agree with the orthology calls? | AlphaFold DB models, PyMOL |
| 11 | Remaining figures | Matplotlib |

Assemblies used:

| Label | Accession | Note |
|---|---|---|
| REF | GCF_900634795.4 | RefSeq reference, annotated |
| NCU | GCA_024678985.1 | |
| LOEWE | GCA_013403625.1 | |
| BGI | GCA_003650155.1 | |

## Requirements

Linux, [conda](https://docs.conda.io) (or mamba), about 10 GB of free disk space and an
internet connection. Three environments are used because the tools cannot all be installed
together:

```bash
conda env create -f envs/kcnj.yml         # main environment
conda env create -f envs/kcnj-qc.yml      # compleasm only
conda env create -f envs/kcnj-pymol.yml   # PyMOL only
```

All versions are pinned in the environment files.

## Running

```bash
conda activate kcnj
bash run_pipeline.sh
```

The script is called from the repository root and only uses relative paths. Tables are
written to `results/`, figures to `figures/`, downloads to `data/`. A step whose output
already exists is skipped, so an interrupted run can simply be started again.

Settings are environment variables:

| Variable | Default | Meaning |
|---|---|---|
| `THREADS` | 8 | CPU threads |
| `RUN_QC` | 1 | `0` skips compleasm (about one hour per genome) |
| `REFRESH_INPUTS` | 0 | `1` rebuilds the query set and the HMMs from the live databases instead of using `reference_data/` |
| `FIG_TITLES` | 0 | `1` draws a title inside each figure (journals usually do not want one) |

Example: `THREADS=16 RUN_QC=0 bash run_pipeline.sh`

With `RUN_QC=0` most of the running time is spent in miniprot and in the four IQ-TREE runs.

## Repository layout

```
run_pipeline.sh      the whole analysis, top to bottom
scripts/             Python scripts called by the pipeline (numbered in order of use)
envs/                conda environments with pinned versions
reference_data/      frozen inputs: query proteins and Pfam HMMs as used in the paper
expected_results/    the tables and trees the pipeline is expected to produce
```

## Scripts

| Script | Purpose | Main output |
|---|---|---|
| `01_filter_queries.py` | one protein per gene, 150–550 aa, from a raw UniProt download | `data/queries/kcnj_queries.faa` |
| `02_map_proteins_to_genes.py` | HMM protein hits → genes; longest isoform per gene | `REF_kcnj_genes.tsv`, `REF_kcnj_canonical.faa` |
| `03_loci_matrix.py` | miniprot alignments → loci; presence, identity and frameshift tables | `presence_matrix.tsv`, `identity_matrix.tsv` |
| `04_pair_distance.py` | distance and relative orientation of two genes in each assembly | `pair_kcnj15_kcnj6.tsv` |
| `05_label_sequences.py` | short labels (`Bs_`, `Dr_`, `Hs_`) for alignment and trees | `kcnj_labeled.faa` |
| `06_compare_trees.py` | sister group of every *B. splendens* gene in two trees | `sister_comparison.tsv` |
| `07_synteny.py` | ten annotated neighbours on each side of focal genes, shared between species | `synteny_scores.tsv` |
| `08_gene_conversion.py` | sliding-window identity in the *kcnj10* group | `gene_conversion_windows.tsv`, Fig. S2 |
| `09_expand_kcnj10.py` | Kir4/Kir5 sequence sets with *Anabas* and medaka blastp hits | `kir45.faa`, `kir45_expanded.faa` |
| `10_gene_structure.py` | exon counts and lengths; chromosome map and exon diagrams | `gene_structure.tsv`, Fig. 1, Fig. 3 |
| `11_motif_matrix.py` | motif presence/absence from the MEME output | `motif_presence.tsv` |
| `12_kcnj15_locus.py` | the *kcnj15* region in the four assemblies (`utr`, `lead`, `cds`) | `variant_genotypes_*.tsv`, `cds_translation.tsv` |
| `13_structure_analysis.py` | pLDDT profiles and RMSD of AlphaFold DB models | `plddt_summary.tsv`, `structural_comparison.tsv`, Fig. S3 |
| `14_structure_figures.pml` | PyMOL renderings | Fig. 5, Fig. S4, Fig. S5 |
| `15_tree_figure.py` | family tree | Fig. 2 |
| `16_synteny_figure.py` | *kcnj16*–*kcnj2* blocks in two species | Fig. 4 |

## Reproducibility notes

- **Frozen inputs.** UniProt and InterPro are updated continuously. The query proteins and
  the HMM profiles used in the paper are therefore shipped in `reference_data/`. With
  `REFRESH_INPUTS=1` the same procedure is run against the current databases, and small
  differences in the query set are to be expected.
- **Checksums.** Every NCBI Datasets package is verified against the md5 sums it contains.
- **Random seeds.** IQ-TREE is run with the seeds of the published analysis, which are
  written in `run_pipeline.sh`. Likelihoods and topology are reproducible for a given
  IQ-TREE version and thread count; support values can differ in the last digit between
  machines.
- **Comparing with the published results.** Tab-separated tables should be identical:

  ```bash
  diff results/02_cross_assembly/presence_matrix.tsv expected_results/cross_assembly/presence_matrix.tsv
  ```

  Trees are best compared by topology rather than by text.
- **Predicted structures** are downloaded from the AlphaFold Protein Structure Database.
  The paper used model version 6. The database may release new versions; the version that
  is downloaded is printed by `scripts/13_structure_analysis.py`.
- **Annotation releases.** The zebrafish, climbing perch and medaka annotations are fetched
  by assembly accession. If NCBI replaces an annotation release, gene symbols in the
  synteny table can change.

## Data sources and licences

Genome assemblies and annotations: NCBI GenBank/RefSeq. Protein sequences: UniProt
(CC BY 4.0). Domain profiles: Pfam via InterPro (CC0). Predicted structures: AlphaFold
Protein Structure Database (CC BY 4.0). The code in this repository is released under the
MIT licence (see `LICENSE`).

## Citation

If you use this code, please cite the accompanying article (reference to be added on
publication).
