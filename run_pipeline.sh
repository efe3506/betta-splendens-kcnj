#!/usr/bin/env bash
set -euo pipefail

THREADS=${THREADS:-8}
RUN_QC=${RUN_QC:-1}
REFRESH_INPUTS=${REFRESH_INPUTS:-0}
QC_ENV=${QC_ENV:-kcnj-qc}
PYMOL_ENV=${PYMOL_ENV:-kcnj-pymol}

ASSEMBLIES="REF NCU LOEWE BGI"
declare -A ACCESSION
ACCESSION[REF]=GCF_900634795.4
ACCESSION[NCU]=GCA_024678985.1
ACCESSION[LOEWE]=GCA_013403625.1
ACCESSION[BGI]=GCA_003650155.1

DANIO=GCF_052040795.1
ANABAS=GCF_900324465.3
MEDAKA=GCF_053564925.1

ID=results/01_identification
CA=results/02_cross_assembly
PH=results/03_phylogeny
ST=results/04_structure

step() {
  printf '\n=== %s ===\n' "$1"
}

fetch() {
  datasets download genome accession "$1" --include "$2" --filename "data/tmp/$1.zip" --no-progressbar
  unzip -q -o "data/tmp/$1.zip" -d "data/tmp/$1"
  (cd "data/tmp/$1" && md5sum -c --quiet md5sum.txt)
}

build_tree() {
  iqtree3 -s "$1" -m MFP -B 1000 --alrt 1000 -T "$THREADS" --seed "$2" --prefix "$3"
}

mkdir -p data/tmp data/genomes data/annotation data/queries figures
mkdir -p results/00_qc $ID $CA $PH $ST results/05_kcnj15_locus results/06_protein_model


step "1. Genome assemblies"
for g in $ASSEMBLIES; do
  if [[ -s data/genomes/$g.fna ]]; then
    continue
  fi
  acc=${ACCESSION[$g]}
  src=data/tmp/$acc/ncbi_dataset/data/$acc
  if [[ $g == REF ]]; then
    fetch "$acc" genome,protein,gff3
    mv "$src"/protein.faa data/annotation/REF_protein.faa
    mv "$src"/genomic.gff data/annotation/REF_genomic.gff
  else
    fetch "$acc" genome
  fi
  mv "$src"/"${acc}"_*_genomic.fna data/genomes/$g.fna
  cp data/tmp/$acc/ncbi_dataset/data/assembly_data_report.jsonl results/00_qc/${g}_assembly_report.jsonl
  samtools faidx data/genomes/$g.fna
done


step "2. Assembly statistics and gene-space completeness"
if [[ ! -s results/00_qc/genome_stats.tsv ]]; then
  seqkit stats -a -T -j "$THREADS" data/genomes/*.fna > results/00_qc/genome_stats.tsv
fi
if [[ $RUN_QC == 1 ]]; then
  if [[ ! -e data/busco_lineages/actinopterygii_odb12.done ]]; then
    conda run --no-capture-output -n "$QC_ENV" \
      compleasm download actinopterygii --odb odb12 -L data/busco_lineages
  fi
  for g in $ASSEMBLIES; do
    if [[ ! -s results/00_qc/compleasm_$g/summary.txt ]]; then
      conda run --no-capture-output -n "$QC_ENV" \
        compleasm run -a data/genomes/$g.fna -o results/00_qc/compleasm_$g \
          -l actinopterygii --odb odb12 -L data/busco_lineages -t "$THREADS"
    fi
  done
fi


step "3. Query proteins and Pfam profiles"
if [[ $REFRESH_INPUTS == 1 ]]; then
  curl -sf 'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28gene%3Akcnj%2A%29%20AND%20%28organism_id%3A9606%20OR%20organism_id%3A7955%29' \
    -o data/queries/uniprot_kcnj_raw.faa
  python3 scripts/01_filter_queries.py data/queries/uniprot_kcnj_raw.faa data/queries/kcnj_queries.faa
  rm -f data/queries/IRK_domains.hmm
  for pf in PF01007 PF08466 PF17655; do
    curl -sf "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/${pf}?annotation=hmm" \
      | gunzip >> data/queries/IRK_domains.hmm
  done
else
  cp reference_data/kcnj_queries.faa reference_data/IRK_domains.hmm data/queries/
fi


step "4. Identification in the reference proteome"
if [[ ! -s $ID/REF_IRK_domtbl.txt ]]; then
  hmmsearch --cut_ga --cpu "$THREADS" --domtblout $ID/REF_IRK_domtbl.txt \
    data/queries/IRK_domains.hmm data/annotation/REF_protein.faa > $ID/REF_IRK_hmmsearch.out
fi
python3 scripts/02_map_proteins_to_genes.py

makeblastdb -in data/queries/kcnj_queries.faa -dbtype prot -out $ID/kcnj_query_db > /dev/null
blastp -query $ID/REF_kcnj_canonical.faa -db $ID/kcnj_query_db -evalue 1e-20 \
    -num_threads "$THREADS" -max_target_seqs 5 \
    -outfmt '6 qseqid sseqid pident length evalue bitscore' \
  | sort -k1,1 -k6,6gr | awk '!seen[$1]++' > $ID/REF_best_hits.tsv


step "5. Cross-assembly validation"
cat $ID/REF_kcnj_canonical.faa data/queries/kcnj_queries.faa > $CA/all_queries.faa
for g in $ASSEMBLIES; do
  if [[ ! -s $CA/${g}_miniprot.gff ]]; then
    miniprot -t "$THREADS" --gff --outs=0.5 data/genomes/$g.fna $CA/all_queries.faa > $CA/${g}_miniprot.gff
  fi
done
python3 scripts/03_loci_matrix.py
python3 scripts/04_pair_distance.py kcnj15 kcnj6


step "6. Phylogeny"
python3 scripts/05_label_sequences.py
if [[ ! -s $PH/kcnj_aln.faa ]]; then
  mafft --localpair --maxiterate 1000 --thread "$THREADS" --quiet $PH/kcnj_labeled.faa > $PH/kcnj_aln.faa
fi
if [[ ! -s $PH/kcnj_trimmed.faa ]]; then
  trimal -in $PH/kcnj_aln.faa -out $PH/kcnj_trimmed.faa -automated1
fi
if [[ ! -s $PH/kcnj.treefile ]]; then
  build_tree $PH/kcnj_trimmed.faa 619787 $PH/kcnj
fi
if [[ ! -s $PH/kcnj_untrimmed.treefile ]]; then
  build_tree $PH/kcnj_aln.faa 104951 $PH/kcnj_untrimmed
fi
python3 scripts/06_compare_trees.py $PH/kcnj.treefile $PH/kcnj_untrimmed.treefile $PH/sister_comparison.tsv


step "7. Synteny, gene conversion test and the Kir4/Kir5 subfamily"
if [[ ! -s data/annotation/DANIO_genomic.gff ]]; then
  fetch $DANIO gff3
  mv data/tmp/$DANIO/ncbi_dataset/data/$DANIO/genomic.gff data/annotation/DANIO_genomic.gff
fi
if [[ ! -s data/annotation/ANABAS_protein.faa ]]; then
  fetch $ANABAS protein
  mv data/tmp/$ANABAS/ncbi_dataset/data/$ANABAS/protein.faa data/annotation/ANABAS_protein.faa
fi
if [[ ! -s data/annotation/MEDAKA_protein.faa ]]; then
  fetch $MEDAKA protein
  mv data/tmp/$MEDAKA/ncbi_dataset/data/$MEDAKA/protein.faa data/annotation/MEDAKA_protein.faa
fi
python3 scripts/07_synteny.py
python3 scripts/08_gene_conversion.py
python3 scripts/09_expand_kcnj10.py

if [[ ! -s $ST/kir45.treefile ]]; then
  mafft --localpair --maxiterate 1000 --quiet $ST/kir45.faa > $ST/kir45_aln.faa
  build_tree $ST/kir45_aln.faa 640120 $ST/kir45
fi
if [[ ! -s $ST/kir45_expanded.treefile ]]; then
  mafft --localpair --maxiterate 1000 --quiet $ST/kir45_expanded.faa > $ST/kir45_expanded_aln.faa
  build_tree $ST/kir45_expanded_aln.faa 31033 $ST/kir45_expanded
fi


step "8. Gene structure and conserved motifs"
python3 scripts/10_gene_structure.py
if [[ ! -s $ST/meme/meme.txt ]]; then
  meme $ID/REF_kcnj_canonical.faa -protein -oc $ST/meme -nmotifs 10 -minw 6 -maxw 50 -mod zoops
fi
python3 scripts/11_motif_matrix.py


step "9. kcnj15 locus"
python3 scripts/12_kcnj15_locus.py utr
python3 scripts/12_kcnj15_locus.py lead
python3 scripts/12_kcnj15_locus.py cds


step "10. Predicted structures"
python3 scripts/13_structure_analysis.py
conda run --no-capture-output -n "$PYMOL_ENV" pymol -cq scripts/14_structure_figures.pml


step "11. Phylogeny and synteny figures"
python3 scripts/15_tree_figure.py
python3 scripts/16_synteny_figure.py

rm -rf data/tmp
step "Finished: tables are in results/, figures are in figures/"
