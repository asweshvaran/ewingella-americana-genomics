#!/usr/bin/env bash
# ==============================================================================
# Pipeline: Comparative Population Genomics and Pan-Genome Analysis of
#           True Ewingella americana
# Author: Asweshvaran R. (Department of Microbiology, Shriram Institute for Industrial Research)
# Journal: Journal of Applied Microbiology (Oxford University Press)
# ==============================================================================

set -euo pipefail

echo "=== Stage 1: Short-Read SRA Assembly of Novel Plant Epiphyte NIAB_247711 ==="
# Download SRR32280131 via fasterq-dump
# fasterq-dump --split-files -O raw_reads/ SRR32280131
# Quality trimming with fastp
# fastp -i raw_reads/SRR32280131_1.fastq -I raw_reads/SRR32280131_2.fastq \
#       -o trimmed_1.fq.gz -O trimmed_2.fq.gz --detect_adapter_for_pe --html fastp_report.html
# De novo assembly using SPAdes v3.15.5 isolate mode
# spades.py --isolate -1 trimmed_1.fq.gz -2 trimmed_2.fq.gz -o spades_niab_247711 -t 8 -m 32
# Filter contigs >= 500 bp
# seqkit seq -m 500 spades_niab_247711/scaffolds.fasta > genomes_fna/NIAB_247711.fna
# Annotate using Prokka v1.14.6
# prokka --outdir prokka_niab_247711 --prefix NIAB_247711 --genus Ewingella --species americana --kingdom Bacteria genomes_fna/NIAB_247711.fna

echo "=== Stage 2: Population-Scale FastANI Taxonomic Demarcation ==="
# Compute all-vs-all FastANI across 36 genomes
# fastANI --ql fastani_genome_list.txt --rl fastani_genome_list.txt -o fastani_all_vs_all.txt -t 8 --matrix
# Compute FastANI of novel isolate NIAB_247711 against type strain CCUG 14506T and B6-1
# fastANI -q genomes_fna/NIAB_247711.fna -r genomes_fna/GCF_008693045.1.fna -o fastani_niab_vs_type.txt

echo "=== Stage 3: High-Stringency Pan-Genome Reconstruction (Panaroo) ==="
# Run Panaroo strict mode across 31 verified true E. americana sensu stricto genomes
panaroo -i gff_panaroo/*.gff \
        -o panaroo_strict_31 \
        --clean-mode strict \
        -a core \
        --aligner mafft \
        --core_threshold 0.98 \
        -t 8 \
        --remove-invalid-genes

echo "=== Stage 4: Core Genome Maximum-Likelihood Phylogenomics (IQ-TREE 2) ==="
# Build phylogenomic tree on 2.63 Mb concatenated core alignment under GTR+F+I+R4 with 1000 UFBoot
iqtree2 -s panaroo_strict_31/core_gene_alignment.aln \
        -m GTR+F+I+R4 \
        -B 1000 \
        -T AUTO \
        --prefix iqtree_strict_31/strict_31

echo "=== Stage 5: Resistome, Mobilome, and Virulome Screening (ABRicate) ==="
# Screen all 36 genomes across CARD, NCBI AMRFinderPlus, ResFinder, PlasmidFinder, and VFDB
for db in ncbi card resfinder plasmidfinder vfdb; do
    echo "Running ABRicate on database: $db"
    abricate --db "$db" genomes_fna/*.fna > "abricate_results/results_${db}.tsv"
    abricate --summary "abricate_results/results_${db}.tsv" > "abricate_results/summary_${db}.tsv"
done

echo "=== Stage 6: Mathematical Pan-Genome Modeling (Heaps' Law & Exponential Decay) ==="
# Fit Heaps' law (n = kappa * N^gamma) and core decay (n = (kappa_c - Omega)*exp(-N/tau) + Omega)
# python3 scripts/fit_heaps_law.py

echo "=== Workflow Execution Complete ==="
