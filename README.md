# Population Genomics of *Ewingella americana*: Taxonomic Reassessment, Clonal Structure, Pan-Genome Diversity and Mobile Genetic Elements

## Overview

*Ewingella americana* (family *Yersiniaceae*, NCBI Taxon ID: 41202) is a bacterial pathogen with a dual ecological profile: it acts as an opportunistic pathogen in clinical settings and serves as the primary etiological agent of internal stipe browning rot in commercially cultivated mushrooms (*Flammulina filiformis* and *Agaricus bisporus*).

Despite its clinical and agricultural importance, genomic investigations have historically been hindered by uncurated sequence deposits, taxonomic misclassifications, and massive clonal redundancy in public repositories. This repository provides the complete analytical pipeline, curated genomic metadata, pan-genome presence/absence matrices, outgroup-verified phylogenomic trees, and resistome/mobilome profiles supporting the definitive population-scale study of the genus ***Ewingella*** ($N = 36$).

---

## Key Scientific Findings

1. **Taxonomic Reassessment and Species Boundaries:**
   - Population-scale FastANI profiling across 36 global assemblies confirmed 31 *E. americana* *sensu stricto* genomes ($\ge 96.39\%$ ANI to type strain `CCUG 14506T`).
   - Polar isolates `E4` (Arctic snow) and `L47` (Antarctic soil) share only ~82.9% ANI with *Ewingella* and harbor the `oqxB9` marker, excluding them from the genus.
   - *De novo* assembly of agricultural plant epiphyte isolate `NIAB_247711` (Cambridge, UK; NCBI SRA `SRR32280131`) revealed **99.09% FastANI** to type strain `20WA0182T` of the proposed taxon ***Ewingella allii***. Furthermore, *E. allii* and cheese isolate *Ewingella docleensis* (`CoE-038-23T`) share 97.74% ANI and 81.0% dDDH with each other, forming a distinct sister-species complex separated from *E. americana* (93.5%–94.0% ANI).

2. **Demasking Clonal Distortion in Pan-Genome Architecture:**
   - High-stringency Panaroo modeling across 31 verified genomes resolved 5,545 orthologous gene families with an open pan-genome ($\gamma = 0.0671 \pm 0.0008, \; R^2 = 0.9965$) and an asymptotic core of 2,658.7 gene families.
   - Provenance tracing revealed that 20 clinical deposits represent a single clonal expansion (`CIP 81.94` / `CIP 91.84` / `HMF3896_M4G`; mean 4.8 core SNPs; separated by 53,765 SNPs [2.06%] from type `CCUG 14506T`), alongside six municipal bioreactor replicates.
   - Dereplicating the cohort into **six authentic biological lineages** expanded the representative core genome from 2,900 (52.30%) to **3,566** gene families (64.31%) and unmasked 1,015 lineage-specific singletons (18.30%), demonstrating that public database redundancy artificially suppresses the core genome by 23%.

3. **Core-Genome Phylogenomics and Basal Purifying Selection:**
   - Maximum-likelihood phylogenomics inferred from 2,634,233 bp concatenated core alignments (2,900 loci; IQ-TREE 2 GTR+F+I+R4, 1,000 UFBoot), validated by a 9-taxon triple-outgroup tree rooted on *Yersinia enterocolitica* 8081, confirmed *E. americana* monophyly (100/100 support).
   - Mushroom pathogen `B6-1` forms the earliest-diverging basal lineage, evolving under strong purifying selection (100% of core protein-coding genes with Nei-Gojobori $K_a/K_s < 1.0$).

4. **Mobilome Architecture and Resolution of Integron Contradictions:**
   - Profiling with standalone IntegronFinder v2.0.6 confirmed that all 31 verified *E. americana* genomes possess **zero complete integrons, zero In0 elements, zero CALIN elements, and zero attC sites**, definitively refuting earlier uncurated web-server reports claiming species-wide class 1 integrons.
   - Clinical lineages maintain a plasmid-free resistome governed by chromosomal AmpC $eta$-lactamase (`blaACT`) and conserved multidrug efflux systems (AcrAB-TolC, EmrB, DinF).
   - Mushroom pathogen `B6-1` uniquely acquired megaplasmid `pB61a` (331.4 kb; `oqxAB`, tellurite/arsenic clusters) and conjugative plasmid `pB61b` (105.3 kb; predicted 32-kb Type-F T4SS), both encoding non-canonical RepA/RepB replicons.

5. **Universal Chitinolysis vs. Pathogen Adaptation Determinants:**
   - Complete chitinolytic cascades (GH18 chitinase ChiD, secondary GH18, chitobiose PTS `chbA`, chitin deacetylase `cdaD`, NagA-Z operon) are conserved across 100% of genomes, indicating ancestral saprophytic chitin-metabolic potential.
   - In contrast, candidate host-colonization determinants in strain `B6-1` include candidate lysozyme inhibitor `LprI`, fimbrial adhesin usher-chaperone clusters (`fimD`/`fimC`), and a Type VI secretion PAAR spike.

---

## Repository Structure

```
ewingella-americana-genomics/
├── README.md                                # Comprehensive documentation and study overview
├── LICENSE                                  # MIT open-source license
├── run_pipeline.sh                          # Master Bash reproduction pipeline
├── metadata_all_ewingella.csv               # Curated metadata for all 36 genomes
├── Table_S1_Genomic_Inventory.csv           # S1: Assembly metrics, accessions, FastANI & QC
├── Table_S2_PanGenome_Matrix.csv            # S2: 5,545 gene family presence/absence matrix (N=31)
├── Table_S3_Pan_Dynamics_Sampling.csv       # S3: Multi-strategy pan-genome sampling comparisons
├── Table_S4_Resistome_Virulome.csv          # S4: ABRicate & VFDB resistome and virulence profiling
├── Table_S5_IntegronFinder.csv              # S5: IntegronFinder v2.0.6 element audit
├── supplementary/
│   └── Supplementary_Tables_S1_S5.xlsx      # Consolidated multi-sheet Excel workbook
├── figures/
│   ├── Figure_1_Workflow.pdf / .png         # Figure 1: Bioinformatic curation workflow
│   ├── Figure_2_FastANI.pdf / .png          # Figure 2: Whole-genome FastANI species demarcation
│   ├── Figure_3_PanGenome.pdf / .png        # Figure 3: Pan-genome dynamics & clonal dereplication
│   ├── Figure_4_Phylogenomics.pdf / .png    # Figure 4: Maximum-likelihood strict core phylogenomics
│   └── Figure_5_Resistome.pdf / .png        # Figure 5: Resistome, mobilome & chitinolysis
├── phylogeny/
│   ├── strict_31_Q1.treefile                # IQ-TREE 2 ML tree (2,634,233 bp nucleotide core)
│   ├── strict_31_Q1.contree                 # Consensus tree with 1,000 UFBoot support
│   ├── core_aa_validation.treefile          # Amino-acid validation tree (120,790 aa)
│   ├── pairwise_core_snps_31.csv            # Pairwise core-genome SNP distance matrix
│   ├── pairwise_pdist_matrix_31.csv         # Pairwise p-distance matrix
│   ├── Table_Core_KaKs_NG86.csv             # Locus-by-locus Nei-Gojobori Ka/Ks selection table
│   └── Table_Phylogenetic_Models_and_Support.csv # IQ-TREE ModelFinder benchmarks
└── scripts/
    ├── 01_prepare_genomes_and_metadata.py   # Genome retrieval and metadata harmonization
    ├── 02_run_fastani.sh                    # FastANI species boundary execution
    ├── 03_run_abricate.sh                   # ABRicate AMR/virulome screening
    ├── 04_analyze_ani_and_qc.py             # ANI distance clustering and QC filtering
    ├── 05_analyze_abricate.py               # Resistome and mobilome reconciliation
    ├── 06_prepare_panaroo_gffs.py           # Annotation preparation for Panaroo
    ├── 07_run_panaroo.sh                    # Panaroo pan-genome pipeline
    ├── 08_run_iqtree.sh                     # IQ-TREE 2 maximum-likelihood phylogenomics
    ├── 09_analyze_pangenome.py              # Heaps' law & rarefaction modeling
    ├── 11_plot_phylogeny_and_summary.py     # Phylogenetic visualization
    └── 04_R_pangenome_and_comparative_genomics.R # R analysis and figure plotting
```

---

## Computational Reproduction

### Prerequisites
- Linux / WSL2 (Ubuntu 20.04+)
- Conda / Mamba environment manager
- Core tools: `fastANI` (v1.33+), `panaroo` (v1.3.4+), `iqtree2` (v2.0.7+), `abricate` (v1.0.1+), `integron_finder` (v2.0.6+), `spades` (v3.15.5+), `prokka` (v1.14.6+)
- Python libraries: `pandas`, `scipy`, `numpy`, `matplotlib`, `seaborn`, `biopython`, `openpyxl`

### Execution
```bash
# Clone the repository
git clone https://github.com/asweshvaran/ewingella-americana-genomics.git
cd ewingella-americana-genomics

# Make execution script runnable
chmod +x run_pipeline.sh

# Run end-to-end pipeline
./run_pipeline.sh
```

---

## Citation

If using the curated assemblies, pan-genome matrices, or methodological workflows from this repository, please cite:

> **Asweshvaran, R. (2026).** Population Genomics of *Ewingella americana*: Taxonomic Reassessment, Clonal Structure, Pan-Genome Diversity and Mobile Genetic Elements. *Antonie van Leeuwenhoek* (Springer Nature), submitted.

---

## Contact

**Asweshvaran R.**  
Department of Microbiology  
Shriram Research Institute  
14–15, Sadarmangala Industrial Area, Whitefield Road,  
Bengaluru 560048, Karnataka, India  
