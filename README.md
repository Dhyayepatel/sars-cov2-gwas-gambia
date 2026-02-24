# GWAS of SARS-CoV-2 Variants in a Gambian Cohort

> MSc Health Data Science Project — University of [Your University], 2024  
> Investigating whether viral genetic variation (SNPs) in SARS-CoV-2 influences patient viral load (CT values) in The Gambia.

---

## Overview

This project applies genome-wide association study (GWAS) methodology to **1,643 SARS-CoV-2 whole-genome sequences** from The Gambia, each matched with RT-PCR cycle threshold (Ct) values as a proxy for viral load. The goal was to identify viral SNPs significantly associated with differences in Ct across infected individuals — and to rigorously assess whether any such associations reflect genuine biology or population structure artefacts.

**Key finding:** A single SNP at position 15,222 initially appeared significant, but lost significance entirely when the genetically divergent A.29 lineage was removed — a textbook demonstration of how population stratification can produce spurious associations in viral GWAS.

---

## Repository Structure

```
├── 01_installation/          # Tool setup: PLINK2, Pangolin, IQ-TREE, Jvarkit, JBrowse
├── 02_pangolin_analysis/     # Lineage assignment and top 10 lineage extraction
├── 03_alignment/             # MAFFT multiple sequence alignment scripts
├── 04_fasta_to_plink/        # VCF generation (Jvarkit), biallelic filtering, PLINK conversion
├── 05_ct_data_cleaning/      # CT count cleaning and PLINK phenotype formatting (R)
├── 06_gwas/                  # QC, LD pruning, simple linear regression (PLINK2)
├── 07_phylogenetic_tree/     # Maximum likelihood tree construction (R: ape, phangorn, ggtree)
├── 08_pca/                   # Population structure analysis (PLINK2 + R)
├── 09_three_model_gwas/      # Three-model framework for population structure control
├── 10_snp_annotation/        # SNP-to-gene mapping and amino acid analysis
├── automate_lineages.sh      # Automated pipeline for all top 10 lineages
└── README.md
```

---

## Methods Summary

### Data
- **1,643 SARS-CoV-2 whole-genome sequences** from routine surveillance in The Gambia
- **Matched Ct values** (RT-PCR cycle threshold) used as a continuous viral load phenotype
- Lineages assigned using **Pangolin**; top 10 most prevalent lineages extracted for lineage-specific analyses

### Pipeline

| Step | Tool | Description |
|------|------|-------------|
| Lineage assignment | Pangolin | Batch processing of 1,643 genomes |
| Multiple sequence alignment | MAFFT (--auto) | Per-lineage and all-lineage alignments |
| Variant calling | Jvarkit (msa2vcf) | FASTA → VCF conversion |
| Biallelic filtering | bcftools | Remove multiallelic sites |
| GWAS format conversion | PLINK2 | VCF → .bed/.bim/.fam |
| Quality control | PLINK2 | SNP missingness <2%, sample missingness <2%, MAF >0.01 |
| LD pruning | PLINK2 | Window=50, step=5, r²=0.2 |
| Phylogenetic tree | ape, phangorn, ggtree (R) | NJ → ML optimisation → bootstrap (n=30) |
| Population structure | PLINK2 + R | PCA, outlier detection, pairwise centroid distances |
| GWAS | PLINK2 | Linear regression (Ct ~ SNP) |
| Three-model framework | PLINK2 | Primary / Sensitivity / Stratified models |
| SNP annotation | bedtools + R | Map SNPs to CDS regions, amino acid analysis |

### Three-Model GWAS Framework

To address population stratification — particularly the genetically divergent A.29 lineage identified in PCA:

- **Model 1 (Primary):** All samples; A.29 binary covariate + PC1/PC2
- **Model 2 (Sensitivity):** Outlier samples removed; same covariates
- **Model 3 (Stratified):** A.29 samples fully excluded; PC1/PC2 only

---

## Key Results

### Lineage Distribution
The 10 most prevalent lineages accounted for the majority of genomes, dominated by Delta sub-variants:

| Lineage | Genomes | Wave |
|---------|---------|------|
| AY.34.1 | 262 | 3rd wave (Delta) |
| B.1.416 | 184 | 1st wave |
| B.1.617.2 | 160 | 3rd wave (Delta) |
| B.1.1.7 | 140 | 2nd wave (Alpha) |
| B.1.525 | 129 | 2nd wave (Eta) |
| B.1 | 133 | Early |
| BA.1.1 | 28 | 4th wave (Omicron) |

### Population Structure
PCA on 1,187 genomes across 36 lineages revealed **A.29 as a clear genetic outlier**, sitting in the lower-left quadrant of PC space with centroid distances of 0.335–0.349 from all other major lineages. All 9 A.29 samples were classified as outliers (distance threshold: mean + 2SD = 0.091).

### GWAS Results

| Analysis | SNPs Tested | Significant SNPs | Notes |
|----------|-------------|-----------------|-------|
| All lineages (simple regression) | 70 | 1 (pos. 15,222) | In pp1ab polyprotein |
| All lineages (Model 1) | 70 | 1 (P=2.71×10⁻⁵) | A.29 covariate included |
| All lineages (Model 3) | 70 | 0 (P=5.62×10⁻³) | A.29 removed — signal lost |
| Top 10 lineages (lineage-specific) | 9–20 per lineage | 0 | No genome-wide hits |
| B.1.416 lineage-specific | 11 | 1 (pos. 1,066) | Synonymous; T→C in ORF1a |

**The SNP at position 15,222 is a population stratification artefact.** Its significance disappeared when A.29 samples were excluded, confirming it reflects lineage frequency differences rather than a genuine SNP-Ct association.

**The B.1.416 SNP (pos. 1,066)** is a synonymous T→C mutation in ORF1ab (codon 267; asparagine unchanged). It is a lineage-specific finding with a large effect size (β=7.32, SE=2.20) but represents a within-lineage signal rather than a population-wide discovery.

### Conclusion
No viral SNPs were found to genuinely influence CT-measured viral load after rigorous population structure control. This suggests that **host factors and technical variables dominate Ct variation**, and that viral genetic effects — if they exist — are small relative to the sample sizes available in this cohort.

---

## Environment & Dependencies

### System
- Windows 10/11 with WSL2 (Ubuntu 24.04)
- R ≥ 4.3

### Bioinformatics Tools
```bash
# Linux tools (install via apt or conda)
PLINK2          # Genetic association analysis
Pangolin        # SARS-CoV-2 lineage assignment (conda: bioconda)
MAFFT           # Multiple sequence alignment
bcftools        # VCF filtering
bedtools        # Genomic interval operations
seqkit          # FASTA manipulation
Jvarkit         # FASTA → VCF (msa2vcf; v20200206)
IQ-TREE         # Phylogenetic inference (optional)
```

### R Packages
```r
install.packages(c("ggplot2", "dplyr", "readr", "qqman", 
                   "ape", "phangorn", "seqinr", "magick",
                   "RColorBrewer", "viridis", "ggrepel",
                   "cluster", "factoextra", "gridExtra",
                   "moments", "nortest", "cowplot"))

# Bioconductor
BiocManager::install(c("Biostrings", "ggtree"))
```

---

## Reproducing the Analysis

> **Note:** All file paths in scripts use the original project directory. Update the `BASE_DIR` variable in `automate_lineages.sh` and path variables in R scripts to match your local setup before running.

```bash
# 1. Clone the repository
git clone https://github.com/yourusername/sars-cov2-gwas-gambia.git
cd sars-cov2-gwas-gambia

# 2. Activate Pangolin environment and assign lineages
conda activate pangolin
bash 02_pangolin_analysis/run_pangolin.sh

# 3. Run automated lineage pipeline (alignment → VCF → PLINK → GWAS → plots)
bash automate_lineages.sh

# 4. Run R scripts for CT cleaning, PCA, three-model GWAS, and plots
Rscript 05_ct_data_cleaning/clean_ct.R
Rscript 08_pca/pca_analysis.R
Rscript 09_three_model_gwas/three_model_plots.R
```

---

## Figures

| Figure | Description |
|--------|-------------|
| Phylogenetic tree | ML tree of all Gambian genomes, coloured by Pangolin lineage |
| Lineage bar plot | Frequency of each lineage across ~1,600 genomes |
| CT distribution | Histogram + QQ plot of cleaned Ct values (n=1,037) |
| PCA — all lineages | PC1 vs PC2 coloured by lineage; A.29 clearly separated |
| PCA — top 10 lineages | Per-lineage PCA panels showing within-lineage structure |
| Manhattan plots | SNP associations for all lineages and each top-10 lineage |
| QQ plots | Genomic inflation assessment (λ values) per analysis |
| Three-model comparison | Combined Manhattan and QQ plots across Models 1–3 |

---

## Limitations
- Low genomic sampling density (~1,643 of ~11,900 reported cases sequenced)
- Urban sampling bias; rural populations underrepresented
- Severity bias: testing prioritised symptomatic cases
- Small per-lineage sample sizes limit statistical power for lineage-specific GWAS
- Ct values are a proxy for viral load, not a direct measure of infectious virus

---

## Citation / Academic Context

This project was submitted in partial fulfilment of the MSc Health Data Science degree (2023–2024). The genomic data derives from surveillance described in:

> Kanteh A, et al. Genomic surveillance of SARS-CoV-2 in The Gambia reveals limited data and global disparities. *Viruses.* 2023;15(3):706.

---

## Acknowledgements

Supervised by Stephan Hue, whose guidance throughout the GWAS workflow and feedback on drafts was invaluable.

AI usage: Claude.ai was used for language refinement. All analysis, code, data cleaning, and conclusions are the author's own work.

---

*MSc Health Data Science | September 2024*

