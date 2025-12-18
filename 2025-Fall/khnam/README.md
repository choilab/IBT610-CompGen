# 🧬 2025-Fall Compuational Genomics Project
**Project Period:** 10/16/2025 ~ 12/18/2025

**Model Organism:** *Staphylococcus epidermidis*

**Researcher:** KIHYUN NAM

---

# Introduction
## Main Goal
This project uses *Staphylococcus epidermidis* a common skin commensal that can act as an opportunistic pathogen as a model organism to build a standardized analysis pipeline.

The pipeline covers the entire workflow, from constructing a clean, high-quality genomic dataset to functional profiling and data visualization.

## Research Detailed
### 1. Construction of a Clean Dataset
*  **QC(Quality Control) & Dereplication:** Apply **CheckM** and **FastANI** to filter out low-quality and redundant genomes.
*  **Bias Reduction:** Build a representative genome set to minimize bias and improve statistical reliability.

#### 1-1. QC (Quality Control)
* Tool: **CheckM**
* Principle: Assess genome completeness and contamination based on single-copy marker genes.
* Criteria: Completeness ≥ **90%**, Contamination ≤ **5%**

#### 1-2. Dereplication
* Tool: **FASTANI**
* Principle: Compute average nucleotide identity (ANI) efficiently using a k-mer–based algorithm instead of whole-genome alignment.
* Method: Genomes with ≥ **99%** ANI are clustered as the same clone, and the highest-quality genome from each cluster is selected as a representative genome.
* Note: Since genomes with < **1%** divergence are treated as the same clone during dereplication, some biologically distinct clones may be filtered out and excluded from the final dataset.


### 2. Genotyping & Population Structure

Assess genetic diversity and population structure among the representative genomes.

#### 2-1. MLST (Multi-Locus Sequence Typing)
* Tool: **mlst (Torsten Seemann)**
* Principle: Assign unique **sequence types (STs)** by comparing sequences of seven housekeeping genes against the PubMLST database

#### 2-2. Pan-genome Analysis
* Tool: **Roary**
* Principle: Classify genes of all strains into core genome (present in ≥ **99%** of strains) and accessory genome, and generate a gene presence/absence matrix.

#### 2-3. GWAS (Genome-Wide Association Study)
* Tool: **Scoary**
* Principle: Identify significant marker genes by statistically testing the association between gene presence/absence (Roary results) and phenotypic traits (e.g., ST-2) using Fisher's exact test


### 3. Functional Profiling
Evaluate metabolic requirements (Auxotroph), AMR, and pathogenic potential from genomic information

#### 3-1. Metabolic Capability (Auxotrophy Analysis)
* Tool: **GapMind** and **KEGG** pathway database
* Principle: Determine the presence of key enzyme genes in biosynthesis pathways to classify strains as prototrophs or auxotrophs
