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
* Principle: **Compute Average Nucleotide Identity (ANI)** efficiently using a k-mer–based algorithm instead of whole-genome alignment.
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
* Method:
  1. Map genome annotations of each strain to KEGG pathways.
  2. Screen for the absence of essential genes for amino acid biosynthesis (e.g., *ilvH* for valine, *hisE* for histidine).
  3. Strains lacking key pathway genes are classified as having auxotrophy for the corresponding amino acids.

#### 3-2. AMR & Virulence Screening
* Tool: **AMRFinderPlus**
* Principle: Use NCBI’s curated databases and Hidden Markov Models (HMMs) to identify antibiotic resistance genes (AMR) and virulence factors with higher accuracy than simple homology searches

#### 3-3. Risk Assessment (Composite Risk Index)
* Tool: **Quantitative Scoring System**
* Principle: Go beyond simple gene presence/absence by quantifying the potential risk of each strain, and perform statistical comparisons between groups (e.g., MRSE vs MSSE).
* Calculation Logic:  $$\text{Risk Index} = \sum (\text{AMR Genes}) + \sum (\text{Virulence Factors})$$
* Statistical Analysis: Statistically test differences in index values between groups and display results with boxplots

#### 3-4. Secondary Metabolite Analysis
* Tool: **antiSMASH**
* Principle: Analyze genome-wide BGC patterns to predict the potential production of secondary metabolites, such as siderophores (iron acquisition) and bacteriocins (competitor inhibition).

### 4. High-Resolution Analysis
Perform high-resolution analyses to track microevolution and transmission routes within populations.

#### 4-1. wgMLST (Whole Genome MLST)
* Tool: **chewBBACA**
* Principle: Perform allele calling across the entire genome (thousands of loci) instead of just seven genes, enabling detection of fine-scale differences between strains.

#### 4-2. Visualization
* Tool: **GrapeTree**, **ITOL**
* Method: Construct a Minimum Spanning Tree (MST) based on wgMLST results, and integrate metadata (e.g., ST, AMR profile, Risk Index) for visualization

---

# Result
## 1. Construction of a Clean Dataset
| All Genome | Annotated |               
|:----------:|:---------:|              
|  5316      |   5038    |

| Complete | Chromosome | Scaffold | Contig |
|:--------:|:----------:|:--------:|:------:|
|   5038   |   4766     |    4719  |  3316  |

All the contents in column are cited from the NCBI database.

In annotated genomes, there are categorized by assembly level.

In complete assembly level, 272 complete Refseq genomes are downloaded at first.

### 1-1. QC & Dereplication
Exclude low-quality genomes with contamination > **5%** or completeness < **90%** (0 genomes removed).

Cluster genetically identical clones at **99%** ANI, removing 148 redundant genomes out of 272.

**Final clean dataset:** Finalized a set of 124 representative genomes.


## 2. Genotyping & Population Structure
### 2-1. MLST Distribution

| FILE | SCHEME | ST | arcC | aroE | gtr | mutS | ... |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| GCF\_001.fna | sepidermidis | 184 | 16 | 1 | 2 | 1 | ... |
| GCF\_002.fna | sepidermidis | 1016 | 1 | 13 | 7 | 6 | ... |
| GCF\_003.fna | sepidermidis | 20 | 1 | 1 | 2 | 2 | ... |
| GCF\_004.fna | sepidermidis | 20 | 1 | 1 | 2 | 2 | ... |
| GCF\_005.fna | sepidermidis | - | 1 | 1 | 2 | - | ... |

| ST | Count of ST |
| :--- | :--- |
| 2 | 18 |
| 5 | 5 |
| 20 | 3 |
| 6 | 3 |
| - | 13 |
| ... | ... |
| **Total** | **124** |

* A total of 69 sequence types (STs) were identified among 124 strains.
* **ST-2** accounted for ~14.5% of strains and was identified as the predominant pathogenic clone.

### 2-2. Pan-genome Analysis

| Category | Gene Count |
| :--- | :--- |
| **Core genes** (99% <= strains <= 100%) | **1522** |
| Soft core genes (95% <= strains < 99%)| 189 |
| Shell genes (15% <= strains < 95%) | 948 |
| Cloud genes (0% <= strains < 15%) | 7355 |
| **Total genes** (0% <= strains <= 100%) | **10014** |

Approximately 85% of all genes are classified as accessory genes

### 2-3. GWAS
* **ST-2 (Pathogenic Type):**
    * A highly significant marker gene, **group_4200 (P < 1.0×10⁻¹³, sensitivity 100%)**, was identified.
    * BLAST analysis: The gene encodes a protein derived from a **bacteriophage (*Caudoviricetes sp.*)**.
  

## 3. Functional Profiling
### Overview

| Functional Category | Gene Count | Proportion(%) | Biological Implication  |
| :--- | :---: | :---: | :--- |
| **Hypothetical / Unknown** | **6,543** | 65.34% | Uncharacterized genes with unknown functions |
| **Enzymes (Metabolism)** | **1,549** | 15.47% | Basic metabolic enzymes essential for survival and growth (e.g., kinases, synthases) |
| **Mobile Elements** | **348** | 3.47% | HGT-related genes driving horizontal gene transfer, including phages and transposons |
| **Transporter / Membrane** | **291** | 2.91% | Genes involved in nutrient acquisition, waste removal, and drug efflux pumps |
| **Translation / Ribosome** | **202** | 2.02% | Ribosomal and tRNA-related genes, mostly part of the core genome |
| **Regulation** | **197** | 1.97% | Factors controlling gene expression in response to environmental changes (e.g., repressors, regulators) |
| **Virulence / Resistance** | **90** | 0.90% | Genes directly related to host attack (toxins) or antibiotic resistance |
| **Other Function** | **794** | 7.92% | Other structural and accessory proteins |
| **Total** | **10,014** | 100.0% | - |

Based on annotations generated by Prokka, functional analysis at the gene level was performed across seven categories.

