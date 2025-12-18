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

#### 3-3. Secondary Metabolite Analysis
* Tool: **antiSMASH**
* Principle: Analyze genome-wide BGC patterns to predict the potential production of secondary metabolites, such as siderophores (iron acquisition) and bacteriocins (competitor inhibition).

#### 3-4. Risk Assessment (Composite Risk Index)
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

### 3-1. Metabolic Auxotrophy
* Summary of Amino Acid Biosynthesis Profile

| Category | Amino Acids | Phenotype | Key Mechanism / Note |
| :--- | :--- | :--- | :--- |
| **Conserved** | **Glu**, **Gln**, **Ala**, **Gly**, **Pro**, **Lys**, **Met**, **Thr**, **Cys** | ✅ **Prototroph**| Pathways linked to Central Metabolism are 100% conserved. |
| **Variable** | **Arg**, **Leu**, **Trp** |  ⚠️ **Variable**  | Strain-dependent presence/absence indicating ongoing evolutionary divergence. |
| **Auxotroph** | **Val**, **Ile**, **His**, **Phe**, **Tyr**, **Ser**, **Asn** |  🔴 **Auxotroph**  | Essential genes (*ilvH*, *hisE*, *tyrB*) are deleted across the entire population. |

| Amino Acid Pathway | Determinant Marker (Gene) | Prevalence (n=124) | Phenotype Prediction | Note (Mechanism) |
| :--- | :--- | :---: | :--- | :--- |
| **Glutamate** | Glutamate synthase (`gltB`) | 124 (100%) | ✅ **Prototroph** | Conserved Core |
| **Glutamine** | Glutamine synthetase (`glnA`) | 124 (100%) | ✅ **Prototroph** | Conserved Core |
| **Alanine** | Alanine dehydrogenase (`ald`) | 124 (100%) | ✅ **Prototroph** | Conserved Core |
| **Glycine** | Glycine hydroxymethyltransferase (`glyA`) | 123 (99%) | ✅ **Prototroph** | Conserved Core |
| **Proline** | Pyrroline-5-carboxylate reductase (`proC`) | 124 (100%) | ✅ **Prototroph** | Conserved Core |
| **Lysine** | Diaminopimelate decarboxylase (`lysA`) | 124 (100%) | ✅ **Prototroph** | Conserved Core |
| **Methionine** | Methionine synthase (`metE`) | 124 (100%) | ✅ **Prototroph** | Conserved Core |
| **Threonine** | Threonine synthase (`thrC`) | 124 (100%) | ✅ **Prototroph** | Conserved Core |
| **Cysteine** | Cysteine synthase (`cysK`) | 123 (99%) | ✅ **Prototroph** | Conserved Core |
| **Tryptophan** | Tryptophan synthase alpha (`trpA`) | 114 (92%) | ⚠️ **Variable** | Deletion confirmed in a small subset of strains |
| **Leucine** | 2-isopropylmalate synthase (`leuA`) | 94 (76%) | ⚠️ **Variable** | Distinct inter-strain variation |
| **Arginine** | Argininosuccinate synthase (`argG`) | 89 (72%) | ⚠️ **Variable** | Pathway decay driven by `argG` deletion |
| **Valine / Isoleucine** | Acetolactate synthase small (`ilvH`) | **0 (0%)** | 🔴 **Auxotroph** | **Large subunit only; functional loss** |
| **Histidine** | P-ribosyl-ATP pyrophosphatase (`hisE`) | **0 (0%)** | 🔴 **Auxotroph** | Complete deletion of intermediate pathway steps |
| **Phenylalanine** | Aromatic-AA aminotransferase (`tyrB`) | **0 (0%)** | 🔴 **Auxotroph** | Deletion of the terminal enzyme (`PPYAT`) |
| **Tyrosine** | Aromatic-AA aminotransferase (`tyrB`) | **0 (0%)** | 🔴 **Auxotroph** | Deletion of enzyme shared with Phe pathway |
| **Serine** | Phosphoserine aminotransferase (`serC`) | **0 (0%)** | 🔴 **Auxotroph** | Complete absence |
| **Asparagine** | Asparagine synthetase (`asnA`) | **0 (0%)** | 🔴 **Auxotroph** | Complete absence |

> **Legend:**
> * ✅ **Prototroph:** Essential gene retained in >95% of strains (Conserved).
> * ⚠️ **Variable:** Gene presence varies between strains (Strain-specific).
> * 🔴 **Auxotroph:** Determinant essential gene deleted in all strains (Absent).



**Visualization**
<img width="930" height="1125" alt="image" src="https://github.com/user-attachments/assets/9385e90d-c565-469a-9a11-f9c0854fe092" />



All 124 analyzed strains showed missing genes for seven amino acids, including **valine (Val)**, **isoleucine(Ile)**, and **histidine(His)**

For **leucine(Leu)**, **arginine(Arg)**, and **tryptophan(Trp)**, some strains are prototrophic while others are auxotrophic

**Note**
Only valine and arginine auxotrophy was confirmed experimentally, indicating a need for validation of computational predictions

### 3-2. AMR & Virulence Screening

| Rank | Gene Symbol | Function / Class | Prevalence (n=124) | Note |
| :--- | :--- | :--- | :---: | :--- |
| **1** | `hld` | Delta-hemolysin (Toxin) | **124 (100.0%)** | Species-specific core virulence factor |
| **2** | `fosB` | Fosfomycin resistance | **122 (98.4%)** | Intrinsic resistance trait |
| **3** | `blaZ` | Penicillin hydrolysis (Beta-lactamase) | **97 (78.2%)** | High prevalence of penicillin resistance |
| **4** | `arsB` | Arsenic efflux pump | **79 (63.7%)** | Heavy metal/Environmental stress resistance |
| **5** | `mecA` | **Methicillin resistance (PBP2a)** | **62 (50.0%)** | Determinant of **MRSE** (Pathogenic marker) |
| **6** | `aac(6')-aph(2'')` | Aminoglycoside resistance | **48 (38.7%)** | Resistance to Gentamicin, etc. |
| **7** | `icaC` | Biofilm formation (ICA operon) | **49 (39.5%)** | Critical for device-related infections |
| **8** | `qacR` | Antiseptic resistance repressor | **43 (34.7%)** | Regulation of antiseptic (QACs) resistance |


* **Key Finding:**
   
   * **`hld` (Toxin):** Present in 100% of strains
  
   * **`mecA` (Methicillin):** Present in 50% of strains, matching precisely with the ST-2 group
  
   * **`fosB`, `blaZ`:** Detected at high frequencies

A total of **47 unique resistance and virulence-associated genes** were identified

**Visualization**
<img width="5197" height="4749" alt="image" src="https://github.com/user-attachments/assets/c4ac9b24-a733-4267-95b9-79205ab21ed5" />

### 3-3. Secondary Metabolite Analysis

**Siderophores & Bacteriocins:** The **siderophore system**, essential for survival, was found in all strains. **Bacteriocins (Ripp-like)** were evenly distributed across MRSE and MSSE, with 16 of 124 isolates carrying these genes


### 3-4. Risk Assessment (Composite Risk Index)

<img width="2400" height="1800" alt="image" src="https://github.com/user-attachments/assets/0dda5e3e-72bb-4f56-8985-eea496826333" />

The **MRSE (resistant)** group showed a significantly higher Risk Index compared to the **MSSE (susceptible)** group (P < 0.001)

## 4. High-Resolution Analysis
### 4-1. wgMLST (Whole Genome MLST)

<img width="593" height="437" alt="image" src="https://github.com/user-attachments/assets/f2ffb5f7-3695-4ac3-a40c-40938c5811a0" />

According to the wgMLST analysis, the ST-2 group—previously classified by conventional MLST—formed a star-burst topology, indicating that this clone underwent rapid clonal expansion within the hospital environment after becoming established.

The red circles indicate that the ST-2 group originated from a single ancestor and rapidly expanded within the hospital environment.

### 4-2 Merged Visualization with Profiling Result

<img width="619" height="582" alt="image" src="https://github.com/user-attachments/assets/d7e048c0-b42f-473d-a90b-839163ae2ade" />

Strains with a red background represent MRSE, while those with a green background represent MSSE. Strains harboring bacteriocins are indicated by blue rings. The virulence index is displayed in the form of bar graphs.
Although no strains in the ST-2 group carried bacteriocin genes, they exhibited an average risk score above 13, indicating that they are potentially hazardous.
