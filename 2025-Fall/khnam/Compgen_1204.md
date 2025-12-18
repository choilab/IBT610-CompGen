# 1. Whole Genome MLST Analysis (wgMLST)
## Goal
Instead of the conventional 7-housekeeping-gene MLST scheme, a whole-genome–based MLST approach was applied.

##  Tools Used
* **chewBBACA (v3.x):** For allele calling and schema creation.
* **Prodigal:** For training file generation.
* **GrapeTree / PHYLOViZ:** For Minimum Spanning Tree (MST) visualization.

## Result
According to the wgMLST analysis, the ST-2 group—previously classified by conventional MLST—formed a star-burst topology, indicating that this clone underwent rapid clonal expansion within the hospital environment after becoming established.

<img width="593" height="437" alt="image" src="https://github.com/user-attachments/assets/e63a1e7d-5f5c-45c9-81f8-49c50f7661cb" />

The red circles indicate that the ST-2 group originated from a single ancestor and rapidly expanded within the hospital environment.

---

# 2. Virulence & Risk Assessment Index
## Goal
The potential hazard of each strain was quantified by calculating a Combined Risk Index based on the accumulation of virulence factors and antibiotic resistance genes.

## Tools used
* **AMRFinderPlus:** Gene detection
* **Python (Pandas):** Custom scoring script

## Result
**MRSE (ST-2):** Showed significantly higher Risk Indices (Avg > 10) compared to MSSE strains

---

# 3. BGC Profiling (antiSMASH)
## Goal
To identify Biosynthetic Gene Clusters (BGCs), specifically focusing on Siderophores (survival) and Bacteriocins (competition).

## Tools used
* **antiSMASH 7.0:** BGC detection

## Result
The siderophore system, which is essential for survival, was present in all strains. Bacteriocins (Ripp-like) were evenly distributed across both MRSE and MSSE groups, with 16 out of 124 isolates carrying these genes.

---

# Integration
## Goal
To merge all datasets (AMR, Virulence, BGC) into a single master table for statistical analysis and visualization (iTOL)

<img width="744" height="697" alt="image" src="https://github.com/user-attachments/assets/efadef37-1b86-47e8-bb9c-66d5e6e688fb" />

Strains with a red background represent MRSE, while those with a green background represent MSSE. Strains harboring bacteriocins are indicated by blue rings. The virulence index is displayed in the form of bar graphs.
Although no strains in the ST-2 group carried bacteriocin genes, they exhibited an average risk score above 13, indicating that they are potentially hazardous.


