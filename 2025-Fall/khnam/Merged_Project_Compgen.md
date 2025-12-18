# 🧬 Staphylococcus epidermidis Genomics Project Log
**Research Period:** Oct 30, 2025 – Ongoing  
**Target:** *Staphylococcus epidermidis* Representative Genomes (Selection from 272 RefSeq)  

---

## 📅 Data Collection
**Date:** 2025. 10. 30

### 🎯 Goal
A statistically robust set of representative strains was established by eliminating redundancy within the extensive NCBI dataset

### 🛠️ Task
1.  **DATA:** Downloaded 272 *S. epidermidis* complete genomes from the NCBI RefSeq database
2.  **Quality Control (QC):** Quality control of the 272 downloaded genomes was performed using the **CheckM**
3.  **Dereplication**
    * **Tool:** **FastANI**
      
    * **Criteria:** Strains exhibiting ≥99% **average nucleotide identity (ANI)** were clustered
      
    * **Result:** A total of **124 representative genomes** were selected from the initial set of 272 assemblies
      
    * **Limitation** Dereplication was performed using a 99% ANI similarity threshold, thus some genomes with less than 1% divergence may have been classified as redundant and consequently removed from the dataset
4.  **Verification:** An ANI heatmap was generated for the 124 genomes to confirm the diversity represented within the dataset
<img width="1052" height="947" alt="image" src="https://github.com/user-attachments/assets/85e111ea-6171-41b9-aa9c-bbf502efcda0" />

---

## 📅 Genotyping
**Date:** 2025. 11. 06

### 🎯 Goal
Determine the sequence types (STs) of the 124 selected isolates and validate their protein-level similarity

### 🛠️ Task
1.  **MLST Analysis**
    * **Tool:** **mlst** (7 housekeeping genes)
    * **Result:** A total of 69 sequence types (STs) were identified, with **ST-2 (18 strains, 14.5%)** emerging as the most predominant pathogenic clone
2.  **AAI Analysis (Proteome Similarity)**
    * **Tool:** **CompareM**
    * **Result:** The average amino acid identity (AAI) among the 124 strains ranged from **90% to 100%**, indicating a high level of genomic similarity
    * **Visualization**
    <img width="1144" height="1125" alt="image" src="https://github.com/user-attachments/assets/d932555d-5950-44b0-9b0f-4df2fafa3cb0" />

    

---

## 📅 Pan-genome Analysis
**Date:** 2025. 11. 13 & 11. 20

### 🎯 Goal
The pan-genome was analyzed to visualize the population structure and to identify ST-2–specific genetic markers

### 🛠️ Task
1.  **Pan-genomic statistical overview (Roary)**
    * **Core Genes:** 1,522 (15%) / **Accessory Genes:** 8,492 (85%)
    
2.  **Clustering**
    * **Analysis:** Dendrogram Based on Gene Presence/Absence Matrix
    * **Result:** **The dataset segregated into two major clades: a **Red Clade** centered on **ST-2** and a **Blue Clade** centered on **ST-5**
      <img width="1114" height="1125" alt="image" src="https://github.com/user-attachments/assets/4fe29fae-11e2-4561-9178-dd6416e62cc7" />


      Notably, some strains were positioned distantly from their corresponding ST groups, indicating intra-ST genomic heterogeneity
3.  **Genome Wide Association Study (GWAS):**
    * **Tool:** **Scoary**
    * **Finding:** dentified a marker gene, **group_4200 (phage-derived protein)**, that distinguishes the ST-2 group (P < 1.0 × 10⁻¹³)

---

## 📅 Functional Profiling
**Date:** 2025. 11. 27 (Refined Analysis)

### 🎯 Goal
Explored key survival strategies in pathogenic strains by examining auxotrophy mechanisms and profiling the current state of antimicrobial resistance (AMR)

### 🛠️ Task
1.  **Auxotrophy**
    * **Method:** Pathway completeness for amino acid biosynthesis was evaluated based on **GapMind/KEGG**
    * **Key Finding:**
        * **Auxotroph:** **Val**, **Ile**, **His**, **Phe**, **Tyr**, **Ser**, **Asn**
        * **Mechanism:** ilvH, hisE, and tyrB genes were found to be absent across all 124 strains
   1-2. **Visualization**
      <img width="930" height="1125" alt="image" src="https://github.com/user-attachments/assets/5403d6b1-f76b-41b8-b7fb-c164ab7df67f" />

        
        However, discrepancies were observed when compared with results from previous experimental auxotrophy assays conducted on *S. epidermidis*, indicating that validation between the two approaches is required

2.  **AMR & Virulence**
    * **Tool:** **AMRFinderPlus**
    * **Key Finding:**
        * **`hld` (Toxin):** Present in 100% of strains
        * **`mecA` (Methicillin):** Present in 50% of strains, matching precisely with the ST-2 group
        * **`fosB`, `blaZ`:** Detected at high frequencies
   2-2 **Visualization**
<img width="1232" height="1125" alt="image" src="https://github.com/user-attachments/assets/86f379d1-a265-420a-9cf8-7633b1f035f2" />


---

## 📅 Advanced Profiling
**Date:** 2025. 12. 04

### 🎯 Goal
**MLST profiling(wgMLST)** and additional functional analyses

### 🛠️ Task
1.  **wgMLST (High-Resolution Phylogeny)**
    * **Tool:** **chewBBACA** -> **GrapeTree**
    * **Topology:** **Star-burst Pattern** (Sea urchin shape)
      <img width="927" height="682" alt="image" src="https://github.com/user-attachments/assets/dbe385e8-bd82-474e-b6e7-42428be7a088" />

    * **Conclusion:** The ST-2 clone originated from a single ancestor and underwent rapid clonal expansion within the hospital environment
2.  **BGC Profiling (Secondary Metabolites)**
    * **Tool:** **antiSMASH 7.0**
    * **Result:** Siderophores are essential, but bacteriocins—key competitive weapons—were not detected within the ST-2 group
3.  **Risk Assessment**
    * **Method:** A Composite Risk Index was calculated by summing the number of virulence factors and AMR genes in each strain
    *  **Result:** The MRSE group exhibited a significantly higher risk level than the MSSE group (P < 0.001)
    * **Visualization:** **Boxplot with Strip plot**
   <img width="2400" height="1800" alt="Figure_Risk_Boxplot" src="https://github.com/user-attachments/assets/561b85bb-2a40-4f7e-90cf-a9d00f58ee76" />
   
4.  **iTOL Integration**
    * **Combined Plot:** Tree + Heatmap + Bar chart
   
   Strains with a red background represent MRSE, while those with a green background represent MSSE. Strains harboring bacteriocins are indicated by blue rings. The virulence index is displayed in the form of bar graphs.
Although no strains in the ST-2 group carried bacteriocin genes, they exhibited an average risk score above 13, indicating that they are potentially hazardous.
    <img width="1163" height="1089" alt="image" src="https://github.com/user-attachments/assets/1bc07aa0-6045-476b-a755-03351934f1d6" />

    

---

To be added
