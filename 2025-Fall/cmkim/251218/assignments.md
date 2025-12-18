# Subcellular Localization Prediction in *Staphylococcus epidermidis*

## Overview
- **Genomes analyzed:** 237 *S. epidermidis* complete genomes  
- **Total genes:** 550,864  
- **Clustering method:** MMseqs2  
  - Sequence identity: **0.99**
  - Coverage: **1.0**
- **Clustered genes obtained:** 24,275  

### Gene cluster categories
- **Core:** presence > 0.95  
- **Accessory:** 0.95 > presence > 0.15  
- **Rare:** presence < 0.15

| Category | Definition | Counts |
|----------|-----------|-----------|
| Core | Present in > 95% of genomes | 673 |
| Accessory | Present in 15–95% of genomes | 3,085 |
| Rare | Present in < 15% of genomes | 20,517 |  

### Gene clustering reuslt  
- **Color:** Core, blue; accessory, orange; rare, green
<img width="80%" height="60%" alt="Clustering" src="https://github.com/user-attachments/assets/48f28e7a-82db-442f-ab7a-a520db5e00f1" />  

---

## Protein localization prediction
Representative sequences from each gene cluster were analyzed using:
- **SignalP** (signal peptide / secretion signal)
- **CWpred** (LPXTG motif or derivatives)
- **TMHMM** (number of predicted transmembrane helices)

The three prediction results were integrated to infer **subcellular localization**.

---

## Integrated localization rules

| Location prediction | SignalP prediction | CWpred prediction type | TMHMM (PredHel) |
|---------------|-------------------|------------------------|-----------------|
| cell wall | SP or LIPO or PILIN | Cell-wall binding protein | ≤ 2 |
| cell wall | SP or LIPO or PILIN | Cell-wall binding protein | > 2 |
| cell wall | OTHER | Cell-wall binding protein | ≤ 2 |
| cell wall | OTHER | Cell-wall binding protein | > 2 |
| secretion | SP or LIPO or PILIN | Globular or Transmembrane protein | ≤ 2 |
| cell membrane | SP or LIPO or PILIN | Globular or Transmembrane protein | > 2 |
| cell membrane | OTHER | Globular or Transmembrane protein | > 2 |
| intracellular | OTHER | Globular or Transmembrane protein | ≤ 2 |  
  
### Predicted subcellular localization counts (% by row-wise ratio)

| MMseqs2 PanClassification | cellwall | intracellular | membrane | secretion |
|---------------------------|----------|---------------|----------|-----------|
| Accessory | 15 (0.49) | 2,441 (79.12) | 483 (15.66) | 146 (4.73) |
| Core | 0 (0) | 541 (80.39) | 119 (17.68) | 13 (1.93) |
| Rare | 474 (2.31) | 16,121 (78.57) | 2,587 (12.61) | 1,335 (6.51) |  
  
### Stacked bar plot showing gene localization prediction by pan-genome classification

<img width="50%" height="30%" alt="PanClassificationByLocation" src="https://github.com/user-attachments/assets/6f88f0c1-01e8-44c7-a201-b66d12f4344f" />
<img width="45%" height="30%" alt="LocationByPanClassficiation" src="https://github.com/user-attachments/assets/716dd20a-db51-49be-b630-5a8f2c354a8f" />


---
## Statistical analysis

A chi-square test of independence was performed to assess the association between  
**MMseqs2 pan-genome classification** and **predicted subcellular localization**.

- **Chi-square statistic:** 122.21  
- **P-value:** 5.60 × 10⁻²⁴

To identify **which subcellular localizations contribute most strongly** to the observed chi-square association,  
**standardized residuals** were calculated for each combination of pan-genome class and predicted location.

### Standardized residuals table

| MMseqs2 PanClassification | cellwall | intracellular | membrane | secretion |
|---------------------------|----------|---------------|----------|-----------|
| Accessory | -5.98 | 0.27 | 3.86 | -3.18 |
| Core | -3.68 | 0.49 | 3.25 | -4.41 |
| Rare | 2.99 | -0.19 | -2.09 | 2.03 |

- visualzed by heatmap  
<img width="60%" height="40%" alt="chi_residual_HM" src="https://github.com/user-attachments/assets/76e680b9-0ab1-4f79-ac6d-5ed1e1688244" />





### Conclusion
The chi-square signal is primarily driven by:
- Enrichment of **surface-associated proteins (cellwall, secretion)** in the **rare gene pool**.
- Enrichment of **membrane proteins** in **core and accessory genes**.

This pattern supports the biological interpretation that **rare genes are more likely to encode host-interacting or surface-exposed functions**, whereas **core genes are biased toward conserved intracellular and membrane functions**.
