# Lactobacillus fermentum Pan-Genome Analysis Report  
### **Defining L. fermentum Phylogroups Through Accessory Gene Clustering and Functional COG Profiling**

---

# 📌 1. Introduction

*Lactobacillus fermentum*은 유제품, 채소 발효식품, 인체 장·구강 등 다양한 발효 환경에서 반복적으로 분리되며, 대체로 **1.8–2.3 Mb** 크기의 유전체를 가진다(Azmal Ali et al., 2018).  
동일 종 내에서도 균주 간 탄수화물 이용성, 에너지 전환, 아미노산 대사 및 스트레스 내성과 같은 **발효 관련 phenotype과 대응 유전자 조성에서 큰 변이**가 보고되어 왔다(Konyanee et al., 2019).

이러한 기능적 변이의 상당 부분은 보존적인 **core genome**만으로는 설명하기 어렵고, 환경 적응과 niche 특이성을 반영하는 **accessory genome 조합 차이**에서 기인한다는 점이 지적되어 왔다(Tatusov et al., 2001).
그럼에도 불구하고 기존 *L. fermentum* 유전체 연구는 주로 ANI 또는 core genome 기반 계통 분석에 집중되어 있으며,accessory genome 패턴과 발효 관련 기능을 통합하여 phylogroup을 정의한 체계적 분석은 제한적이다.

따라서 본 연구에서는 공개된 *L. fermentum* 유전체 **155개**를 ANI 기준으로 필터링하여 **80개의 high-quality genome**을 선별한 뒤 pan-genome을 구축하고, accessory gene의 존재/부재 패턴을 기반으로 phylogroup을 정의하였다.  
이어 각 균주의 accessory 유전자에 **COG 기반 기능 주석**을 부여하고, 발효와 직접적으로 연관된 **C (Energy), E (Amino acid), G (Carbohydrate)** 카테고리의 조성을 정량화하여 phylogroup 간 발효 대사 시그니처를 비교하였다.

또한 동일 phylogroup 내에서 관찰되는 **outlier strain**의 경우, 발효와 관련된 특정 COG 기능 유전자의 소실 또는 획득이해당 균주가 다른 위치에 분리되는 현상을 설명할 수 있는지를 분석함으로써,
L. fermentum phylogroup이 단순한 유전체 유사도가 아니라 발효 기능 차이를 반영한 기능적 그룹임을 확인하고자 하였다.

<img width="2000" height="1080" alt="introduction" src="https://github.com/user-attachments/assets/1ba507ea-6724-4b91-9e32-7d4ba39389f8" />

> **Figure 1.** Overview of Functional Divergence and Phylogroup Structure in L. fermentum  

---

# 🧬 2. Materials and Methods

## 2.1 Genome Sequencing Data Acquisition

NCBI RefSeq 데이터베이스에서 *Lactobacillus fermentum*으로 분류된 모든 유전체(총 **155개**)를 FASTA 형식으로 다운로드 하였다. Genome assembly level이 “Complete Genome”인 시퀀스만 1차 선별하였으며, 이후 fastANI를 이용해 균주 간 쌍별 ANI를 계산하였다.

Figure 2.와 같이 ANI 95% 미만 값을 보이는 종 내 이질적 균주 (예: AP017974.1 등)와 중복 유전체를 제거하고, complete한 genome 수준의 유전체만을 남겨 최종 **81개 균주**를 downstream 분석에 사용하였다.

<img width="3200" height="2800" alt="ani_heatmap" src="https://github.com/user-attachments/assets/12dec12e-fa07-4ee5-b7c0-6763e4dad9fe" />

> **Figure 2.** ANI heatmap

**Workflow**  
- **Input**:  
  - [NCBI RefSeq *L. fermentum* genomes (FASTA)](https://www.ncbi.nlm.nih.gov/refseq/)
- **Tool**:  
  - [fastANI](https://github.com/ParBLiSS/FastANI)
- **Output**:   
  - [ANI_matrix.tsv](https://github.com/minjBaek02/file_collection/blob/5fa516d39e5c16e8ae79d74fc3bd6a41baf088dc/%EA%B9%83%ED%97%88%EB%B8%8C%20%EC%97%85%EB%A1%9C%EB%93%9C%20%ED%8C%8C%EC%9D%BC/fastani_lactobacillus_fermentum_.matrix)

---

## 2.2 Genome Annotation and Pan-Genome Construction

각 균주의 주석 형식을 통일하기 위해 Prokka (v1.14.6)를 사용하여 각 균주별 GFF 파일과 단백질 서열(FAA)을 생성하였다.

생성된 Prokka GFF 전체를 Roary (v3.13.0)에 입력하여 pan-genome을 구축하였으며, gene_presence_absence.* 및 gene_presence_absence.Rtab 파일을 통해 core 및 accessory 유전자 집합을 정의하였다.  
전체 균주의 **95% 이상에서 발견되는 유전자**를 core genome으로, 그 외(shell + cloud)를 accessory genome으로 정의하였다.

Accessory genome 기반 계통수에서 비정상적으로 긴 branch를 형성한 **CP033371.1**을 추가로 제거하였고,  
최종적으로 **80개 균주에 대한 presence/absence matrix**를 확정하였다.

**Workflow**  
- **Input**:  
  - [Final 81 genome FASTA]
- **Tool**:  
  - [Prokka v1.14.6](https://github.com/tseemann/prokka)  
  - [Roary v3.13.0](https://github.com/sanger-pathogens/Roary)
- **Output**:  
  - [Prokka GFF / FAA]
  - [gene_presence_absence.csv](https://github.com/minjBaek02/file_collection/blob/90180ecefe398d81ddba80568d430fc4c187ce96/%EA%B9%83%ED%97%88%EB%B8%8C%20%EC%97%85%EB%A1%9C%EB%93%9C%20%ED%8C%8C%EC%9D%BC/Galaxy1052-%5BRoary%20on%20data%201033%2C%20data%201021%2C%20and%20others%20Gene%20Presence%20Absence%5D.csv)  
  - [accessory_matrix_80strains.tsv](https://github.com/minjBaek02/file_collection/blob/d6fbe99682eeba875db3a349534f11ef114b4562/%EA%B9%83%ED%97%88%EB%B8%8C%20%EC%97%85%EB%A1%9C%EB%93%9C%20%ED%8C%8C%EC%9D%BC/accessory_presence_absence_renamed.tsv)

---

## 2.3 Phylogrouping 

Accessory genome 기반 phylogroup 정의를 위해 Roary에서 생성된 accessory gene 계통수 (`accessory_tree.nwk`)를 사용하였다.

계통수 기반 거리 정보를 이용해 계층적 군집화를 수행하고, **cut height = 0.35** 기준으로 phylogroup을 정의하였다.  
최종 phylogroup 정보는 [`phylogroup_h35.tsv`](https://github.com/minjBaek02/file_collection/blob/90180ecefe398d81ddba80568d430fc4c187ce96/%EA%B9%83%ED%97%88%EB%B8%8C%20%EC%97%85%EB%A1%9C%EB%93%9C%20%ED%8C%8C%EC%9D%BC/phylogroup_h35.tsv)로 정리하였다.

**Workflow**  
- **Input**:  
  - [accessory_tree.nwk](https://github.com/minjBaek02/file_collection/blob/90180ecefe398d81ddba80568d430fc4c187ce96/%EA%B9%83%ED%97%88%EB%B8%8C%20%EC%97%85%EB%A1%9C%EB%93%9C%20%ED%8C%8C%EC%9D%BC/Accessory%20Binary%20Genes%20.Newick.nhx)
- **Tool**:
  - [iTOL (Interactive Tree Of Life)](https://itol.embl.de/)
- **Output**:  
  - [phylogroup_h35.tsv](https://github.com/minjBaek02/file_collection/blob/90180ecefe398d81ddba80568d430fc4c187ce96/%EA%B9%83%ED%97%88%EB%B8%8C%20%EC%97%85%EB%A1%9C%EB%93%9C%20%ED%8C%8C%EC%9D%BC/phylogroup_h35.tsv) 

Accessory genome 기반 phylogrouping을 수행하기 전에, CP033371.1이 트리에서 지나치게 긴 단독 branch를 형성하는 extreme outlier임을 확인하였다. (ANI 분석에서는 이상 없음, Roary 단계에서 gene clustering 혹은 annotation mismatch 가능성)

이로 인해 거리 기반 분석(H35 분할 포함)을 왜곡할 정도의 비정상적 값 생성하기에 downstream 분석의 정확성을 위해
CP033371.1을 제외하고 총 80 strain 기준으로 H35 phylogrouping 및 후속 C/G/E 분석을 수행하였다.

---

## 2.4 Functional Categorization (COG Profiling)

Roary로부터 얻은 pan-genome 단백질 서열 집합(또는 Prokka FAA 전체 병합본)을 eggNOG-mapper v2.1.9에 입력하여 eggNOG 5.0 데이터베이스 기반의 기능 주석과 COG 카테고리를 할당하였다. ('egg-nog-result.tabular')

이 중 에너지 대사 [C], 탄수화물 대사 및 수송 [G], 아미노산 대사 및 수송 [E] 세 카테고리를 선택하여  
strain × C/G/E 카운트 및 비율 매트릭스를 구성하였다.

CP033371.1을 제외한 **80개 균주**에 대해 `cog_by_strain_CGE_h35_no_CP033371.tsv` 파일을 생성하여  
phylogroup별 C/G/E 조성 비교 및 시각화에 사용하였다.

**Workflow**  
- **Input**:  
  - [Pan-genome FAA](data/prokka/)
- **Tool**:  
  - [eggNOG-mapper v2.1.9](https://github.com/eggnogdb/eggnog-mapper)
- **Output**:
  -  [cog_by_strain_CGE_h35_no_CP033371.tsv](https://github.com/minjBaek02/file_collection/blob/90180ecefe398d81ddba80568d430fc4c187ce96/%EA%B9%83%ED%97%88%EB%B8%8C%20%EC%97%85%EB%A1%9C%EB%93%9C%20%ED%8C%8C%EC%9D%BC/cog_by_strain_CGE_h35_no_CP033371.tsv)
  -  [egg-nog-result.tabular](https://drive.google.com/file/d/1QG7--Y-G5xNm89FKU1c2Bf3uvwBT9V7K/view?usp=drive_link)

---

## 2.5 Statistics

Phylogroup 간 기능 차이를 평가하기 위해 각 strain의 C/G/E 유전자 수 및 비율을 phylogroup 정보와 결합한 뒤,  
카이제곱 검정과 Kruskal–Wallis 검정을 통해 phylogroup 간 분포 차이의 유의성을 평가하였다.  
또한 유전체 크기 및 총 유전자 수 차이에 따른 편향을 보정하기 위해 Fisher의 정확 검정을 이용한 enrichment 분석을 추가로 수행하였다.

마지막으로 strain × C/G/E 비율 매트릭스를 입력으로 PCA 분석을 수행하여 주요 성분(PC1, PC2) 상에서 strain을 투영하고,  
phylogroup별 색상 구분을 통해 기능적 프로파일의 군집 패턴을 시각화하였다.

**Workflow**  
- **Input**:   
   -  [cog_by_strain_CGE_h35_no_CP033371.tsv](https://github.com/minjBaek02/file_collection/blob/90180ecefe398d81ddba80568d430fc4c187ce96/%EA%B9%83%ED%97%88%EB%B8%8C%20%EC%97%85%EB%A1%9C%EB%93%9C%20%ED%8C%8C%EC%9D%BC/cog_by_strain_CGE_h35_no_CP033371.tsv)
   -  [phylogroup_h35.tsv](https://github.com/minjBaek02/file_collection/blob/90180ecefe398d81ddba80568d430fc4c187ce96/%EA%B9%83%ED%97%88%EB%B8%8C%20%EC%97%85%EB%A1%9C%EB%93%9C%20%ED%8C%8C%EC%9D%BC/phylogroup_h35.tsv) 
- **Tool**:  
  - [Biopython](https://biopython.org/)  
  - [SciPy](https://scipy.org/)  
  - [scikit-learn](https://scikit-learn.org/)
- **Output**:  
  - [CGE_statistics.tsv](https://github.com/minjBaek02/file_collection/blob/7921a4415cf2962ba541b6143adc52028ba16572/%EA%B9%83%ED%97%88%EB%B8%8C%20%EC%97%85%EB%A1%9C%EB%93%9C%20%ED%8C%8C%EC%9D%BC/C_G_E_statistical_results.tsv)  
  - [PCA_scores.csv](https://github.com/minjBaek02/file_collection/blob/7921a4415cf2962ba541b6143adc52028ba16572/%EA%B9%83%ED%97%88%EB%B8%8C%20%EC%97%85%EB%A1%9C%EB%93%9C%20%ED%8C%8C%EC%9D%BC/CGE_PCA_scores.csv)
and so on...

---

# 📊 3. Results

## Results I. Accessory Genome-based Phylogrouping and Tree Structure

### 3.1 Phylogenetic Tree Based on Accessory Genome
Accessory 유전자의 존재/부재 패턴을 기반으로 구축한 accessory genome 계통수에서, patristic distance를 이용한 계층적 군집화 결과 *Lactobacillus fermentum* 균주들은 **cut height = 0.35 (H35)**에서 **다섯 개의 주요 phylogroup**으로 명확하게 분리되었다.  

따라서, method 2.4와 같이 phylogenetic tree를 구축하고, 각 strain에 대해 COG C/E/G 카테고리 유전자 수를 매핑하여  
**계통 구조와 발효 기능 조성 간의 연관성**을 평가하였다.

<img width="1000" height="800" alt="image" src="https://github.com/user-attachments/assets/234d12b6-65d1-47e3-a589-7c7ecc681bea" />

> **Figure 3.** Accessory gene presence/absence 기반 h 0.35 cut phylogenetic tree.  
> 각 strain 옆에는 COG C (에너지 대사), G (탄수화물 대사/수송), E (아미노산 대사/수송) 카테고리 유전자 수가 함께 시각화되어 있다.

#### 🔎 Key Observations

1. **CP076082.1은 GROUP1에 속하지만, 내부에서도 가장 외곽(branch tip)에 위치하는 outlier strain으로 확인됨**  
   → accessory genome 조성에서 GROUP1 평균 패턴과의 이탈 가능성 시사  

2. Tree 상의 서로 다른 구역에서 **C/G/E functional profile이 일관되게 다른 패턴**을 보임  
   → phylogroup 분리가 단순한 계통학적 거리뿐 아니라 **기능적 조성 차이와 연관**되어 있을 가능성 제기  

이러한 관찰이 단순한 계통 구조상의 효과인지, 혹은 실제로 **phylogroup 간 기능 조성 차이가 통계적으로 유의한지**를 검증하기 위해  후속 기능 비교 및 통계 분석을 수행하였다.

---

## Results II. Functional Divergence Among Phylogroups

### 3.3 Statistical Validation of Functional Divergence

#### 3.3.1 Group-level Total Counts (Chi-square Test)

Phylogroup 간 C/G/E 카테고리별 **총 유전자 수 분포**를 비교하기 위해 카이제곱 검정을 수행한 결과, 세 기능 카테고리 모두에서 극도로 유의한 차이가 관찰되었다.

- Energy metabolism (C): *p* ≈ 1.2 × 10⁻²⁷⁶  
- Carbohydrate metabolism/transport (G): *p* ≈ 5.1 × 10⁻²¹⁷  
- Amino acid metabolism/transport (E): *p* < 0.001  

<img width="40%" alt="cog_group_counts_h35_heatmap" src="https://github.com/user-attachments/assets/d024d185-b164-4468-a3b6-d5caefd101eb" />

> **Figure 4A.** Phylogroup × C/G/E 총 유전자 수 count heatmap

이는 관찰된 기능 조성 차이가 우연에 의한 변동이 아니라, 각 phylogroup의 진화 과정에서 고착된 **구조적 기능 차이**임을 강하게 시사한다

---

#### 3.3.2 Strain-level Distribution (Kruskal–Wallis Test)

개별 균주 수준에서 C/G/E 유전자 수 분포를 비교하기 위해 비모수 Kruskal–Wallis 검정을 적용한 결과, 모든 카테고리에서 phylogroup 간 분포 차이가 통계적으로 유의하였다.

- C: *p* ≈ 6.4 × 10⁻⁶  
- G: *p* ≈ 3.6 × 10⁻⁵  
- E: *p* ≈ 3.7 × 10⁻¹¹

<img
  width="1560"
  height="960"
  alt="cog_CGE_h35_boxplot_smallfliers"
  src="https://raw.githubusercontent.com/minjBaek02/file_collection/6d92eb33ab26ca0789dd45a316c029d4aec6d0d7/깃허브%20업로드%20파일/cog_CGE_h35_boxplot_smallfliers.png"
/>
  
> **Figure 4B.** Phylogroup별 strain-level C/G/E 분포 (boxplot)

특히 E 카테고리에서 phylogroup 간 분포 분리가 가장 뚜렷하게 나타났으며, 이는 strain별 C/G/E 비율을 나타낸 boxplot에서  
**Group 1과 저기능 phylogroup 간 중앙값 차이**로 명확히 확인된다.

---

#### 3.3.3 Genome-size Adjusted Enrichment (Fisher’s Exact Test)

유전체 크기 차이에 따른 편향을 보정하기 위해, 각 phylogroup에서 C/G/E 카테고리가 차지하는 비율을 기준으로 2×2 contingency table을 구성하고 Fisher’s exact test를 수행하였다.

<img width="1560" height="960" alt="fisher_CGE_enrichment_h35_barplot" src="https://github.com/user-attachments/assets/daf6f4c0-8f05-49d4-9899-ae4fc93d1a64" />

> **Figure 4C.** Genome-size 보정 후 phylogroup별 C/G/E enrichment 분석 결과

그 결과 **Group 1**은 C/G/E 기능 유전자에 대해 **Odds Ratio = 1.297** (*p* = 4.9 × 10⁻²⁰)을 보여, 전체 유전자 수를 고려하더라도 발효 관련 기능이 실제로 **유의하게 농축(enriched)**된 phylogroup으로 나타났다.

반대로 **Groups 2, 4, 5**는 동일 분석에서 C/G/E 카테고리가 유의하게 결손(depleted)된 패턴을 보였다.
→ 기능 차이는 genome size artifact가 아닌 **구조적 특성**임을 확인하였다.

---

## 🧬Results III. Deep-dive into an Outlier Strain (CP076082.1)

### 4.1 Anatomical Analysis of the Outlier: A Gap Within Group 1

Figure 3.과 같이 **CP076082.1**은 Group 1의 하위 클러스터에 속함에도 불구하고, 동일 phylogroup 내 다른 균주들과 비교하여 **비정상적으로 긴 가지(long branch)**를 형성하며 계통수의 외곽에 위치하였다.

흥미롭게도 ANI (Average Nucleotide Identity) 분석에서는 CP076082.1이 Group 1 균주들과 **높은 유전체 유사도**를 유지하고 있음이 확인되었으나, **accessory gene presence/absence 패턴에서는 뚜렷하게 이탈된 조성**을 보였다.

이는 CP076082.1이 비교적 최근의 진화 과정에서 **급격한 악세사리 유전자 조성 변화 또는 유전체 재구성(genomic rearrangement)**을 경험했을 가능성이 있어보였다.

---

### 4.2  Functional Divergence Analysis of CP076082.1

CP076082.1의 계통적 이탈 원인을 규명하기 위해, C/G/E 카테고리별 유전자 수를 기반으로 **Z-score 분석 및 PCA**를 수행하였다.  

#### 4.2.1 Functional Boxplot Analysis  
각 phylogroup(Group 1) 내에서 C/G/E 기능군 유전자 수의 분포를 boxplot으로 나타내고, 그 위에 CP076082.1의 상대적 위치를 z-score로 정규화하여 나타냈다. 

<img width="3000" height="1800" alt="CP0760821_CGE_boxplot_zscore" src="https://github.com/user-attachments/assets/1946ddf9-9e7f-4455-86d0-dd30c4215e11" />

> **Figure 5A.** Z-score–normalized boxplots of C/G/E functional gene counts within phylogroup Group 1.

- C (에너지 대사) 및 G (탄수화물 대사/수송) 유전자 수는  Group 1의 분포 범위 내에 안정적으로 위치하였다.
- 반면, **E 카테고리는 Group 1 내에서 유일하게 하위 5% 미만**에 해당하는 급격한 감소를 보였다.

이는 CP076082.1의 기능적 이탈이 전반적인 대사 붕괴가 아닌, **아미노산 대사 기능에 선택적으로 집중된 현상**임을 의미한다.

#### 4.2.2 Principal Component Analysis (PCA)  
또한, PCA 분석으로 C/G/E를 동시에 고려하며, strain들이 이 조합으로 어떻게 퍼지는지 시각화하였다.

<img width="2400" height="1800" alt="GROUP1_only_CGE_PCA_plot" src="https://github.com/user-attachments/assets/71875945-e470-49ef-a1cd-72e3392d0920" />

> **Figure 5B.** PCA of C/G/E Functional Profiles Reveals Within-Phylogroup Divergence in L. fermentum Group 1

PCA plot 상에서 CP076082.1을 Group 1의 중심부에서 분리시키는 주요 loading factor는 **E 카테고리 유전자 결손**으로 확인되었다.

즉, 해당 균주는 Group 1이 공유하는 기본적인 대사 골격(C/G)은 유지하면서도, **아미노산 대사 기능만 선택적으로 약화된 ‘기능적 변이체(functional variant)’**로 해석된다.

---

### 4.3 Genomic Evidence: Strategy of “Loss and Gain”

단순한 유전자 수 차이를 넘어, 실제로 어떤 대사 경로(pathway)가 변화했는지를 분석한 결과 CP076082.1의 유전체에서는 **‘대사적 최적화(metabolic optimization)’ 전략**이 뚜렷하게 관찰되었다.

#### 4.3.1 Major Pathway Loss: The Cost of Adaptation

Group 1 균주들이 공통적으로 보유하는 accessory gene 중, CP076082.1에서만 선택적으로 소실된 주요 경로는 다음과 같다.

- **Histidine biosynthesis operon (hisA/B/C/D/G/H/Z/K)**  
  → 히스티딘 생합성에 필요한 전체 오페론이 완전히 결손됨
- **Methionine metabolism genes (metI, metB, metE)**  
  → 메티오닌 생합성 및 재생 경로 약화
- **Translational factors (lepA, tuf)**  
  → 단백질 합성 효율과 관련된 보조 인자의 결손

이러한 유전자 소실 패턴은, 아미노산이 풍부한 환경에 적응한 균주가 **내생적 아미노산 합성 능력을 포기하는 대신 에너지 비용을 절감한** 전형적인 *genome streamlining* 사례로 해석될 수 있을 듯하다.

#### 4.3.2 Unique Gene Gain: Niche Specialization

반면, CP076082.1에서만 관찰되는 독특한 accessory gene 획득은 **특정 탄소원 활용 능력의 강화**를 시사한다.

- **malK**: Maltose transporter  
- **ugpA**: Glycerol-3-phosphate transporter  
- **araQ**: Arabinose transporter / regulator  

이는 해당 균주가 아미노산 자급 능력을 축소하는 대신,  **외부 탄수화물 자원을 보다 효율적으로 이용하는 방향으로 특화**되었음을 보여준다.

---

### 4.4 Summary of Results II

종합하면, CP076082.1의 계통적 이탈은 무작위적인 유전자 손실이 아니라,

> **아미노산 자급 능력(E)을 포기하는 대신,**  
> **특정 탄수화물 자원(G)을 효율적으로 활용하기 위한**  
> **전략적 유전체 재구성의 결과**

로 해석된다.

이러한 선택적 기능 재편성은 accessory genome presence/absence matrix에 직접 반영되었으며, 그 결과 계통 트리 상에서 **독립적인 장가지(long-branch) 분기**를 형성하게 된 것이다.

---

# 🧩 4. Integrated Conclusion

본 연구는 크게 두가지 축으로 분석을 진행했다.  
① **H35 phylogroup 간 C/G/E 기능 조성 차이**,  
② **동일 그룹 내 outlier(CP076082.1)의 기능·구조 이탈 분석** — 을 종합하면 다음과 같다.

## 핵심 결론

1. **C/G/E 기능군은 L. fermentum의 발효 대사를 규정하는 핵심 functional axis이며,**  
   phylogroup(H35) 간 이 축에서 통계적으로 유의한 차이가 존재한다.  
   → 즉, accessory genome 기반 phylogroup은 단순 gene 조성 차이가 아니라  
     **각 그룹의 고유한 발효 대사 signature**를 반영한다.

2. **CP076082.1은 E category 핵심 경로 결손 + 특이 carbon transporter 획득이라는 재구성을 통해**  
   동일 그룹 내에서도 기능적으로 이탈된 패턴을 보였다.  
   → 즉, **군집 분리와 outlier 형성 모두가 C/G/E 기능축의 구조적 변화와 직접적으로 연결된다.**
 
## ⭐ Final Statement

**“Accessory genome 기반 phylogroup은 L. fermentum의 발효 대사 signature를 명확히 구분하는 지표이며,  
CP076082.1과 같은 outlier의 분리 또한 C/G/E 기능축 재구성으로 설명된다.”**
