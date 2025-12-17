# 📌 1. Introduction

*Lactobacillus fermentum*은 유제품, 채소 발효식품, 인체 장·구강 등 다양한 발효 환경에서 반복적으로 분리되며, 대체로 **1.8–2.3 Mb** 크기의 유전체를 가진다(Azmal Ali et al., 2018).  
동일 종 내에서도 균주 간 탄수화물 이용성, 에너지 전환, 아미노산 대사 및 스트레스 내성과 같은 **발효 관련 phenotype과 대응 유전자 조성에서 큰 변이**가 보고되어 왔다(Konyanee et al., 2019).

이러한 기능적 변이의 상당 부분은 보존적인 **core genome**만으로는 설명하기 어렵고, 환경 적응과 niche 특이성을 반영하는 **accessory genome 조합 차이**에서 기인한다는 점이 지적되어 왔다(Tatusov et al., 2001).
그럼에도 불구하고 기존 *L. fermentum* 유전체 연구는 주로 ANI 또는 core genome 기반 계통 분석에 집중되어 있으며,accessory genome 패턴과 발효 관련 기능을 통합하여 phylogroup을 정의한 체계적 분석은 제한적이다.

따라서 본 연구에서는 공개된 *L. fermentum* 유전체 **155개**를 ANI 기준으로 필터링하여 **80개의 high-quality genome**을 선별한 뒤 pan-genome을 구축하고, accessory gene의 존재/부재 패턴을 기반으로 phylogroup을 정의하였다.  
이어 각 균주의 accessory 유전자에 **COG 기반 기능 주석**을 부여하고, 발효와 직접적으로 연관된 **C (Energy), E (Amino acid), G (Carbohydrate)** 카테고리의 조성을 정량화하여 phylogroup 간 발효 대사 시그니처를 비교하였다.

또한 동일 phylogroup 내에서 관찰되는 **outlier strain**의 경우, 발효와 관련된 특정 COG 기능 유전자의 소실 또는 획득이해당 균주가 다른 위치에 분리되는 현상을 설명할 수 있는지를 분석함으로써,
L. fermentum phylogroup이 단순한 유전체 유사도가 아니라 발효 기능 차이를 반영한 기능적 그룹임을 확인하고자 하였다.

---

# 🧬 2. Materials and Methods

## 2.1 Genome Sequencing Data Acquisition

NCBI RefSeq 데이터베이스에서 *Lactobacillus fermentum*으로 분류된 모든 유전체(총 **155개**)를 FASTA 형식으로 다운로드 하였다. Genome assembly level이 “Complete Genome”인 시퀀스만 1차 선별하였으며, 이후 fastANI를 이용해 균주 간 쌍별 ANI를 계산하였다.

ANI 95% 미만 값을 보이는 종 내 이질적 균주 (예: AP017974.1 등)는 분석에서 제외하였고, 중복 및 저품질 유전체를 제거한 후  
최종 **81개 균주**를 downstream 분석에 사용하였다.

**Workflow**  
- **Input**:  
  - [NCBI RefSeq *L. fermentum* genomes (FASTA)](https://www.ncbi.nlm.nih.gov/refseq/)
- **Tool**:  
  - [fastANI](https://github.com/ParBLiSS/FastANI)
- **Output**:  
  - [final_81_strains.tsv](data/metadata/final_81_strains.tsv)  
  - [ANI_matrix.tsv](data/metadata/ANI_matrix.tsv)

---

## 2.2 Genome Annotation and Pan-Genome Construction

각 균주의 주석 형식을 통일하기 위해 Prokka (v1.14.6)를 사용하여 ORF를 재예측하고 기능을 주석화하였다.  
각 균주별 GFF 파일과 단백질 서열(FAA)을 생성하였다.

생성된 Prokka GFF 전체를 Roary (v3.13.0)에 입력하여 pan-genome을 구축하였으며, gene_presence_absence.* 및 gene_presence_absence.Rtab 파일을 통해 core 및 accessory 유전자 집합을 정의하였다.  
전체 균주의 **95% 이상에서 발견되는 유전자**를 core genome으로, 그 외(shell + cloud)를 accessory genome으로 정의하였다.

Accessory genome 기반 계통수에서 비정상적으로 긴 branch를 형성한 **CP033371.1**을 추가로 제거하였고,  
최종적으로 **80개 균주에 대한 presence/absence matrix**를 확정하였다.

**Workflow**  
- **Input**:  
  - [Final 81 genome FASTA](data/genomes/)
- **Tool**:  
  - [Prokka v1.14.6](https://github.com/tseemann/prokka)  
  - [Roary v3.13.0](https://github.com/sanger-pathogens/Roary)
- **Output**:  
  - [Prokka GFF / FAA](data/prokka/)  
  - [gene_presence_absence.csv](data/roary/gene_presence_absence.csv)  
  - [accessory_matrix_80strains.tsv](data/roary/accessory_matrix_80strains.tsv)

---

## 2.3 Functional Categorization (COG Profiling)

Roary로부터 얻은 pan-genome 단백질 서열 집합(또는 Prokka FAA 전체 병합본)을 eggNOG-mapper v2.1.9에 입력하여 eggNOG 5.0 데이터베이스 기반의 기능 주석과 COG 카테고리를 할당하였다.

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
  - [cog_by_strain_CGE_h35_no_CP033371.tsv](data/cog/cog_by_strain_CGE_h35_no_CP033371.tsv)

---

## 2.4 Phylogrouping and Statistical Analysis

Accessory genome 기반 phylogroup 정의를 위해 Roary에서 생성된 accessory gene 계통수 (`accessory_tree.nwk`)를 사용하였다.

계통수 기반 거리 정보를 이용해 계층적 군집화를 수행하고, **cut height = 0.35** 기준으로 phylogroup을 정의하였다.  
최종 phylogroup 정보는 [`phylogroup_h35.tsv`](data/metadata/phylogroup_h35.tsv)로 정리하였다.

Phylogroup 간 기능 차이를 평가하기 위해 각 strain의 C/G/E 유전자 수 및 비율을 phylogroup 정보와 결합한 뒤,  
카이제곱 검정과 Kruskal–Wallis 검정을 통해 phylogroup 간 분포 차이의 유의성을 평가하였다.  
또한 유전체 크기 및 총 유전자 수 차이에 따른 편향을 보정하기 위해 Fisher의 정확 검정을 이용한 enrichment 분석을 추가로 수행하였다.

마지막으로 strain × C/G/E 비율 매트릭스를 입력으로 PCA 분석을 수행하여 주요 성분(PC1, PC2) 상에서 strain을 투영하고,  
phylogroup별 색상 구분을 통해 기능적 프로파일의 군집 패턴을 시각화하였다.

**Workflow**  
- **Input**:  
  - [accessory_tree.nwk](data/roary/accessory_tree.nwk)  
  - [cog_by_strain_CGE_h35_no_CP033371.tsv](data/cog/cog_by_strain_CGE_h35_no_CP033371.tsv)
- **Tool**:  
  - [Biopython](https://biopython.org/)  
  - [SciPy](https://scipy.org/)  
  - [scikit-learn](https://scikit-learn.org/)
- **Output**:  
  - [phylogroup_h35.tsv](data/metadata/phylogroup_h35.tsv)  
  - [CGE_statistics.tsv](results/statistics/CGE_statistics.tsv)  
  - [PCA_plot.png](figures/PCA_plot.png)

---

# 📊 3. Results

## Results I. Accessory Genome-based Phylogrouping and Tree Structure

### 3.1 Accessory Genome-based Phylogrouping (H35)

Accessory 유전자의 존재/부재 패턴을 기반으로 구축한 accessory genome 계통수에서,  
patristic distance를 이용한 계층적 군집화 결과 *Lactobacillus fermentum* 균주들은  
**다섯 개의 주요 phylogroup**으로 명확하게 분리되었다.

DynamicTreeCut 알고리즘을 적용하여 여러 cut height를 비교한 결과,  
**cut height = 0.35 (H35)**에서 군집 경계가 가장 안정적으로 재현되었으며,  
이후 모든 구조적·기능적 분석은 해당 phylogroup 레이블을 기준으로 수행하였다.

분석 초기 단계에서 **CP033371.1**은 accessory tree 상에서  
비정상적으로 긴 branch를 형성하는 extreme outlier로 확인되었고,  
patristic distance 구조를 왜곡할 가능성이 있어  
phylogroup 정의 및 후속 분석에서 제외하였다.

> **Figure 2.** Accessory genome 기반 계통수 및 H35 기준 phylogroup 분리 결과  
> (accessory_tree.nwk, phylogroup_h35.tsv)

---

### 3.2 Phylogenetic Tree Based on Accessory Genome

Accessory 유전자는 core genome에 비해 변동성이 크므로,  
phylogroup 정의에 사용되는 cut-off의 안정성이 중요하다.  
본 연구에서는 여러 cut height를 비교한 결과,  
**H35 threshold가 가장 안정적으로 phylogroup을 분해하는 기준**으로 판단되었다.

Accessory gene presence/absence matrix를 기반으로 phylogenetic tree를 구축하고,  
각 strain에 대해 COG C/E/G 카테고리 유전자 수를 매핑하여  
**계통 구조와 발효 기능 조성 간의 연관성**을 평가하였다.

> **Figure 3.** Accessory gene presence/absence 기반 phylogenetic tree.  
> 각 strain 옆에는 COG C (에너지 대사), G (탄수화물 대사/수송),  
> E (아미노산 대사/수송) 카테고리 유전자 수가 함께 시각화되어 있다.

#### 🔎 Key Observations

1. **CP076082.1은 GROUP1에 속하지만, 내부에서도 가장 외곽(branch tip)에 위치하는 outlier strain으로 확인됨**  
   → accessory genome 조성에서 GROUP1 평균 패턴과의 이탈 가능성 시사  

2. Tree 상의 서로 다른 구역에서 **C/G/E functional profile이 일관되게 다른 패턴**을 보임  
   → phylogroup 분리가 단순한 계통학적 거리뿐 아니라  
   → **기능적 조성 차이와 연관**되어 있을 가능성 제기  

이러한 관찰이  
단순한 계통 구조상의 효과인지,  
혹은 실제로 **phylogroup 간 기능 조성 차이가 통계적으로 유의한지**를 검증하기 위해  
후속 기능 비교 및 통계 분석을 수행하였다.

---

## Results II. Functional Divergence Among Phylogroups

### 3.3 Statistical Validation of Functional Divergence

#### 3.3.1 Group-level Total Counts (Chi-square Test)

Phylogroup 간 C/G/E 카테고리별 **총 유전자 수 분포**를 비교하기 위해  
카이제곱 검정을 수행한 결과,  
세 기능 카테고리 모두에서 극도로 유의한 차이가 관찰되었다.

- Energy metabolism (C): *p* ≈ 1.2 × 10⁻²⁷⁶  
- Carbohydrate metabolism/transport (G): *p* ≈ 5.1 × 10⁻²¹⁷  
- Amino acid metabolism/transport (E): *p* < 0.001  

이는 관찰된 기능 조성 차이가  
우연에 의한 변동이 아니라,  
각 phylogroup의 진화 과정에서 고착된  
**구조적 기능 차이**임을 강하게 시사한다.

> **Figure 4A.** Phylogroup × C/G/E 총 유전자 수 count heatmap

---

#### 3.3.2 Strain-level Distribution (Kruskal–Wallis Test)

개별 균주 수준에서 C/G/E 유전자 수 분포를 비교하기 위해  
비모수 Kruskal–Wallis 검정을 적용한 결과,  
모든 카테고리에서 phylogroup 간 분포 차이가 통계적으로 유의하였다.

- C: *p* ≈ 6.4 × 10⁻⁶  
- G: *p* ≈ 3.6 × 10⁻⁵  
- E: *p* ≈ 3.7 × 10⁻¹¹  

특히 E 카테고리에서  
phylogroup 간 분포 분리가 가장 뚜렷하게 나타났으며,  
이는 strain별 C/G/E 비율을 나타낸 boxplot에서  
**Group 1과 저기능 phylogroup 간 중앙값 차이**로 명확히 확인된다.

> **Figure 4B.** Phylogroup별 strain-level C/G/E 분포 (boxplot)

---

#### 3.3.3 Genome-size Adjusted Enrichment (Fisher’s Exact Test)

유전체 크기 차이에 따른 편향을 보정하기 위해,  
각 phylogroup에서 C/G/E 카테고리가 차지하는 비율을 기준으로  
2×2 contingency table을 구성하고  
Fisher’s exact test를 수행하였다.

그 결과 **Group 1**은  
C/G/E 기능 유전자에 대해 **Odds Ratio = 1.297**  
(*p* = 4.9 × 10⁻²⁰)을 보여,  
전체 유전자 수를 고려하더라도  
발효 관련 기능이 실제로 **유의하게 농축(enriched)**된 phylogroup으로 나타났다.

반대로 **Groups 2, 4, 5**는  
동일 분석에서 C/G/E 카테고리가  
유의하게 결손(depleted)된 패턴을 보였다.

> **Figure 4C.** Genome-size 보정 후 phylogroup별 C/G/E enrichment 분석 결과

---

## Results III. Deep-dive into an Outlier Strain (CP076082.1)

### 4.1 Anatomical Analysis of the Outlier: A Gap Within Group 1

Accessory genome 기반 계통 분석에서 **CP076082.1**은  
Group 1의 하위 클러스터에 속함에도 불구하고,  
동일 phylogroup 내 다른 균주들과 비교하여  
**비정상적으로 긴 가지(long branch)**를 형성하며 계통수의 외곽에 위치하였다.

ANI 분석에서는 Group 1과 높은 유전체 유사성을 유지하고 있음에도 불구하고,  
accessory gene presence/absence 패턴에서는  
뚜렷하게 이탈된 조성을 보였다.

---

### 4.2 Selective Functional Decay in Category E

C/G/E 카테고리 기반 Z-score 분석 및 PCA 결과,  
CP076082.1의 계통적 이탈은  
**아미노산 대사(E) 기능의 선택적 결손**에 의해 주도됨이 확인되었다.

> **Figure 5A.** Group 1 내 CP076082.1의 C/G/E boxplot  
> **Figure 5B.** C/G/E 기반 PCA plot

---

### 4.3 Genomic Evidence of Loss-and-Gain Strategy

CP076082.1은 아미노산 생합성 관련 경로(his-operon, methionine metabolism)를  
선택적으로 소실한 반면,  
특정 탄수화물 수송 및 이용 유전자(malK, ugpA, araQ)를 획득하였다.

이는  
**아미노산 자급 능력을 포기하는 대신  
환경 내 탄수화물 자원을 효율적으로 활용하는  
전략적 유전체 재구성**으로 해석된다.

---

### 4.4 Summary of Outlier Analysis

CP076082.1은  
Group 1의 기본적인 대사 골격을 유지하면서도,  
아미노산 대사 기능(E)을 선택적으로 축소한  
**기능적 변이체(functional variant)**로 정의된다.

이러한 기능 재편성은  
accessory genome matrix에 반영되어  
계통 트리 상에서 독립적인 장가지(long-branch) 분기를 형성하였다.
