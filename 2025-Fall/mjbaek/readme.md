## 📌 Introduction

*Lactobacillus fermentum*은 유제품, 채소 발효식품, 인체 장·구강 등 다양한 발효 환경에서 반복적으로 분리되며, 대체로 **1.8–2.3 Mb** 크기의 유전체를 가진다(Azmal Ali et al., 2018).  
동일 종 내에서도 균주 간 탄수화물 이용성, 에너지 전환, 아미노산 대사 및 스트레스 내성과 같은 **발효 관련 phenotype과 대응 유전자 조성에서 큰 변이**가 보고되어 왔다(Konyanee et al., 2019).

이러한 기능적 변이의 상당 부분은 보존적인 **core genome**만으로는 설명하기 어렵고, 환경 적응과 niche 특이성을 반영하는 **accessory genome 조합 차이**에서 기인한다는 점이 지적되어 왔다(Tatusov et al., 2001).
그럼에도 불구하고 기존 *L. fermentum* 유전체 연구는 주로 ANI 또는 core genome 기반 계통 분석에 집중되어 있으며,accessory genome 패턴과 발효 관련 기능을 통합하여 phylogroup을 정의한 체계적 분석은 제한적이다.

따라서 본 연구에서는 공개된 *L. fermentum* 유전체 **155개**를 ANI 기준으로 필터링하여 **80개의 high-quality genome**을 선별한 뒤 pan-genome을 구축하고, accessory gene의 존재/부재 패턴을 기반으로 phylogroup을 정의하였다.  
이어 각 균주의 accessory 유전자에 **COG 기반 기능 주석**을 부여하고, 발효와 직접적으로 연관된 **C (Energy), E (Amino acid), G (Carbohydrate)** 카테고리의 조성을 정량화하여 phylogroup 간 발효 대사 시그니처를 비교하였다.

또한 동일 phylogroup 내에서 관찰되는 **outlier strain**의 경우, 발효와 관련된 특정 COG 기능 유전자의 소실 또는 획득이해당 균주가 다른 위치에 분리되는 현상을 설명할 수 있는지를 분석함으로써,
L. fermentum phylogroup이 단순한 유전체 유사도가 아니라 발효 기능 차이를 반영한 기능적 그룹임을 확인하고자 하였다.


# 🧬Materials and Methods

## 1.1 Genome Sequencing Data Acquisition

NCBI RefSeq 데이터베이스에서 *Lactobacillus fermentum*으로 분류된  
모든 유전체(총 **155개**)를 FASTA 형식으로 다운로드하였다.  
Genome assembly level이 “Complete Genome”인 시퀀스만 1차 선별하였으며,  
이후 fastANI를 이용해 균주 간 쌍별 ANI를 계산하였다.

ANI 95% 미만 값을 보이는 종 내 이질적 균주  
(예: AP017974.1 등)는 분석에서 제외하였고,  
중복 및 저품질 유전체를 제거한 후  
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

## 1.2 Genome Annotation and Pan-Genome Construction

각 균주의 주석 형식을 통일하기 위해  
Prokka (v1.14.6)를 사용하여 ORF를 재예측하고 기능을 주석화하였다.  
각 균주별 GFF 파일과 단백질 서열(FAA)을 생성하였다.

생성된 Prokka GFF 전체를 Roary (v3.13.0)에 입력하여  
pan-genome을 구축하였으며,  
gene_presence_absence.* 및 gene_presence_absence.Rtab 파일을 통해  
core 및 accessory 유전자 집합을 정의하였다.  
전체 균주의 **95% 이상에서 발견되는 유전자**를 core genome으로,  
그 외(shell + cloud)를 accessory genome으로 정의하였다.

Accessory genome 기반 계통수에서  
비정상적으로 긴 branch를 형성한 **CP033371.1**을 추가로 제거하였고,  
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

## 1.3 Functional Categorization (COG Profiling)

Roary로부터 얻은 pan-genome 단백질 서열 집합  
(또는 Prokka FAA 전체 병합본)을 eggNOG-mapper v2.1.9에 입력하여  
eggNOG 5.0 데이터베이스 기반의 기능 주석과 COG 카테고리를 할당하였다.

이 중 에너지 대사 [C],  
탄수화물 대사 및 수송 [G],  
아미노산 대사 및 수송 [E]  
세 카테고리를 선택하여  
strain × C/G/E 카운트 및 비율 매트릭스를 구성하였다.

CP033371.1을 제외한 **80개 균주**에 대해  
`cog_by_strain_CGE_h35_no_CP033371.tsv` 파일을 생성하여  
phylogroup별 C/G/E 조성 비교 및 시각화에 사용하였다.

**Workflow**  
- **Input**:  
  - [Pan-genome FAA](data/prokka/)
- **Tool**:  
  - [eggNOG-mapper v2.1.9](https://github.com/eggnogdb/eggnog-mapper)
- **Output**:  
  - [cog_by_strain_CGE_h35_no_CP033371.tsv](data/cog/cog_by_strain_CGE_h35_no_CP033371.tsv)

---

## 1.4 Phylogrouping and Statistical Analysis

Accessory genome 기반 phylogroup 정의를 위해  
Roary에서 생성된 accessory gene 계통수  
(`accessory_tree.nwk`)를 사용하였다.

계통수 기반 거리 정보를 이용해 계층적 군집화를 수행하고,  
**cut height = 0.35** 기준으로 phylogroup을 정의하였다.  
최종 phylogroup 정보는  
[`phylogroup_h35.tsv`](data/metadata/phylogroup_h35.tsv)로 정리하였다.

Phylogroup 간 기능 차이를 평가하기 위해  
각 strain의 C/G/E 유전자 수 및 비율을 phylogroup 정보와 결합한 뒤,  
카이제곱 검정과 Kruskal–Wallis 검정을 통해  
phylogroup 간 분포 차이의 유의성을 평가하였다.  
또한 유전체 크기 및 총 유전자 수 차이에 따른 편향을 보정하기 위해  
Fisher의 정확 검정을 이용한 enrichment 분석을 추가로 수행하였다.

마지막으로 strain × C/G/E 비율 매트릭스를 입력으로  
PCA 분석을 수행하여 주요 성분(PC1, PC2) 상에서 strain을 투영하고,  
phylogroup별 색상 구분을 통해  
기능적 프로파일의 군집 패턴을 시각화하였다.

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

