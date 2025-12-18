# Lactobacillus fermentum Pan-Genome Analysis Report  
### **Defining L. fermentum Phylogroups Through Accessory Gene Clustering and Functional COG Profiling**

---

## 📌 Background  
*Lactobacillus fermentum*은 다양한 발효환경에서 분리되며, 탄수화물 대사·에너지 전환·아미노산 대사 등 발효 관련 기능에서 균주 간 큰 변이를 보인다.  
이 변이의 상당 부분은 **accessory genome 조합 차이**에서 기인하므로,  
발효 핵심 기능에 해당하는 COG(C, E, G)를 중심으로 분석하면 **기능적으로 의미 있는 phylogroup 정의**가 가능하다.

---

## 🎯 Objective  
본 분석에서는 다음을 목표로 한다:

1. **155 → ANI filtering → 81개**의 high-quality genome 기반 pan-genome 구축  
2. Accessory gene 기반 phylogroup 정의  
3. 각 그룹 간 COG(C/E/G) 기능 차이를 평가  
4. Outlier strain 분석을 통해 특정 COG의 결손이 phylogroup 패턴을 설명하는지 검증  

➡️ 최종적으로 *L. fermentum* phylogroup이 **단순 유사도 기반 그룹이 아니라, 발효 기능적 특성을 반영한 그룹**임을 규명한다.

---

# 1. 🧪 ANI-based Genome Filtering

## 1.1 Initial Dataset  
- NCBI에서 **총 155개 genome** 다운로드  
- 이후 모든 분석의 품질을 보장하기 위해 ANI 기반 genome filtering 수행  


## 1.2 Pairwise ANI Processing and ANI Heatmap   
<img width="3600" height="3000" alt="ani_heatmap" src="https://github.com/user-attachments/assets/12dec12e-fa07-4ee5-b7c0-6763e4dad9fe" />

- 155개 전체 Heatmap으로 계통적 패턴 및 outlier 시각적 확인  
- Pairwise ANI/AAI 계산 → 중복 제거  
- Outlier(AP017974.1) 탐지  
- NCBI “Complete” 수준만 유지 → **최종 81개 strain** 확정   

---

# 2. 🧬 Pan-Genome Construction

### 🔧 **Tools Used**
- **Prokka** (https://github.com/tseemann/prokka)  
  → *표준화된 genome annotation을 위해 사용*  
- **Roary** (https://github.com/sanger-pathogens/Roary)  
  → *core–accessory 구조 및 presence/absence matrix 생성용*  
- **EggNOG-mapper** (https://github.com/eggnogdb/eggnog-mapper)  
  → *각 유전자의 COG 기능 주석을 얻기 위해 사용*

---

## 2.1 Genome Annotation  
81개 genome을 **Prokka**로 균일하게 annotation하여 Roary 입력용 GFF 생성.

## 2.2 Roary Pan-Genome  
- Core / Accessory / Singleton 구조 도출  
- 분석 중 CP033371.1이 이상치로 확인되어 제외  
- 최종 gene_presence_absence matrix 확보  

## 2.3 Accessory Gene Set 추출  
- Core cutoff = 95%  
- 95% 미만 빈도 유전자를 모두 accessory로 정의(shell + cloud)  
➡️ 이 accessory genome이 **phylogroup clustering의 주요 feature**가 된다.  

---

# 3. 🔬 Functional Categorization of Accessory Genes

Accessory gene에 대해 EggNOG 기반 COG 주석을 수행하고,  
이 중 발효 핵심 카테고리인 **C(에너지), G(탄수화물), E(아미노산/탄수화물 조절)**에 집중하였다.

> EggNOG의 KO/EC 정보는 기능 비교에 직접적이지 않지만  
> **COG은 발효 미생물의 기능 분류가 명확하여 phylogroup 비교에 적합**하다.

---

## 3.1 Workflow  
1. EggNOG 결과에서 COG category 추출  
2. Accessory presence/absence matrix와 merge  
3. Strain별 C/G/E category count 생성  

➡️ 이렇게 얻은 **기능 signature(C/G/E)**는 이후 그룹 간 발효 대사 특성 비교의 핵심 자료가 된다.

---

# 🌐 4. Accessory Genome–Based Phylogrouping and Functional Interpretation

##  4.1 Pre-processing: Removal of CP033371.1 Prior to H35 Phylogrouping

<img width="1920" height="912" alt="image" src="https://github.com/user-attachments/assets/3cd9a848-692d-4cfd-93f9-8cf1a2c3e972" />

Accessory genome 기반 phylogrouping을 수행하기 전에,  
**CP033371.1이 트리에서 지나치게 긴 단독 branch를 형성하는 extreme outlier**임을 확인하였다.

- ANI 분석에서는 이상 없음  
- Roary 단계에서 gene clustering 혹은 annotation mismatch 가능성  
- 이로 인해 **거리 기반 분석(H35 분할 포함)을 왜곡할 정도의 비정상적 값** 생성

따라서 downstream 분석의 정확성을 위해  
**CP033371.1을 제외하고 총 80 strain 기준으로 H35 phylogrouping 및 후속 C/G/E 분석을 수행하였다.**

---

##  4.2 Phylogenetic Tree

Accessory 유전자 조성은 변동성이 크므로,  
**H35 threshold가 가장 안정적으로 phylogroup을 분해하는 cutoff**로 판단되었다.

Accessory gene presence/absence 기반 phylogenetic tree를 구성하고,  
strain별 C/G/E functional count를 매핑하여 **계통 구조–기능 연관성**을 평가하였다.

<img width="1183" height="828" alt="image" src="https://github.com/user-attachments/assets/234d12b6-65d1-47e3-a589-7c7ecc681bea" />


### 🔎 주요 관찰
1. **CP076082.1이 GROUP1 내부에서도 외곽 branch에 위치하는 outlier로 확인됨**  
2. Tree 구역별 **C/G/E functional profile이 상이한 패턴을 보임**

이 현상이 계통학적 요인인지 실제 기능 조성 차이인지 확인하기 위해 후속 분석을 수행하였다.

---

# 🔬 5. Functional Divergence Analysis of CP076082.1

##  5.1 C/G/E Functional Profile: Selective Reduction in Category E

Strain-level C/G/E count 비교 결과, CP076082.1은 GROUP1 대비 다음과 같은 특징을 보였다.

<img width="3000" height="1800" alt="CP0760821_CGE_boxplot_zscore" src="https://github.com/user-attachments/assets/1946ddf9-9e7f-4455-86d0-dd30c4215e11" />


- **C:** GROUP1과 유사  
- **G:** 경미한 감소  
- **E:** GROUP1에서 유일하게 **크게 감소 ⬇️**

→ 즉, **CP076082.1은 E 기능만 선택적으로 약화된 기능적 이탈형**이며, 이는 tree 상의 분리와 일치한다.

---

## 🧭 5.2 PCA Analysis: E Category Drives the Separation

 Strain별 C/G/E category count 파일로 변수(C/G/E)를 하나의 좌표축으로 축약해서, 서로 다른 strain들이 기능적으로 얼마나 비슷하거나 다른지를 시각화하는 PCA 분석 후 plot으로 나타내었다.

<img width="2400" height="1800" alt="GROUP1_only_CGE_PCA_plot" src="https://github.com/user-attachments/assets/71875945-e470-49ef-a1cd-72e3392d0920" />

C/G/E matrix 기반 PCA 결과:

- CP076082.1은 GROUP1 중심 cluster에서 벗어난 외곽 위치  
- **PC1·PC2 loading 분석(GROUP1_only_CGE_PCA_scores 파일)으로 봤을 때, E category가 분리를 결정하는 핵심 요소**

→ **E 기능 감소가 분리의 직접 원인**임을 통계적으로 확인하였다.

---

## 🧬 5.3 Accessory Genome Composition  
### (1) Major Loss of E-Category Pathways
 또한, 해당 strain만 가지고 있는 유전자, 아닌 유전자가 어떤 것이고, E 카테고리와 관련이 있는지 확인하기 위해 accessory_presence_absence_renamed.tsv 파일에서 악세사리 중의 악세사리 유전자를 탐색했다.
 
GROUP1 strain의 80% 이상이 공통 보유하는 accessory gene 중, CP076082.1에서만 결손된 기능성 경로는 다음과 같다.

- his-operon (hisA/B/C/D/G/H/Z/K): 히스티딘 생합성 전체  
- 메티오닌 대사 경로 (metI, metB, metE)  
- 성장·번역 관련 요소 (lepA, tuf)

→ 단순한 count 감소가 아니라 **핵심 대사 경로 단위의 실제 기능 결손**임을 의미한다.

### (2) Acquisition of Unique Transporters

CP076082.1 전용 unique accessory gene에서는 특정 carbon source 활용을 강화하는 transporter들이 확인되었다.

- malK — maltose transporter  
- ugpA — G3P transporter  
- araQ — arabinose transporter/regulator  

→ **E 기능 약화와 동시에 niche-specialized carbon metabolism 재구성**이 일어났음을 시사한다.

---

# 📊 6. Statistical Evaluation of C/G/E Functional Differences Between Groups

## 📌 6.1 Chi-square Test: Group-level Total Counts

<img width="40%" alt="cog_group_counts_h35_heatmap" src="https://github.com/user-attachments/assets/d024d185-b164-4468-a3b6-d5caefd101eb" />


phylogroup 간 **C/G/E total gene count** 동일성 검정:

- **C:** p ≈ 1.2e⁻²⁷⁶  
- **G:** p ≈ 5.1e⁻²¹⁷  
- **E:** 매우 작음  

→ 모든 기능군에서 **그룹 간 차이가 강하게 유의**하였다.

---

## 📌 6.2 Kruskal–Wallis Test: Strain-level Distribution

<img width="2400" height="1200" alt="cog_CGE_h35_boxplot_smallfliers" src="https://github.com/user-attachments/assets/f79f3baf-2cf0-427e-bf90-14eeab05ea80" />

strain-level C/G/E 분포 비교 결과:

- **C:** p ≈ 6.4e⁻⁶  
- **G:** p ≈ 3.6e⁻⁵  
- **E:** p ≈ 3.7e⁻¹¹  

→ strain-level에서도 **그룹 간 차이가 유지**됨을 확인하였다.

---

## 📌 6.3 Fisher’s Exact Test: Genome-size–Adjusted Ratio Differences

<img width="1560" height="960" alt="fisher_CGE_enrichment_h35_barplot" src="https://github.com/user-attachments/assets/daf6f4c0-8f05-49d4-9899-ae4fc93d1a64" />


Accessory genome size 차이를 보정한 비율 검정 결과:

- **Group 1:** OR = 1.297, p = 4.9e⁻20 → **enriched ⬆️**  
- **Group 2, 4, 5:** depleted  
- **Group 3:** neutral  

→ 기능 차이는 genome size artifact가 아닌 **구조적 특성**임을 확인하였다.

---

# 🧩 7. Integrated Conclusion

본 연구의 두 축 분석 —  
① **H35 phylogroup 간 C/G/E 기능 조성 차이**,  
② **동일 그룹 내 outlier(CP076082.1)의 기능·구조 이탈 분석** — 을 종합하면 다음과 같다.

## ✨ 핵심 결론

1. **C/G/E 기능군은 L. fermentum의 발효 대사를 규정하는 핵심 functional axis이며,**  
   phylogroup(H35) 간 이 축에서 통계적으로 유의한 차이가 존재한다.  
   → 즉, accessory genome 기반 phylogroup은 단순 gene 조성 차이가 아니라  
     **각 그룹의 고유한 발효 대사 signature**를 반영한다.

2. **CP076082.1은 E category 핵심 경로 결손 + 특이 carbon transporter 획득이라는 재구성을 통해**  
   동일 그룹 내에서도 기능적으로 이탈된 패턴을 보였다.  
   → 즉, **군집 분리와 outlier 형성 모두가 C/G/E 기능축의 구조적 변화와 직접적으로 연결된다.**

---

# ⭐ Final Statement

**“Accessory genome 기반 phylogroup은 L. fermentum의 발효 대사 signature를 명확히 구분하는 지표이며,  
CP076082.1과 같은 outlier의 분리 또한 C/G/E 기능축 재구성으로 설명된다.”**
