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

1. **145 → ANI filtering → 81개**의 high-quality genome 기반 pan-genome 구축  
2. Accessory gene 기반 phylogroup 정의  
3. 각 그룹 간 COG(C/E/G) 기능 차이를 평가  
4. Outlier strain 분석을 통해 특정 COG의 결손이 phylogroup 패턴을 설명하는지 검증  

➡️ 최종적으로 *L. fermentum* phylogroup이 **단순 유사도 기반 그룹이 아니라, 발효 기능적 특성을 반영한 그룹**임을 규명한다.

---

# 1. 🧪 ANI-based Genome Filtering

## 1.1 Initial Dataset  
- NCBI에서 **총 155개 genome** 다운로드  
- 이후 모든 분석의 품질을 보장하기 위해 ANI 기반 genome filtering 수행  

## 1.2 Pairwise ANI Processing  
- Pairwise ANI/AAI 계산 → 중복 제거  
- Outlier(AP017974.1) 탐지  
- NCBI “Complete” 수준만 유지 → **최종 81개 strain** 확정  

## 1.3 ANI Heatmap (📊 위치 제안: Filtering 근거 제공)  
- 155개 전체 Heatmap으로 계통적 패턴 및 outlier 시각적 확인  
- Filtering 후 81개의 균주가 더 응집된 패턴을 형성함  

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

## 2.4 EggNOG Functional Annotation  
Accessory gene 리스트를 EggNOG-mapper 결과와 병합하여 **각 유전자의 COG category** 확보.

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

# 4. Accessory Genome–Based Phylogrouping and Functional Interpretation
## 4.1 Pre-processing: Removal of CP033371.1 Prior to H35 Phylogrouping

Accessory genome 기반 phylogrouping을 진행하기에 앞서,  
**CP033371.1은 트리에서 지나치게 긴 단독 branch를 형성하는 extreme outlier**로 나타났다.

- ANI 분석에서는 이상이 없었으나,  
- Roary 단계에서 gene clustering 혹은 annotation mismatch가 발생한 것으로 추정되며,  
- 이로 인해 **전체 그룹 내 통계 분석(H35 분할 포함)을 왜곡할 수준의 비정상적 거리값**을 보였다.

따라서 downstream 분석의 정확성을 유지하기 위해  
**CP033371.1은 사전에 제외하고, 총 80개의 strain을 기준으로 H35 phylogrouping 및 후속 C/G/E 통계 분석을 수행하였다.**

---

## 4.2 phylogenetic tree
Accessory 유전자 조성은 strain 간 변동성이 높기 때문에, **H35 threshold는 가장 안정적으로 phylogroup이 분해되는 cutoff**로 판단되어 본 연구의 tree construction에 사용되었다.
Accessory gene presence/absence 정보를 기반으로 phylogenetic tree를 구성하였고, 여기에 strain별 C/G/E functional count를 매핑하여 계통 구조와 기능적 특성의 연관성을 평가하였다.

이 분석에서 두 가지 핵심 관찰이 나타났다:

1. **CP076082.1이 동일 phylogroup(GROUP1) 내에서 외곽 branch에 위치하는 outlier로 확인됨**
2. **Tree 상 각 구역에서 C/G/E functional profile이 서로 다른 패턴을 나타냄**

이러한 분리 현상이 단순 계통학적 요인인지, 실제 기능적 조성 차이에 기반한 것인지 확인하기 위해 후속 분석을 수행하였다.

---

# 5. Functional Divergence Analysis of CP076082.1

## 5.1 C/G/E Functional Profile: Selective Reduction in Category E
Strain-level C/G/E category count 비교 결과, CP076082.1은 GROUP1 평균 대비 다음과 같은 특징을 보였다.

- **C (Energy production & conversion):** GROUP1과 유사  
- **G (Carbohydrate transport & metabolism):** 경미한 감소  
- **E (Amino acid/carbohydrate-associated metabolism):** GROUP1에서 유일하게 **크게 감소**

즉, CP076082.1은 C/G 기반 발효 기능은 유지하면서도 **E category 기능만 선택적으로 약화**된 기능적 이탈형임이 확인되었다. 이는 phylogroup tree에서 중심 cluster로부터 멀어진 패턴과 일치한다.

---

## 5.2 PCA Analysis: E Category Drives the Separation
C/G/E matrix에 대한 PCA에서도 CP076082.1은 GROUP1 cluster 중심에서 벗어난 외곽 위치를 차지했다.  
PC1·PC2 loading 분석 결과 **E category 유전자군이 분리를 결정하는 가장 큰 요소**로 기여하였다.

따라서 CP076082.1이 phylogroup에서 분리된 직접적 요인은 E 기능 감소임이 통계적으로도 확인되었다.

---

## 5.3 Accessory Genome Composition: Loss of Core E Pathways + Acquisition of Unique Transporters

### (1) Major Loss of E-Category Pathways  
GROUP1의 80% 이상이 공통 보유하는 accessory gene 중, CP076082.1에서만 결손된 기능성 경로는 다음과 같다.

- his-operon (hisA/B/C/D/G/H/Z/K) — 히스티딘 생합성 전체 과정  
- 메티오닌 대사 경로 (metI, metB, metE)  
- 성장·번역 관련 요소 (lepA, tuf)  

→ 이는 E category 감소가 단순 count 변화가 아니라 **핵심 경로 수준의 실제 결손**임을 의미한다.

### (2) Acquisition of Unique Transporters  
반대로 CP076082.1에서만 등장하는 unique accessory gene 분석 결과, 특정 carbohydrate utilization을 강화하는 transporter가 확인되었다.

- malK — maltose transporter  
- ugpA — G3P transporter  
- araQ — arabinose transporter/regulator  

→ 이는 E category 약화와 동시에 **특정 carbon source 중심으로 niche-specialized 재구성**이 일어났음을 시사한다.

---

# 6. Statistical Evaluation of C/G/E Functional Differences Between Groups

앞서 phylogenetic tree 시각화에서 확인된 **(2) C/G/E functional pattern의 그룹 간 차이**가  
실제 통계적으로 유의한지 확인하기 위해 세 단계의 검정을 수행하였다.

---

## 6.1 Chi-square Test: Group-level Total Counts
각 phylogroup이 보유한 **C/G/E total gene count** 분포가 동일한지 검정하였다.

- **C:** p ≈ 1.2e⁻²⁷⁶  
- **G:** p ≈ 5.1e⁻²¹⁷  
- **E:** p ≈ 매우 작음  

→ 세 기능군 모두 그룹 간 total gene count 차이가 **강하게 유의**하였다.  
→ 즉, phylogroup 간 분리는 accessory genome 구조뿐 아니라 **기능적 조성 차이**도 반영한다.

---

## 6.2 Kruskal–Wallis Test: Strain-level Distribution
각 strain의 C/G/E gene count 분포를 그룹 간 비교하였다.

- **C:** p ≈ 6.4e⁻⁶  
- **G:** p ≈ 3.6e⁻⁵  
- **E:** p ≈ 3.7e⁻¹¹  

→ strain-level에서도 그룹 간 분포 차이가 존재하였다.  
→ 즉, functional differentiation은 **군집의 구조적 특성**이다.

---

## 6.3 Fisher’s Exact Test: Genome-size–Adjusted Ratio Differences
Group 간 accessory genome size가 서로 다르기 때문에,  
단순 유전자 수 차이가 아닌 **비율 차이**가 유지되는지 Fisher 검정을 수행하였다.

- **Group 1:** OR = 1.297, p = 4.9e⁻20 → C/G/E 비율이 유의하게 **enriched**  
- **Group 2, 4, 5:** depleted  
- **Group 3:** 중립  

→ Genome size를 보정해도 C/G/E 비율 차이는 유지되었다.  
→ 기능적 차이는 단순 genome size artifact가 아니라 **그룹 고유의 구조적 특성**임이 확인되었다.

---

## Interpretation of the Three Tests
세 단계의 통계 검정은 모두 일관된 결론을 지지한다.

- **그룹 간 C/G/E total 수는 유의하게 다르며**,  
- **strain-level에서도 분포가 다르고**,  
- **genome size 보정 후에도 C/G/E 비율 차이는 유지된다.**

이는 트리에서 관찰된 기능적 패턴(각 phylogroup의 distinct metabolic signature)이  
단순 우연이나 데이터 artifact가 아니라,  
**실제 biological differentiation**임을 강하게 뒷받침한다.

---

# 7. Integrated Conclusion

본 연구에서 얻은 두 축의 분석 결과—  
① **H35 phylogroup 간 C/G/E 기능군의 유의미한 차이**,  
② **동일 그룹 내 outlier(CP076082.1)의 기능·구조적 이탈 분석**—을 통합하면  
다음 결론을 도출할 수 있다.

1. **C/G/E 카테고리는 L. fermentum의 발효 대사를 규정하는 핵심 기능축**이며,  
   phylogroup(H35) 간 이 기능군에서 통계적으로 유의한 차이가 존재한다.  
   → 즉, accessory genome 기반 phylogroup은 단순한 유전자 조성 차이가 아니라,  
     **발효 대사 기능 조합(functional signature)의 차이를 반영한다.**

2. 동일 그룹 내 outlier(CP076082.1) 분석에서도  
   **E category 감소 및 고유 transporter 획득이라는 기능적 재구성**이  
   해당 strain이 phylogroup 중심에서 벗어난 원인임이 확인되었다.  
   → 즉, **군집 분리 + outlier 형성 모두가 C/G/E 기능 변화와 직접적으로 연결**된다.

---

## Final Statement
따라서 다음의 결론이 도출된다.

**“발효 관련 핵심 기능군(C/G/E)은 phylogroup을 구분하는 데 충분한 기능적 지표이며,  
각 그룹은 서로 다른 발효 대사 signature를 가진다.  
CP076082.1과 같은 outlier의 분리 또한 이러한 기능적 축의 재구성에서 기인한다.”**

