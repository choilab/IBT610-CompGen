
# 🧫 Lactobacillus fermentum 유전체 비교 분석

## 1. 연구 배경 (Background)
*Lactobacillus fermentum*은 그람양성, 이형발효(heterofermentative) 유산균으로,  
발효식품, 식물, 인체 장내 등 다양한 환경에서 발견됩니다.  
산성 및 담즙 환경에서도 생존력이 높아 **프로바이오틱 균주로서의 잠재력**이 큽니다.

 **총 394개의 공개된 유전체 서열**중에서 중복 데이터를 제거한 후 **155개의 완전(Complete) 유전체**를 대상으로 분석을 진행했습니다.
 * nr set != 100%

## 지난 주 요약
| **ANI ≥ 95% 및 AAI ≥ 95%** | 높은 유사도를 대부분 나타냈지만, 특정 id를 가진 sample에서는 80% 초반의 유사도를 가지는 값이 있었습니다.


따라서, 이번 과제를 통해 **동일 종으로 신뢰할 수 있는 클러스터만 유지**하고, 어떻게 그룹화 할 지 고민했습니다. (불일치 케이스는 별도로 추가 분석할 수도 있음)

---

## 2. 데이터 전처리 과정 (Data Processing Pipeline)

### 1️⃣ 데이터 구조 통일
- ANI 및 AAI 결과를 `GenomeA`, `GenomeB`, `Value` 형식의 **pairwise 구조**로 변환  
- 대칭쌍(A–B, B–A) 및 중복값 제거

### 2️⃣ 공통쌍 추출 및 병합
- `pd.merge()`를 이용하여 **inner join** 수행  
- 공통된 조합(`GenomeA`, `GenomeB`)만 유지  
- 전체 155개 유전체의 가능한 조합: **11,935쌍**  
  - ANI 결과: 11,935쌍  
  - AAI 결과: 11,781쌍 (일부 유전체 누락-하나의 fasta파일이 aai분석에서 문제가 되는 파일이었기에 제거함)

### 3️⃣ 불일치 쌍 탐지 (abs_diff, flag)
- `abs_diff` = ANI와 AAI의 절대 차이  
- `flag` = 기준 불일치 시 표기  
  - 예: ANI ≥ 95% & AAI < 95%  
  - 또는 반대 경우  
→ 불일치 쌍만 따로 추출하여 **low consistency** 그룹으로 분류

### 4️⃣ Low-ANI 유전체 탐지
- **AP017974.1** 유전체가 낮은 ANI 쌍에서 **308회 이상 반복 등장**,  
  이는 특이적 변이 혹은 계통적으로 다른 분리주임을 시사. 다른 유전체들도 2회 정도 low-ANI 쌍에서 나타났는데, 이는 해당 id와 비교할 때 낮은 값을 가지게 된 것이기에 AP017974.1가 문제가 됨.
> 파일:https://github.com/igchoi/IBT610-CompGen/blob/c7f2158e09483995fcbe89c653a3632d4fc43d5f/2025-Fall/mjbaek/low_ani_outlier_pairs.csv
---

## 3. 메타데이터 요약 (Metadata Summary)

| 항목 | 요약 |
|------|------|
| **종(Species)** | 모두 *Limosilactobacillus fermentum* |
| **Assembly Level** | 모두 Complete Genome 수준 |
| **Assembly Name** | 총 82개의 고유한 어셈블리 이름(총 155개 중 MLST를 분석할 수 있는 ACCESSION number를 가진것은 81개였음) |
| **BioSample / BioProject** | 다수 결측 (빈값 많음) |
| **Release Date** | 대부분 미기재 |
| **고유 BioSample ID 수** | 82개 (중복 제거 후) |

> ⚠️ 중복된 BioSampleAccession은 동일 시료가 여러 프로젝트나 시점에서 재등록되었음을 의미합니다.
> 파일: https://github.com/igchoi/IBT610-CompGen/blob/e72b45430ba6900e04f5af59094fcd7c1de298dd/2025-Fall/mjbaek/biosample_metadata_155.csv
---

## 5. MLST (Multilocus Sequence Typing) 분석

보통의 MLST분석은 다음과 같습니다.

| 구분         | 이름                                                      | 특징                                                                                                            | URL                                                                                    |
| ---------- | ------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------- |
| **웹 기반**   | **PubMLST**                                             | 가장 표준적인 플랫폼. 각 종(species)별로 curated된 MLST scheme과 allele database를 보유. 자동으로 ST 지정 가능.                         | [https://pubmlst.org/](https://pubmlst.org/)                                           |
|            | **EnteroBase**                                          | Enterobacteriaceae, *Listeria*, *Salmonella* 등 대규모 comparative genomics에 적합. 자동 QC 및 phylogenetic tree 기능 포함. | [https://enterobase.warwick.ac.uk/](https://enterobase.warwick.ac.uk/)                 |
|            | **BIGSdb (Bacterial Isolate Genome Sequence Database)** | MLST, cgMLST, wgMLST 등 고해상도 타이핑 가능. PubMLST와 연동.                                                              | [https://bigsdb.pasteur.fr/](https://bigsdb.pasteur.fr/)                               |
| **로컬 실행형** | **mlst (by Torsten Seemann)**                           | NCBI 등에서 genome FASTA를 받아서 CLI로 빠르게 ST 자동 지정. (`brew install mlst` 혹은 conda 설치 가능)                            | [https://github.com/tseemann/mlst](https://github.com/tseemann/mlst)                   |
|            | **ARIBA, SRST2**                                        | raw read 기반으로 allele 매칭 수행 (assembly 불필요). NGS pipeline과 연동 가능.                                               | [https://github.com/sanger-pathogens/ariba](https://github.com/sanger-pathogens/ariba) |
|            | **MLSTcheck**                                           | genome assembly 파일로 MLST scheme 자동 인식 및 ST 부여. python 기반.                                                     | [https://github.com/tseemann/mlst_check](https://github.com/tseemann/mlst_check)       |

웹기반의 PubMLST, 로컬 실행형의 mlst를 사용했으나, Lactobacillus fermentum에 대한 tool은 제공하지 않고 있었습니다.
따라서, tool을 이용한 분석이 불가능해 (참고: Lactobacillus fermentum의 경우 PubMLST에 Lactobacillus 속 일부 종(예: L. plantarum, L. casei)은 MLST scheme이 구축되어 있지만,
L. fermentum은 종종 표준화된 MLST scheme이 부족하거나 **연구자 맞춤형(custom MLST)**를 사용합니다.-from gpt)gpt의 말을 참고하여 논문을 통해 house keeping gene들을 찾았습니다.
> 참고 논문: https://pmc.ncbi.nlm.nih.gov/articles/PMC4437502/

논문을 통해 11개의 housekeeping gene 서열을 얻었고, 이를 기반으로 mlst분석을 시행했습니다.

### 🧬 분석 대상 유전자
MLST에 사용된 *housekeeping gene* 11개 : clpX, dnaA, dnaK, groEL, murC, murE, pepX, pyrG, recA, rpoB, uvrC

### 🧠 분석 목적
- 각 유전체 내 **핵심 보존 유전자 서열의 변이(돌연변이)**를 파악  
- 이를 통해 **계통적 다양성과 집단 구조**를 밝힘

### 🧪 분석 절차
1. 가지고있는 각 유전체(FASTA)에서 위의 11개 유전자 서열 검색  
2. 동일 서열에 같은 **allele 번호** 부여  
3. 유전자 조합에 따라 고유한 **Sequence Type (ST)** 정의  
4. 최종적으로 GenomeID × 11개 유전자 매트릭스 생성 및 CSV 저장
> 파일: https://github.com/igchoi/IBT610-CompGen/blob/a39f93c2de325c4220dec52770bd688e9e527986/2025-Fall/mjbaek/analysis_result.csv

### 💾 결과 요약
- 총 추출된 유전자 서열: **426개**  
  (이론상 704개 = 64개 유전체 × 11유전자보다 적음)
  
- 일부 유전체에서 유전자를 찾지 못한 이유:
  1. **NCBI 데이터 갱신/폐기**로 인한 누락  
  2. **어노테이션 불일치** (유전자명 표기 차이, 미기재 등)

---

## 6. ANI/AAI – MLST-Metadata 통합 분석 결과

### 📌 주요 관찰점 및 후속 연구 방향
- Low ANI 쌍에 반복적으로 등장한 **AP017974.1 (BioSample: SAMD00073748)**  → 제거하고 전체 분석을 앞으로 진행하나, 동일한 st가 있어 어떤 차이가 있는지 추후 분석할 수 있을듯 함. 
- 전체적으로 **ST 다양성이 매우 높게 나타남**, 즉 *L. fermentum* 내에서도 유전적 이질성이 큼.
- **ID는 다르지만 ST가 동일한 strain**이 일부 관찰됨 → 이는 **서로 다른 분리주이지만 유전자형(MLST 유전자 조합)이 유사**함을 의미함.->ST는 같지만 ID가 다른 strain 간의 genome alignment | SNP 및 유전자 gain/loss 확인 
- **MLST 유형과 메타데이터(분리원: 인체/식품/식물 등) 연계 분석**
 -> L.fermentum이 발효에 관여하다보니, metadata 정보에서 milk, kimchi 유래에서 분리된 정보가 많았다. 따라서, 발효식품별로 서열을 분리해보고 공통적인 특징이 있거나 서열에서 가지는 돌연변이 서열이 같은지 mlst를 참고해 프로젝트의 목표처럼 그룹화 할 수 있을 듯하다.

---


