# Integrated Pangenome & Surfaceome Analysis of *Lactiplantibacillus plantarum*

> **Project Duration:** 2025.10.16 ~ 2025.11.27 (Current)  
> **Target:** *Lactiplantibacillus plantarum* Representative Genomes  

## 📖 Project Overview
본 프로젝트는 NCBI 데이터베이스에 등록된 방대한 *L. plantarum* 유전체 데이터를 정제하여 **대표 유전체(Representative Genomes)**를 선별하고, 이들의 **Pangenome(범유전체)** 구조와 숙주 상호작용에 핵심적인 **Surfaceome(표면 단백질체)**을 규명하는 것을 목표로 합니다.

## 📅 Project Log & Timeline (Daily Progress)

### 1. Data Acquisition & Initial Statistics (2025.10.16)
**"Lactobacillus 유전체 현황 파악 및 데이터 수집"**
NCBI Assembly Database에서 *Lactobacillus* 속 전체에 대한 통계를 확보하고, 분석 대상인 *L. plantarum*의 전체 유전체를 다운로드했습니다.

| Category | Count | Description |
| :--- | :--- | :--- |
| Total assemblies | 5,427 | *Lactobacillus* 전체 (Draft + Complete) |
| **Target Download** | **2,976** | ***L. plantarum* 종 전체 데이터 확보** |
| RefSeq Annotated | 2,326 | PGAP 주석 완료된 고품질 데이터 |

---

### 2. Analysis Design & Method Setup (2025.10.30)
**"분석 방법론 정립"**
- **ANI vs AAI:** 유전체 유사도 측정을 위한 FastANI 및 AAI(Amino Acid Identity) 알고리즘 비교 및 선정.
- **Model Genome:** 참조 균주(Reference strain) 선정 기준 확립.
- **Genetic Variation:** SNP, InDel, MLST 등 미세 변이 분석 파이프라인 구상.

---

### 3. Genome Selection & QC (2025.11.04)
**"고품질 대표 균주 선별 (Genome Reduction)"**
초기 2,976개 중 Chromosome-level 359개를 1차 선별 후, 중복을 제거하여 최종 분석 세트를 확정했습니다.

- **Quality Control:** CheckM2 (Completeness > 98%, Contamination < 1%) 및 N50 값 기반 필터링.
- **Clustering (FastANI):** 359개 균주에 대해 All-vs-All ANI 계산 수행.
- **Threshold Decision:**
    - ANI 99.9% 기준 적용 시: **214개 그룹** 형성 (최종 선택).
    - *참고: SKANI 100% 기준 적용 시 195개 그룹 형성.*

https://github.com/igchoi/IBT610-CompGen/blob/main/2025-Fall/msan/result/fastani_heatmap_99.5-100.0.png

### 4. Pipeline Development (2025.11.13)
**"AAI 및 cgMLST 분석 환경 구축"**
선별된 214개 균주의 정밀 분석을 위해 Conda 기반의 독립된 분석 환경을 구축하고 스크립트를 작성했습니다.

- **Environment 1 (`aai_env`):** `CompareM`을 이용한 AAI 분석 및 매트릭스 시각화.
- **Environment 2 (`chewbbaca_env`):** `chewBBACA`를 이용한 cgMLST(core genome MLST) 스키마 생성 및 Allele Calling.

---

### 5. Pangenome Profiling (2025.11.20)
**"Roary를 이용한 유전자 다양성 분석"**
214개 대표 균주에 대해 Prokka Annotation 후 Roary를 수행하여 유전자 풀(Gene pool)을 분류했습니다.

| 유전자 그룹 | 기준 (Strains %) | 유전자 수 | 의미 |
| :--- | :--- | :--- | :--- |
| **Core** | 99% ≤ | **1,173** | 종(Species)의 핵심 기능 유지 |
| **Shell** | 15% ≤ < 95% | 1,907 | 환경 적응 및 변이 |
| **Cloud** | < 15% | **17,274** | 균주 특이적 희귀 유전자 (다양성의 원천) |
| **Total** | - | **20,970** | 전체 유전자 군집 (Cluster) |

---

### 6. Functional Integration & Visualization (2025.11.27 - Current)
**"Surfaceome 예측 및 최종 데이터 시각화"**
`lp_final.py` 파이프라인을 통해 Pangenome 데이터와 기능 분석을 통합하고 시각화했습니다.

#### A. Pangenome Expansion (Heaps' Law)
유전체 수가 늘어날수록 신규 유전자가 계속 발견되는 **Open Pangenome** ($\lambda = 0.378$) 특성을 확인했습니다.
![[Heaps Law](2025-Fall/msan/result/02_heaps_law.png)](https://github.com/igchoi/IBT610-CompGen/blob/main/2025-Fall/msan/result/02_heaps_law.png)

#### B. Phylogenetic Clustering (Mash)
Mash distance 기반의 계통수를 통해 214개 균주 간의 유전적 거리를 시각화했습니다.
![[Mash Clustering](images/03_mash_clustering.png)](https://github.com/igchoi/IBT610-CompGen/blob/main/2025-Fall/msan/result/03_mash_clustering.png)

#### C. Surfaceome Prediction (SignalP + TMHMM으로 재분석을 사용하여 업데이트 예정)(데이터 삭제 예정)
SignalP 및 TMHMM 로직을 적용하여 세포 표면 단백질(Secreted, Membrane, Cell wall)을 예측하고, CAR(Core/Accessory/Rare) 카테고리별 분포를 분석했습니다.
![[Surfaceome](images/05_surfaceome.png)](https://github.com/igchoi/IBT610-CompGen/blob/main/2025-Fall/msan/result/05_surfaceome.png)

---

cat << 'EOF' > ~/quality_control_analysis_kr.md
# Pangenome 품질 관리 및 Surfaceome 분석

## 📋 개요

본 문서는 세균 유전체 데이터에 대한 품질 관리 파이프라인과 Surfaceome/Secretome 분석 과정을 설명합니다.

---

## 1. 데이터 품질 관리 (Data Quality Control)

유전체 데이터의 기술적 오류를 배제하고 생물학적 유의성을 높이기 위해 추가적인 **이상치 탐지(Outlier Detection)** 작업을 수행했습니다.

### 🔹 1.1 FASTA 품질 검사 및 필터링

| 단계 | 설명 |
|------|------|
| **문제 식별** | 초기 데이터셋에 전체 시퀀스가 병합된 파일(`combined_all`) 및 품질 저하 샘플 포함 확인 |
| **기준 적용** | 서열 수 및 길이 분포 분석 후 평균에서 **±2 표준편차(SD)** 벗어나는 샘플 식별 및 제거 |

#### 분석된 품질 지표:
- 샘플당 서열 수
- 평균 서열 길이
- 길이 분포 패턴

### 🔹 1.2 Roary 이상치 제거

**검증**: Gene Presence/Absence 매트릭스에서 유전자 수가 비정상적으로 적거나 많은 **샘플 9개 제거**.

#### 결과 요약 (198 → 189 균주):

| 지표 | 제거 전 | 제거 후 | 변화 |
|------|---------|---------|------|
| 총 샘플 수 | 198 | 189 | -9 |
| Core Genes (100%) | 1,173 | 1,223 | **+50** |
| Pangenome 견고성 | - | - | ✅ 향상 |

#### 주요 발견:
- ✅ 이상치 제거 후 Core gene 수가 **50개 증가**
- ✅ 클러스터링 히트맵 분석 결과, 노이즈 제거로 **균주 간 패턴이 명확해짐**
- ✅ **Pangenome 견고성(Robustness) 향상** 확인

![Fasta Group Analysis](https://raw.githubusercontent.com/igchoi/IBT610-CompGen/main/2025-Fall/msan/result/06_fasta_group_analysis.png)

![Fasta Group Heatmap](https://raw.githubusercontent.com/igchoi/IBT610-CompGen/main/2025-Fall/msan/result/07_fasta_group_heatmap.png)

---

## 2. Surfaceome 및 Secretome 분석

숙주 상호작용의 핵심 인자인 **'세포 밖으로 분비되는 단백질(Secretome)'**을 선별하기 위해 구조 예측 파이프라인을 구축했습니다.

### 🔧 사용 도구

| 도구 | 버전 | 목적 |
|------|------|------|
| **SignalP** | 6.0 | 신호 펩타이드 예측 |
| **TMHMM** | 2.0 | 막관통 헬릭스 예측 |


# Secretome 선별 의사 코드
분비_단백질 = 단백질.필터(
    (SignalP_예측 != "OTHER") &  # 신호 펩타이드 보유
    (TMHMM_헬릭스 == 0)          # 막관통 도메인 없음
)

🔹 2.3 최종 결과
출력물	설명
📂 Final_Secreted_Proteins.xlsx	선별된 분비 단백질 목록
📊 Secretome_Heatmap.png	균주별 분비 단백질 분포 히트맵
통계:
범주	수량	비율
분석된 전체 CDS	~655,000	100%
신호 펩타이드 양성	32,653	~5%
최종 분비 단백질 (SP+ & TM-)	TBD	-

![Type Distribution Bar](https://raw.githubusercontent.com/igchoi/IBT610-CompGen/main/2025-Fall/msan/result/02_Type_Distribution_Bar.png)
![Signal Peptide Presence](https://raw.githubusercontent.com/igchoi/IBT610-CompGen/main/2025-Fall/msan/result/03_Signal_Peptide_Presence.png)

### 5. Surfaceome (2025.12.17)

## 📋 1. 분석 개요 (Overview)
이번 주 **SignalP 6.0**과 **DeepTMHMM**을 결합한 전략을 통해, 214개 균주의 전체 단백질 중 표면 단백질(Surfaceome) 후보군 **32,653개**를 최종 선별하고 구조적 특성에 따라 분류했습니다.

*   **Total Surfaceome Candidates:** 32,653 proteins
*   **Method:** Hybrid Classification (SignalP + DeepTMHMM)
*   **Goal:** 정확한 막 구조 예측을 통한 카테고리화

---

## 📊 2. 분석 결과 시각화 (Key Visualization)

분석 결과, Surfaceome은 **Complex Topology(복잡한 구조)**와 **Lipoprotein(지질단백질)**이 전체의 약 89%를 차지하는 것으로 나타났습니다.

![Integrated Surfaceome Analysis](https://github.com/igchoi/IBT610-CompGen/blob/main/2025-Fall/msan/result/%20Surfaceome_Final.png?raw=true)

---

## 🔬 3. 상세 분류 결과 (Detailed Statistics)

이미지 데이터를 바탕으로 집계된 최종 수치입니다. **Complex Topology** 그룹이 가장 큰 비중을 차지하며, 이는 다수의 TM Helix를 포함하거나 복합적인 신호를 가진 단백질군으로 해석됩니다.

| 순위 | 카테고리 (Category) | 단백질 수 (Count) | 비율 (Percentage) | 특징 및 해석 |
|:---:|:---|:---:|:---:|:---|
| **1** | **Complex Topology** | **17,040** | **52.2%** | 가장 큰 비중을 차지. 복잡한 막 관통 구조를 가지거나 다양한 위상(Topology)이 혼재된 그룹 |
| **2** | **Lipoprotein** | **11,918** | **36.5%** | SignalP가 LIPO로 예측한 그룹. 세포막 표면에 지질 앵커로 부착됨 (주요 백신 타겟) |
| **3** | **Multi-pass Membrane** | **2,842** | **8.7%** | 명확하게 여러 번 막을 관통하는 구조 (TM Helix 다수 보유) |
| **4** | **Other Surface Protein** | **853** | **2.6%** | 기타 표면 단백질로 분류된 소수 그룹 |
| **Total** | **Surfaceome** | **32,653** | **100%** | (분석된 214개 균주 합계) |

---

## 📂 4. 데이터셋 다운로드 (Data Assets)

분석에 사용된 Raw Data 및 교차 분석 테이블입니다.

*   📄 **Summary Crosstab (SignalP vs Category):** [signalp_vs_category_crosstab.csv](https://github.com/igchoi/IBT610-CompGen/blob/main/2025-Fall/msan/result/signalp_vs_category_crosstab.csv)
*   📄 **Full Classification List:** [surfaceome_hybrid_classified.csv](https://github.com/igchoi/IBT610-CompGen/blob/main/2025-Fall/msan/result/surfaceome_hybrid_classified.csv)

---

### 1. Directory Setup
스크립트 내 경로(`BASE_DIR`)를 사용자의 환경에 맞게 수정해야 합니다.
```text
/mnt/c/Users/AN/Desktop/gDrive/Study/Leacture/IBT610/
├── gff/                  # Input GFF3 files
├── Prokka_faa/           # Protein sequences (.faa)
├── Roary/                # Roary output csv
├── Result_sample/        # Genome assemblies (.fna)
├── Prokka_gbk/           # GenBank files
└── output_12600k_4060_v2/# (Created automatically)


