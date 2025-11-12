"""
Lactobacillus plantarum AAI 분석 파이프라인 (12600K 최적화)
- 스레드 경합 방지
- 타임아웃 30분 설정
- 진행 상황 실시간 모니터링
"""

import os
import sys
import shutil
import subprocess
import time
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram

# ========== 설정 ==========
# 사용자 경로 설정 (실제 경로로 수정 필요!)
FNA_DIR = r"C:\Users\AN\Desktop\gDrive\Study\Leacture\IBT610\Result_sample"
FAA_DIR = r"C:\Users\AN\Desktop\gDrive\Study\Leacture\IBT610\Prokka_faa"
OUTPUT_NAME = "output_lplantarum_aai"

# 작업 설정
GROUP_SIZE = 25
CPU_THREADS = 16
ANI_THRESHOLD = 95
TIMEOUT_SECONDS = 2700   # 30분 타임아웃 추가

# MLST 설정 (선택적)
MLST_SCHEME = "lplantarum"
MLST_TOOL = "mlst"  # 'mlst' 또는 None으로 비활성화

# Linux 경로 (WSL)
FNA_DIR_LINUX = "/mnt/c/Users/AN/Desktop/gDrive/Study/Leacture/IBT610/Result_sample"
FAA_DIR_LINUX = "/mnt/c/Users/AN/Desktop/gDrive/Study/Leacture/IBT610/Prokka_faa"

# ========== 경로 변환 ==========
def is_wsl():
    """WSL 환경인지 확인"""
    try:
        with open('/proc/version', 'r') as f:
            return 'microsoft' in f.read().lower()
    except:
        return False

def get_working_paths():
    """운영체제에 맞는 경로 반환"""
    if is_wsl() or os.path.exists('/mnt/c'):
        return FNA_DIR_LINUX, FAA_DIR_LINUX
    else:
        return FNA_DIR, FAA_DIR

# ========== 체크포인트 관리 ==========
def checkpoint_file(step_name, output_dir):
    """체크포인트 파일 경로 반환"""
    return os.path.join(output_dir, f"checkpoint_{step_name}.flag")

def is_step_done(step_name, output_dir):
    """해당 단계가 완료되었는지 확인"""
    cp_file = checkpoint_file(step_name, output_dir)
    return os.path.exists(cp_file)

def mark_step_done(step_name, output_dir):
    """단계 완료 표시"""
    cp_file = checkpoint_file(step_name, output_dir)
    with open(cp_file, 'w') as f:
        f.write(f"Completed at: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

# ========== FASTA 파일 검증 ==========
def validate_fasta_file(file_path):
    """FASTA 파일 형식 검증"""
    try:
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if not first_line.startswith('>'):
                return False, f"첫 줄이 '>'로 시작하지 않음: {first_line[:50]}"
            
            second_line = f.readline().strip()
            if not second_line:
                return False, "두 번째 줄이 비어있음"
            
            if len(second_line) < 10:
                return False, f"서열이 너무 짧음: {len(second_line)}bp"
        
        file_size = os.path.getsize(file_path)
        if file_size < 100:
            return False, f"파일 크기가 너무 작음: {file_size}bytes"
            
        return True, "OK"
    except Exception as e:
        return False, f"파일 읽기 실패: {str(e)}"

# ========== 이미 처리된 파일 확인 ==========
def get_processed_files(work_dir):
    """이미 처리된 파일 목록 확인"""
    processed = set()
    if not os.path.exists(work_dir):
        return processed
        
    for d in os.listdir(work_dir):
        if d.startswith("group_"):
            group_dir = os.path.join(work_dir, d)
            # 결과 파일 확인
            result_file = os.path.join(group_dir, "output_aai", "aai", "aai_summary.tsv")
            if os.path.exists(result_file) and os.path.getsize(result_file) > 100:
                # 이 그룹의 입력 파일들 목록
                for f in os.listdir(group_dir):
                    if f.endswith('.fasta'):
                        processed.add(f)
    return processed

# ========== MLST 도구 자동 탐지 (안전 버전) ==========
def find_mlst():
    """MLST 도구 자동 탐지 - 환경 충돌 시 None 반환"""
    print("\n" + "="*70)
    print("[MLST 도구 탐지]")
    print("="*70)
    
    # MLST 비활성화 설정
    if MLST_TOOL is None:
        print("⚠️  MLST_TOOL이 None으로 설정됨 → MLST 분석 건너뜀")
        return None
    
    # which 명령어로 찾기
    try:
        result = subprocess.run(['which', 'mlst'], capture_output=True, text=True, check=True)
        mlst_path = result.stdout.strip()
        if os.path.exists(mlst_path):
            print(f"✓ MLST 발견: {mlst_path}")
            
            # 스킴 확인
            try:
                scheme_result = subprocess.run(['mlst', '--list'], capture_output=True, text=True, timeout=10)
                available_schemes = scheme_result.stdout.strip().split('\n')
                if MLST_SCHEME in available_schemes:
                    print(f"✓ 스킴 '{MLST_SCHEME}' 사용 가능")
                    return mlst_path
                else:
                    print(f"⚠️  스킴 '{MLST_SCHEME}' 없음")
                    return None
            except:
                print("⚠️  MLST 스킴 확인 실패")
                return None
    except:
        print("❌ MLST 도구를 찾을 수 없습니다!")
        print("   설치 방법 (Python 3.13 충돌 주의):")
        print("   → 별도 환경 생성: conda create -n mlst_env python=3.11 mlst")
        print("   → 활성화: conda activate mlst_env")
        print("   → 수동 실행: mlst --scheme lplantarum *.fasta")
        return None

# ========== 1단계: Prokka-FNA 매칭 ==========
def match_prokka_fna(output_dir):
    """Prokka .fasta 파일 매칭 (통합 파일만)"""
    step_name = "01_matching"
    if is_step_done(step_name, output_dir):
        match_file = os.path.join(output_dir, "matched_samples.tsv")
        if os.path.exists(match_file):
            return pd.read_csv(match_file, sep="\t"), *get_working_paths()
        return None, None, None
    
    print("\n" + "="*70)
    print("[1단계] Prokka-FNA 매칭 (통합 파일만)")
    print("="*70)
    
    fna_dir, faa_dir = get_working_paths()
    
    if not os.path.exists(fna_dir):
        print(f"❌ FNA 디렉토리 없음: {fna_dir}")
        return None, None, None
    
    if not os.path.exists(faa_dir):
        print(f"❌ FAA 디렉토리 없음: {faa_dir}")
        return None, None, None
    
    fna_files = [f for f in os.listdir(fna_dir) if f.endswith('.fna')]
    faa_files = [f for f in os.listdir(faa_dir) if f.endswith('.fasta') and '_part' not in f]
    
    print(f"\n📂 파일 탐색:")
    print(f"   - FNA 파일: {len(fna_files)}개")
    print(f"   - FAA 파일 (통합만): {len(faa_files)}개")
    
    # GCF/GCA ID 추출
    fna_gcf_map = {}
    for f in fna_files:
        if f.startswith('GCF_') or f.startswith('GCA_'):
            gcf_id = f.split('.')[0]
            fna_gcf_map[gcf_id] = f
    
    faa_gcf_map = {}
    for f in faa_files:
        if f.startswith('GCF_') or f.startswith('GCA_'):
            prefix = f.replace('.fasta', '')
            faa_gcf_map[prefix] = f
    
    print(f"   - FNA ID: {len(fna_gcf_map)}개")
    print(f"   - FAA ID: {len(faa_gcf_map)}개")
    
    # 매칭 및 검증
    matched = []
    validation_errors = []
    for fna_gcf, fna_file in fna_gcf_map.items():
        if fna_gcf in faa_gcf_map:
            faa_file = faa_gcf_map[fna_gcf]
            faa_path = os.path.join(faa_dir, faa_file)
            
            if not os.path.exists(faa_path):
                validation_errors.append(f"{faa_file}: 파일이 존재하지 않음")
                continue
            
            is_valid, msg = validate_fasta_file(faa_path)
            if is_valid:
                matched.append({
                    'GCF_ID': fna_gcf,
                    'FNA_file': fna_file,
                    'FAA_file': faa_file
                })
            else:
                validation_errors.append(f"{faa_file}: {msg}")
    
    if validation_errors:
        print(f"\n⚠️  검증 오류 ({len(validation_errors)}개):")
        for err in validation_errors[:5]:
            print(f"   - {err}")
        if len(validation_errors) > 5:
            print(f"   ... 외 {len(validation_errors)-5}개")
    
    if not matched:
        print("\n❌ 매칭되고 검증된 파일이 없습니다!")
        return None, None, None
    
    matched_df = pd.DataFrame(matched)
    os.makedirs(output_dir, exist_ok=True)
    
    match_file = os.path.join(output_dir, "matched_samples.tsv")
    matched_df.to_csv(match_file, sep="\t", index=False)
    
    print(f"\n✅ 매칭 완료: {len(matched)}개 파일")
    print(f"📄 저장: {match_file}")
    
    mark_step_done(step_name, output_dir)
    return matched_df, fna_dir, faa_dir

# ========== 2단계: FAA 파일 복사 ==========
def prepare_faa_working_directory(matched_df, faa_dir, fna_dir):
    """매칭된 FAA 파일만 작업 디렉토리에 복사"""
    output_dir = os.path.join(fna_dir, OUTPUT_NAME)
    step_name = "02_copy_faa"
    if is_step_done(step_name, output_dir):
        work_dir = os.path.join(faa_dir, "matched_prokka_faa")
        return work_dir if os.path.exists(work_dir) else None
    
    print("\n" + "="*70)
    print("[2단계] 매칭된 FAA 파일 복사")
    print("="*70)
    
    work_dir = os.path.join(faa_dir, "matched_prokka_faa")
    os.makedirs(work_dir, exist_ok=True)
    
    print(f"\n📁 작업 디렉토리: {work_dir}")
    print(f"🔄 {len(matched_df)}개 파일 복사 중...")
    
    copied = 0
    for idx, row in matched_df.iterrows():
        src = os.path.join(faa_dir, row['FAA_file'])
        dst = os.path.join(work_dir, row['FAA_file'])
        
        if os.path.exists(src):
            if not os.path.exists(dst):
                shutil.copy2(src, dst)
            copied += 1
        else:
            print(f"   ⚠️  파일 없음: {src}")
    
    print(f"✅ 복사 완료: {copied}개 파일")
    mark_step_done(step_name, output_dir)
    
    return work_dir

# ========== CompareM 자동 탐지 ==========
def find_comparem():
    """CompareM 실행 파일 경로 자동 탐지"""
    print("\n" + "="*70)
    print("[CompareM 탐지]")
    print("="*70)
    
    # 1. 현재 conda 환경
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if conda_prefix:
        possible_path = os.path.join(conda_prefix, 'bin', 'comparem')
        if os.path.exists(possible_path):
            print(f"✓ 발견 (conda): {possible_path}")
            return possible_path
    
    # 2. which 명령어
    try:
        result = subprocess.run(['which', 'comparem'], capture_output=True, text=True, check=True)
        comparem_path = result.stdout.strip()
        if os.path.exists(comparem_path):
            print(f"✓ 발견: {comparem_path}")
            return comparem_path
    except:
        pass
    
    # 3. PATH에서 찾기
    comparem_path = shutil.which('comparem')
    if comparem_path:
        print(f"✓ 발견 (PATH): {comparem_path}")
        return comparem_path
    
    print("❌ CompareM을 찾을 수 없습니다!")
    print("   설치: conda install -c bioconda comparem")
    return None

# ========== 3단계: 그룹별 AAI 분석 (개선) ==========
def group_aai_analysis(comparem_path, work_dir, output_dir):
    """그룹별 AAI 분석 (타임아웃 + 스레드 경합 방지)"""
    step_name = "03_aai_analysis"
    if is_step_done(step_name, output_dir):
        return True, work_dir
    
    print("\n" + "="*70)
    print("[3단계] 그룹별 AAI 분석")
    print("="*70)
    
    faa_files = [f for f in os.listdir(work_dir) if f.endswith('.fasta')]
    total_files = len(faa_files)
    
    if total_files == 0:
        print("❌ 작업 디렉토리에 FAA 파일이 없습니다!")
        return False, work_dir
    
    processed_files = get_processed_files(work_dir)
    remaining_files = [f for f in faa_files if f not in processed_files]
    
    print(f"\n📊 처리 상태:")
    print(f"   - 전체 파일: {total_files}개")
    print(f"   - 이미 처리됨: {len(processed_files)}개")
    print(f"   - 남은 파일: {len(remaining_files)}개")
    
    if not remaining_files:
        print("✅ 모든 파일이 이미 처리되었습니다!")
        mark_step_done(step_name, output_dir)
        return True, work_dir
    
    groups = [remaining_files[i:i+GROUP_SIZE] for i in range(0, len(remaining_files), GROUP_SIZE)]
    total_groups = len(groups)
    
    print(f"\n📊 분석 설정:")
    print(f"   - 처리할 파일: {len(remaining_files)}개")
    print(f"   - 그룹 수: {total_groups}개")
    print(f"   - 그룹 크기: {GROUP_SIZE}개")
    print(f"   - CPU 스레드: {CPU_THREADS}개")
    print(f"   - 타임아웃: {TIMEOUT_SECONDS//60}분")
    
    success_count = 0
    failed_groups = []
    
    for gi, files in enumerate(groups, 1):
        group_dir = os.path.join(work_dir, f"group_{gi}")
        os.makedirs(group_dir, exist_ok=True)
        
        # 파일 복사
        for fname in files:
            src = os.path.join(work_dir, fname)
            dst = os.path.join(group_dir, fname)
            if not os.path.exists(dst):
                shutil.copy2(src, dst)
        
        aai_out = os.path.join(group_dir, "output_aai")
        os.makedirs(aai_out, exist_ok=True)
        
        print(f"\n[그룹 {gi}/{total_groups}] 실행 중...")
        print(f"   ├── 파일: {len(files)}개 (첫 파일: {files[0]})")
        print(f"   ├── CPU: {CPU_THREADS} threads")
        ts = time.time()
        
        cmd = [
            comparem_path, "aai_wf",
            "--proteins", "--file_ext", "fasta",
            "--cpus", str(CPU_THREADS),
            group_dir, aai_out
        ]
        
        try:
            # CPU 사용량 모니터링 백그라운드 프로세스 시작
            monitor_process = subprocess.Popen([
                'bash', '-c', 
                f'while true; do top -bn1 | grep "Cpu(s)"; sleep 10; done'
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=TIMEOUT_SECONDS
            )
            
            # 모니터링 종료
            monitor_process.terminate()
            
            result_file = os.path.join(aai_out, "aai", "aai_summary.tsv")
            if os.path.exists(result_file) and os.path.getsize(result_file) > 100:
                te = time.time()
                elapsed_min = (te - ts) / 60
                print(f"   └── ✓ 성공 ({elapsed_min:.1f}분)")
                success_count += 1
            else:
                print(f"   └── ⚠️  결과 없음 (코드: {result.returncode})")
                failed_groups.append(gi)
                
        except subprocess.TimeoutExpired:
            monitor_process.terminate()
            print(f"   └── ✗ 타임아웃 ({TIMEOUT_SECONDS//60}분 초과)")
            failed_groups.append(gi)
        except Exception as e:
            try:
                monitor_process.terminate()
            except:
                pass
            print(f"   └── ✗ 오류: {str(e)}")
            failed_groups.append(gi)
    
    print(f"\n{'─'*70}")
    print(f"📊 그룹별 분석: {success_count}/{total_groups} 성공")
    if failed_groups:
        print(f"⚠️  실패 그룹: {', '.join(map(str, failed_groups))}")
    print(f"{'─'*70}")
    
    if success_count > 0:
        mark_step_done(step_name, output_dir)
    
    return success_count > 0, work_dir

# ========== 4단계: 결과 통합 ==========
def merge_aai_results(work_dir, output_dir):
    """AAI 결과 통합"""
    step_name = "04_merge_results"
    if is_step_done(step_name, output_dir):
        merge_file = os.path.join(output_dir, "aai_summary_merged.tsv")
        if os.path.exists(merge_file):
            return pd.read_csv(merge_file, sep="\t")
        return None
    
    print("\n" + "="*70)
    print("[4단계] AAI 결과 통합")
    print("="*70)
    
    result_files = []
    for d in os.listdir(work_dir):
        if d.startswith("group_"):
            result_file = os.path.join(work_dir, d, "output_aai", "aai", "aai_summary.tsv")
            if os.path.exists(result_file) and os.path.getsize(result_file) > 100:
                result_files.append(result_file)
    
    if not result_files:
        print("❌ 통합할 결과 없음")
        return None
    
    print(f"\n📂 발견된 결과: {len(result_files)}개 그룹")
    
    aai_parts = []
    for idx, rf in enumerate(result_files, 1):
        try:
            df = pd.read_csv(rf, sep="\t")
            df.columns = [c.strip() for c in df.columns]
            
            g1_col = [c for c in df.columns if any(x in c for x in ["Genome A", "genome1", "Genome_A"])][0]
            g2_col = [c for c in df.columns if any(x in c for x in ["Genome B", "genome2", "Genome_B"])][0]
            aai_col = [c for c in df.columns if any(x in c for x in ["Mean AAI", "aai", "AAI"])][0]
            
            df2 = df[[g1_col, g2_col, aai_col]].copy()
            df2.columns = ["Genome_A", "Genome_B", "AAI"]
            aai_parts.append(df2)
            
        except Exception as e:
            print(f"[{idx}] ✗ 실패: {str(e)}")
    
    if not aai_parts:
        print("\n❌ 읽을 수 있는 결과 없음")
        return None
    
    all_aai = pd.concat(aai_parts, axis=0).drop_duplicates()
    merge_file = os.path.join(output_dir, "aai_summary_merged.tsv")
    all_aai.to_csv(merge_file, sep="\t", index=False)
    
    print(f"\n✅ 통합 완료:")
    print(f"   - 총 비교 쌍: {len(all_aai):,}개")
    print(f"   - 파일: {merge_file}")
    print(f"   - AAI 범위: {all_aai['AAI'].min():.2f}% ~ {all_aai['AAI'].max():.2f}%")
    
    mark_step_done(step_name, output_dir)
    return all_aai

# ========== 5단계: 매트릭스 및 시각화 ==========
def build_matrix_and_visualize(aai_df, output_dir):
    """매트릭스 생성 및 시각화"""
    step_name = "05_matrix_viz"
    if is_step_done(step_name, output_dir):
        mat_file = os.path.join(output_dir, "aai_matrix.csv")
        if os.path.exists(mat_file):
            return pd.read_csv(mat_file, index_col=0)
        return None
    
    print("\n" + "="*70)
    print("[5단계] 매트릭스 및 시각화")
    print("="*70)
    
    names = sorted(set(aai_df["Genome_A"]).union(set(aai_df["Genome_B"])))
    print(f"\n📊 매트릭스: {len(names)} × {len(names)}")
    
    mat = pd.DataFrame(np.nan, index=names, columns=names)
    for _, row in aai_df.iterrows():
        mat.at[row["Genome_A"], row["Genome_B"]] = row["AAI"]
        mat.at[row["Genome_B"], row["Genome_A"]] = row["AAI"]
    
    np.fill_diagonal(mat.values, 100)
    mat = mat.fillna(100)
    
    # 매트릭스 저장
    mat_file = os.path.join(output_dir, "aai_matrix.csv")
    mat.to_csv(mat_file)
    print(f"✓ 매트릭스 저장: {mat_file}")
    
    # 히트맵
    print("\n🎨 히트맵 생성...")
    try:
        plt.figure(figsize=(max(12, len(names)*0.3), max(10, len(names)*0.3)))
        sns.heatmap(mat, cmap='viridis', square=True, cbar_kws={'label': 'AAI (%)'})
        plt.title("Lactobacillus plantarum AAI Heatmap", fontsize=16, pad=20)
        plt.tight_layout()
        
        heatmap_file = os.path.join(output_dir, "aai_heatmap.png")
        plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ 히트맵: {heatmap_file}")
    except Exception as e:
        print(f"✗ 히트맵 실패: {str(e)}")
    
    # 덴드로그램
    print("\n🌳 덴드로그램 생성...")
    try:
        dist = 100 - mat.values
        np.fill_diagonal(dist, 0)
        Z = linkage(squareform(dist, checks=False), method='average')
        
        plt.figure(figsize=(max(16, len(names)*0.2), 8))
        dendrogram(Z, labels=mat.index, leaf_rotation=90, leaf_font_size=8)
        plt.title("Lactobacillus plantarum AAI Dendrogram", fontsize=16, pad=20)
        plt.ylabel("Distance (100 - AAI%)", fontsize=12)
        plt.tight_layout()
        
        dendro_file = os.path.join(output_dir, "aai_dendrogram.png")
        plt.savefig(dendro_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ 덴드로그램: {dendro_file}")
    except Exception as e:
        print(f"✗ 덴드로그램 실패: {str(e)}")
    
    mark_step_done(step_name, output_dir)
    return mat

# ========== 6단계: Low ANI 필터링 ==========
def filter_low_ani(mat, output_dir):
    """Low ANI 필터링"""
    step_name = "06_filter_low_ani"
    if is_step_done(step_name, output_dir):
        return
    
    print("\n" + "="*70)
    print(f"[6단계] Low ANI 필터링 (< {ANI_THRESHOLD}%)")
    print("="*70)
    
    low_pairs = []
    for i in mat.index:
        for j in mat.columns:
            if i < j and mat.at[i, j] < ANI_THRESHOLD:
                low_pairs.append({"Genome_A": i, "Genome_B": j, "AAI": mat.at[i, j]})
    
    low_file = os.path.join(output_dir, "low_ani_filtered.tsv")
    
    if low_pairs:
        df_low = pd.DataFrame(low_pairs).sort_values("AAI")
        df_low.to_csv(low_file, sep="\t", index=False)
        print(f"\n✓ Low ANI 쌍: {len(df_low)}개")
        print(f"📄 저장: {low_file}")
        for _, row in df_low.head(5).iterrows():
            print(f"   - {row['Genome_A']}: {row['Genome_B']}: {row['AAI']:.2f}%")
    else:
        pd.DataFrame(columns=["Genome_A", "Genome_B", "AAI"]).to_csv(low_file, sep="\t", index=False)
        print(f"\n✓ Low ANI 쌍 없음 (모두 {ANI_THRESHOLD}% 이상)")
    
    mark_step_done(step_name, output_dir)

# ========== 7단계: 자동 메타데이터 생성 ==========
def generate_metadata(base_dir, output_dir, matched_df):
    """파일명 기반 메타데이터 자동 생성"""
    step_name = "07_metadata"
    if is_step_done(step_name, output_dir):
        meta_file = os.path.join(output_dir, "metadata.tsv")
        if os.path.exists(meta_file):
            return pd.read_csv(meta_file, sep="\t")
        return None
    
    print("\n" + "="*70)
    print("[7단계] 자동 메타데이터 생성")
    print("="*70)
    
    # 파일명에서 정보 추출
    metadata = []
    for _, row in matched_df.iterrows():
        gcf_id = row['GCF_ID']
        fna_file = row['FNA_file']
        
        # GCF ID 파싱
        parts = gcf_id.split('_')
        assembly_version = parts[1] if len(parts) > 1 else "unknown"
        
        # 기본 메타데이터 생성
        metadata.append({
            "Sample": gcf_id,
            "FNA_File": fna_file,
            "FAA_File": row['FAA_file'],
            "Assembly_Version": assembly_version,
            "Species": "Lactobacillus plantarum",
            "Source": "NCBI RefSeq",
            "Created_Date": time.strftime("%Y-%m-%d")
        })
    
    meta_df = pd.DataFrame(metadata)
    meta_file = os.path.join(output_dir, "metadata.tsv")
    meta_df.to_csv(meta_file, sep="\t", index=False)
    
    print(f"✓ 메타데이터 생성: {len(meta_df)}개 샘플")
    print(f"📄 저장: {meta_file}")
    
    # 요약
    try:
        meta_sum = meta_df.describe(include="all").T
        meta_sum.to_csv(os.path.join(output_dir, "metadata_summary.tsv"), sep="\t")
        print(f"✓ 메타데이터 요약: metadata_summary.tsv")
    except Exception as e:
        print(f"⚠️  요약 실패: {str(e)}")
    
    mark_step_done(step_name, output_dir)
    return meta_df

# ========== 8단계: 자동 MLST 분석 (안전 버전) ==========
def run_mlst_analysis(base_dir, output_dir, work_dir):
    """MLST 도구로 자동 분석 (실패 시에도 계속 진행)"""
    step_name = "08_mlst"
    if is_step_done(step_name, output_dir):
        mlst_file = os.path.join(output_dir, "mlst_results.tsv")
        if os.path.exists(mlst_file):
            return pd.read_csv(mlst_file, sep="\t")
        return None
    
    # MLST 도구 확인 (없으면 체크포인트만 생성하고 리턴)
    mlst_path = find_mlst()
    if not mlst_path:
        print("\n⚠️  MLST 분석을 건너뜁니다")
        mark_step_done(step_name, output_dir)
        return None
    
    print("\n" + "="*70)
    print("[8단계] 자동 MLST 분석")
    print("="*70)
    
    # FAA 파일 확인
    if not os.path.exists(work_dir):
        print(f"❌ 작업 디렉토리 없음: {work_dir}")
        mark_step_done(step_name, output_dir)
        return None
    
    fasta_files = [f for f in os.listdir(work_dir) if f.endswith('.fasta')]
    if not fasta_files:
        print("❌ 분석할 FAA 파일 없음")
        mark_step_done(step_name, output_dir)
        return None
    
    # MLST 실행
    mlst_results = []
    print(f"\n🧬 MLST 분석 시작 ({len(fasta_files)}개 파일)...")
    
    for idx, fname in enumerate(fasta_files, 1):
        fpath = os.path.join(work_dir, fname)
        prefix = fname.replace('.fasta', '')
        
        try:
            cmd = ['mlst', '--scheme', MLST_SCHEME, '--quiet', fpath]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            
            if result.returncode == 0:
                # 결과 파싱
                parts = result.stdout.strip().split('\t')
                if len(parts) >= 3:
                    mlst_results.append({
                        "Sample": prefix,
                        "Scheme": parts[0],
                        "ST": parts[1],
                        "Allele_Profile": "\t".join(parts[2:])
                    })
                    print(f"   [{idx}/{len(fasta_files)}] ✓ {prefix}: ST={parts[1]}")
                else:
                    print(f"   [{idx}/{len(fasta_files)}] ⚠️  {prefix}: 형식 오류")
            else:
                print(f"   [{idx}/{len(fasta_files)}] ⚠️  {prefix}: 실패 (반환코드: {result.returncode})")
                
        except subprocess.TimeoutExpired:
            print(f"   [{idx}/{len(fasta_files)}] ⚠️  {prefix}: 타임아웃 (120초)")
        except Exception as e:
            print(f"   [{idx}/{len(fasta_files)}] ⚠️  {prefix}: {str(e)}")
    
    # 결과 저장
    if mlst_results:
        mlst_df = pd.DataFrame(mlst_results)
        mlst_file = os.path.join(output_dir, "mlst_results.tsv")
        mlst_df.to_csv(mlst_file, sep="\t", index=False)
        
        print(f"\n✓ MLST 완료: {len(mlst_results)}개 샘플")
        print(f"📄 저장: {mlst_file}")
        
        # ST 요약
        st_counts = mlst_df['ST'].value_counts()
        st_summary_file = os.path.join(output_dir, "st_summary.tsv")
        st_counts.to_csv(st_summary_file, sep="\t", header=['Count'])
        print(f"📄 ST 요약: {st_summary_file}")
        
        # ST별 샘플 목록
        st_samples = mlst_df.groupby('ST')['Sample'].apply(list)
        st_samples_file = os.path.join(output_dir, "st_samples.tsv")
        with open(st_samples_file, 'w') as f:
            f.write("ST\tSamples\tCount\n")
            for st, samples in st_samples.items():
                f.write(f"{st}\t{','.join(samples)}\t{len(samples)}\n")
        print(f"📄 ST별 샘플: {st_samples_file}")
        
    else:
        print("\n⚠️  MLST 결과가 없습니다")
    
    mark_step_done(step_name, output_dir)
    return mlst_df if mlst_results else None

# ========== 메인 실행 ==========
def main():
    """메인 실행"""
    print("\n" + "="*70)
    print("  Lactobacillus plantarum AAI 분석 파이프라인")
    print("  (환경 충돌 해결: MLST 선택적 실행)")
    print("="*70)
    print(f"  Python: {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}")
    print(f"  MLST 스킴: {MLST_SCHEME}")
    print(f"  ANI 임계값: {ANI_THRESHOLD}%")
    print(f"  MLST 도구: {'사용' if MLST_TOOL else '비활성화'}")
    print("="*70)
    
    start_time = time.time()
    
    # 출력 디렉토리
    fna_dir, faa_dir = get_working_paths()
    output_dir = os.path.join(fna_dir, OUTPUT_NAME)
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\n📁 출력 디렉토리: {output_dir}")
    
    try:
        # 1단계: 매칭
        matched_df, _, _ = match_prokka_fna(output_dir)
        if matched_df is None:
            print("\n❌ 매칭 실패로 종료")
            sys.exit(1)
        
        work_dir = os.path.join(faa_dir, "matched_prokka_faa")
        
        # 2단계: 파일 복사
        work_dir = prepare_faa_working_directory(matched_df, faa_dir, fna_dir)
        if work_dir is None:
            print("\n❌ 파일 복사 실패")
            sys.exit(1)
        
        # 3단계: CompareM 찾기
        comparem_path = find_comparem()
        if not comparem_path:
            print("\n❌ CompareM 없음")
            sys.exit(1)
        
        # 4단계: AAI 분석
        success, work_dir = group_aai_analysis(comparem_path, work_dir, output_dir)
        if not success:
            print("\n⚠️  AAI 분석 부분 실패, 결과 통합 시도")
        
        # 5단계: 결과 통합
        aai_df = merge_aai_results(work_dir, output_dir)
        if aai_df is None:
            print("\n❌ 통합할 결과 없음")
            sys.exit(1)
        
        # 6단계: 시각화
        mat = build_matrix_and_visualize(aai_df, output_dir)
        
        # 7단계: Low ANI 필터링
        filter_low_ani(mat, output_dir)
        
        # 8단계: 자동 메타데이터 생성
        generate_metadata(fna_dir, output_dir, matched_df)
        
        # 9단계: 자동 MLST 분석 (선택적)
        if MLST_TOOL:
            try:
                mlst_result = run_mlst_analysis(fna_dir, output_dir, work_dir)
                if mlst_result is None:
                    print("\n⚠️  MLST 결과 없음 → 주 분석은 정상 완료")
            except Exception as e:
                print(f"\n⚠️  MLST 분석 오류 (무시): {str(e)}")
        
        # 완료 요약
        elapsed = time.time() - start_time
        mins, secs = divmod(int(elapsed), 60)
        
        print("\n" + "="*70)
        print("  🎉 분석 완료!")
        print("="*70)
        print(f"\n⏱️  실행 시간: {mins}분 {secs}초")
        print(f"📁 결과 위치: {output_dir}")
        print(f"\n📊 처리 샘플: {len(matched_df)}개")
        
        print(f"\n{'─'*70}")
        print("📂 생성된 핵심 파일")
        print(f"{'─'*70}")
        
        core_files = [
            "matched_samples.tsv",
            "aai_summary_merged.tsv",
            "aai_matrix.csv",
            "aai_heatmap.png",
            "aai_dendrogram.png",
            "metadata.tsv"
        ]
        
        for cf in core_files:
            full_path = os.path.join(output_dir, cf)
            if os.path.exists(full_path):
                print(f"   ✓ {cf}")
            else:
                print(f"   ✗ {cf}")
        
        print(f"\n{'─'*70}")
        print("🔧 MLST 관련 (선택적)")
        print(f"{'─'*70}")
        print(f"   ✓ low_ani_filtered.tsv")
        if MLST_TOOL and os.path.exists(os.path.join(output_dir, "mlst_results.tsv")):
            print(f"   ✓ mlst_results.tsv")
            print(f"   ✓ st_summary.tsv")
        else:
            print(f"   ⚠ mlst_results.tsv (MLST 도구 없음)")
            print(f"\n💡 MLST 분석을 원한다면:")
            print(f"   conda create -n mlst_env python=3.11 mlst")
            print(f"   conda activate mlst_env")
            print(f"   cd {work_dir}")
            print(f"   mlst --scheme {MLST_SCHEME} *.fasta > {os.path.join(output_dir, 'mlst_results.tsv')}")
        
        print(f"\n{'─'*70}")
        print("✅ 주 분석 (AAI) 정상 완료!")
        print("="*70)
        
    except KeyboardInterrupt:
        print("\n\n⚠️  사용자 중단")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ 에러: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
