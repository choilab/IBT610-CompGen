import os
import glob
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations
from collections import defaultdict

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO

# --- 사용자 설정 (USER CONFIGURATION) ---
# fastANI 실행 파일의 절대 경로를 입력하세요. (필수)
# 예: "C:/Users/user/fastani.exe" 또는 "/home/user/bin/fastANI"
FASTANI_PATH = "/root/bin/fastANI"

# 입력 및 출력 폴더 설정
INPUT_DIR = os.path.join("input", "genomes")
OUTPUT_DIR = "output"

# 사용할 CPU 코어 수 (0이면 가능한 모든 코어 사용)
CPU_CORES = os.cpu_count() or 1

# --- 스크립트 내부 상수 (DO NOT EDIT BELOW) ---
QUALITY_FILE = os.path.join(OUTPUT_DIR, "genome_quality.csv")
ANI_FILE = os.path.join(OUTPUT_DIR, "all_pairs.tsv")
PROGRESS_FILE = os.path.join(OUTPUT_DIR, 'progress.txt')

# --- 헬퍼 함수: 한글 폰트 설정 ---
def set_korean_font():
    """matplotlib에서 한글이 깨지지 않도록 폰트를 설정합니다."""
    try:
        import matplotlib.font_manager as fm
        font_paths = [
            "C:/Windows/Fonts/malgun.ttf",                        # Windows
            "/usr/share/fonts/truetype/nanum/NanumGothic.ttf",   # Ubuntu
            "/System/Library/Fonts/AppleGothic.ttf"              # macOS
        ]
        for fp in font_paths:
            if os.path.exists(fp):
                plt.rcParams['font.family'] = fm.FontProperties(fname=fp).get_name()
                return
        # 적절한 폰트가 없으면 기본값으로 유지
    except Exception as e:
        print(f"폰트 설정 중 오류 발생: {e} (기본 폰트로 진행)")

# --- 1단계 함수: 게놈 품질 평가 ---
def calculate_assembly_stats(fasta_file):
    """주어진 FASTA 파일로부터 N50, L50, GC 함량 등 통계를 계산합니다."""
    lengths, gc_count, total_len = [], 0, 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        l = len(seq)
        lengths.append(l)
        gc_count += seq.count('G') + seq.count('C')
        total_len += l
        
    if not lengths:
        return {"N50": 0, "L50": 0, "Total_Length": 0, "Seq_Count": 0, "GC_Content": 0}

    lengths.sort(reverse=True)
    cum_sum = np.cumsum(lengths)
    
    half_total = total_len / 2
    n50_val = next((l for i, l in enumerate(lengths) if cum_sum[i] >= half_total), 0)
    l50_val = next((i + 1 for i, l in enumerate(lengths) if cum_sum[i] >= half_total), 0)
    gc_content = (gc_count / total_len * 100) if total_len else 0
    
    return {
        "N50": n50_val, "L50": l50_val, "Total_Length": total_len,
        "Seq_Count": len(lengths), "GC_Content": round(gc_content, 2)
    }

# --- 2단계 함수: 병렬 fastANI 실행 ---
def print_progress(done, total, start_time):
    """콘솔과 파일에 진행 상황을 출력합니다."""
    elapsed = time.time() - start_time
    speed = done / elapsed if elapsed > 0 else 0
    eta = (total - done) / speed if speed > 0 else 0
    pct = int(100 * done / total)
    
    elapsed_str = f"{int(elapsed//60)}분 {int(elapsed%60)}초"
    eta_str = f"{int(eta//60)}분 {int(eta%60)}초"
    
    msg = f"[{pct:3d}%] ({done}/{total}) | 경과: {elapsed_str} | ETA: {eta_str}"
    print(msg, end="\r", flush=True)
    
    with open(PROGRESS_FILE, 'w', encoding='utf-8') as f:
        f.write(f"진행률: {pct}% ({done}/{total})\n")
        f.write(f"경과 시간: {elapsed_str}\n")
        f.write(f"예상 남은 시간: {eta_str}\n")

def run_fastani_parallel(fna_files, max_workers=CPU_CORES):
    """ThreadPoolExecutor를 사용하여 fastANI를 병렬로 실행합니다."""
    pairs = list(combinations(fna_files, 2))
    total_pairs = len(pairs)
    results = []
    
    print(f"\n[2/4] fastANI 병렬 분석 시작 (총 {total_pairs}쌍, {max_workers}코어 사용)")
    start_time = time.time()

    def fastani_worker(query, ref):
        """단일 fastANI 작업을 수행하는 함수."""
        q_name = os.path.splitext(os.path.basename(query))[0]
        r_name = os.path.splitext(os.path.basename(ref))[0]
        tmp_out = os.path.join(OUTPUT_DIR, f"tmp_{q_name}_vs_{r_name}.txt")
        cmd = [FASTANI_PATH, "--query", query, "--ref", ref, "--threads", "1", "--output", tmp_out]
        
        try:
            subprocess.run(cmd, check=True, timeout=600, capture_output=True, text=True)
            if os.path.exists(tmp_out):
                with open(tmp_out, 'r') as f:
                    line = f.readline()
                    if line:
                        parts = line.strip().split("\t")
                        if len(parts) >= 3:
                            return [q_name, r_name, float(parts[2])]
        except subprocess.CalledProcessError as e:
            # 오류 발생 시 stderr를 출력하여 디버깅에 도움
            print(f"\nfastANI 오류: {q_name} vs {r_name}\n{e.stderr}")
        except Exception as e:
            print(f"\nfastANI 작업 실패: {q_name} vs {r_name} ({str(e)})")
        finally:
            if os.path.exists(tmp_out):
                os.remove(tmp_out)
        return None

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(fastani_worker, q, r) for q, r in pairs}
        for i, future in enumerate(as_completed(futures)):
            result = future.result()
            if result:
                results.append(result)
            print_progress(i + 1, total_pairs, start_time)
            
    print("\n✓ fastANI 분석 완료")
    return pd.DataFrame(results, columns=['Genome1', 'Genome2', 'ANI'])

# --- 3단계 함수: 시각화 ---
def save_heatmap(matrix, file_path, title):
    """주어진 행렬을 히트맵으로 저장합니다."""
    plt.figure(figsize=(min(30, max(12, len(matrix)//2)), min(30, max(12, len(matrix)//2))))
    sns.heatmap(matrix.astype(float), cmap='viridis', annot=False)
    plt.title(title, fontsize=16)
    plt.tight_layout()
    plt.savefig(file_path, dpi=300)
    plt.close()

# --- 4단계 함수: 클러스터링 및 요약 ---
def build_clusters(df_ani, threshold):
    """ANI 값을 기준으로 그래프를 만들고 연결된 컴포넌트(클러스터)를 찾습니다."""
    edges = defaultdict(set)
    for _, row in df_ani.iterrows():
        if row['ANI'] >= threshold:
            edges[row['Genome1']].add(row['Genome2'])
            edges[row['Genome2']].add(row['Genome1'])
            
    visited = set()
    clusters = []
    all_nodes = set(df_ani['Genome1']).union(df_ani['Genome2'])
    
    for node in all_nodes:
        if node in visited:
            continue
        
        group = set()
        queue = [node]
        visited.add(node)
        
        while queue:
            current_node = queue.pop(0)
            group.add(current_node)
            for neighbor in edges.get(current_node, set()):
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        
        if group:
            clusters.append(sorted(list(group)))
            
    return clusters

def get_representative(cluster, quality_scores):
    """클러스터 내에서 품질 점수가 가장 높은 멤버를 대표로 선정합니다."""
    if not quality_scores:
        return cluster[0]
    return max(cluster, key=lambda genome: quality_scores.get(genome, 0))

# --- 메인 실행 함수 ---
def main():
    """파이프라인의 모든 단계를 순차적으로 실행합니다."""
    try:
        # 0. 준비 단계: 폴더 생성 및 폰트 설정
        os.makedirs(INPUT_DIR, exist_ok=True)
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        set_korean_font()
        
        fasta_files = glob.glob(os.path.join(INPUT_DIR, '*.fna')) + glob.glob(os.path.join(INPUT_DIR, '*.fasta'))
        if not fasta_files:
            print(f"오류: '{INPUT_DIR}' 폴더에 게놈 파일(*.fna, *.fasta)이 없습니다.")
            return

        # 1. 게놈 품질 평가
        print("[1/4] 게놈 품질 평가 시작...")
        stats_list = []
        for i, f in enumerate(fasta_files):
            genome_name = os.path.splitext(os.path.basename(f))[0]
            stats = calculate_assembly_stats(f)
            stats['Genome'] = genome_name
            stats_list.append(stats)
            print(f"  ({i+1}/{len(fasta_files)}) {genome_name} | N50={stats['N50']:,}")
        
        df_quality = pd.DataFrame(stats_list)
        # Z-score 기반 품질 점수 계산
        for col in ['N50', 'Total_Length']:
            if df_quality[col].std() > 0:
                df_quality[f'{col}_Z'] = (df_quality[col] - df_quality[col].mean()) / df_quality[col].std()
            else:
                df_quality[f'{col}_Z'] = 0
        df_quality['Quality_Score'] = (df_quality['N50_Z'] + 0.5 * df_quality['Total_Length_Z'] - 0.2 * df_quality['L50']).rank(ascending=False)
        df_quality.to_csv(QUALITY_FILE, index=False)
        print(f"✓ 게놈 품질 평가 완료. 결과 저장: {QUALITY_FILE}")

        # 2. 병렬 fastANI 실행
        df_ani = run_fastani_parallel(fasta_files, max_workers=CPU_CORES)
        if df_ani.empty:
            print("경고: fastANI 결과가 없습니다. fastANI 경로 및 입력 파일을 확인하세요.")
            return
        df_ani.to_csv(ANI_FILE, sep='\t', index=False, header=False)

        # 3. 히트맵 생성
        print("\n[3/4] 시각화 (히트맵) 생성 시작...")
        all_genomes = sorted(pd.unique(df_ani[['Genome1', 'Genome2']].values.ravel()))
        ani_matrix = pd.DataFrame(np.identity(len(all_genomes)) * 100, index=all_genomes, columns=all_genomes)
        for _, row in df_ani.iterrows():
            ani_matrix.loc[row['Genome1'], row['Genome2']] = row['ANI']
            ani_matrix.loc[row['Genome2'], row['Genome1']] = row['ANI']
        
        # 전체 히트맵
        heatmap_full_path = os.path.join(OUTPUT_DIR, 'heatmap_full_ani.png')
        save_heatmap(ani_matrix, heatmap_full_path, "전체 게놈 ANI 유사도 히트맵")
        print(f"  - 전체 히트맵 저장 완료: {heatmap_full_path}")
        
        # 구간별 히트맵
        bins = np.arange(98, 100.1, 0.5)
        labels = [f"{bins[i]:.1f}-{bins[i+1]:.1f}" for i in range(len(bins)-1)]
        df_ani['ANI_Group'] = pd.cut(df_ani['ANI'], bins=bins, labels=labels, include_lowest=True, right=False)
        
        for grp_label in labels:
            subset = df_ani[df_ani['ANI_Group'] == grp_label]
            if subset.empty: continue
            sub_genomes = sorted(pd.unique(subset[['Genome1', 'Genome2']].values.ravel()))
            sub_matrix = ani_matrix.loc[sub_genomes, sub_genomes]
            fname = os.path.join(OUTPUT_DIR, f'heatmap_ani_{grp_label}.png')
            save_heatmap(sub_matrix, fname, f'ANI {grp_label}% 구간 히트맵')
        print("✓ 구간별 히트맵 생성 완료")

        # 4. 클러스터링 및 요약
        print("\n[4/4] ANI 임계값별 클러스터링 및 요약 시작...")
        quality_scores = dict(zip(df_quality['Genome'], df_quality['Quality_Score']))
        
        thresholds = sorted(set(np.round(np.arange(98.0, 100.01, 0.5), 2)) | set(np.round(np.arange(99.0, 100.01, 0.1), 2)))
        summary_rows = []
        for th in thresholds:
            clusters = build_clusters(df_ani, th)
            reps = [get_representative(c, quality_scores) for c in clusters]
            sizes = [len(c) for c in clusters]
            summary_rows.append({
                "ANI_threshold": th, "Num_groups": len(clusters),
                "Rep_samples": ";".join(reps), "Group_sizes": ";".join(map(str, sizes))
            })
            print(f"  - ANI {th}%: {len(clusters)}개 그룹 생성")

        # 결과 저장
        df_summary = pd.DataFrame(summary_rows)
        summary_csv_path = os.path.join(OUTPUT_DIR, "ANI_thresholds_clusters.csv")
        summary_txt_path = os.path.join(OUTPUT_DIR, "ANI_thresholds_clusters.txt")
        df_summary.to_csv(summary_csv_path, index=False)
        with open(summary_txt_path, "w", encoding="utf-8") as f:
            f.write("ANI 임계값별 클러스터링 요약\n")
            f.write("="*30 + "\n")
            for row in summary_rows:
                reps_preview = ", ".join(row['Rep_samples'].split(';')[:5])
                f.write(f"ANI 임계값: {row['ANI_threshold']}%\n")
                f.write(f"  - 그룹 수: {row['Num_groups']}\n")
                f.write(f"  - 대표 샘플 (최대 5개): {reps_preview}...\n")
                f.write(f"  - 각 그룹 크기: {row['Group_sizes']}\n\n")
        print(f"✓ 클러스터링 요약 완료. 결과 저장: {summary_csv_path}, {summary_txt_path}")
        
        print("\n\n>>> 모든 분석 파이프라인이 성공적으로 완료되었습니다. 'output' 폴더를 확인하세요. <<<")

    except FileNotFoundError:
        print(f"\n치명적 오류: fastANI 실행 파일({FASTANI_PATH})을 찾을 수 없습니다.")
        print("스크립트 상단의 'FASTANI_PATH' 변수에 올바른 절대 경로를 입력했는지 확인해주세요.")
    except Exception as e:
        print(f"\n예상치 못한 오류가 발생했습니다: {str(e)}")
        print("파일 권한, 패키지 설치 여부, 메모리 등을 확인해주세요.")

if __name__ == "__main__":
    main()