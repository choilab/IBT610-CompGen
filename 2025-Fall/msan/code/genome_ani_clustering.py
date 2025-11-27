"""
게놈 ANI 기반 클러스터링 통합 파이프라인
- fastANI/skANI를 이용한 게놈 유사도 분석
- 대표 균주 자동 선정
- 품질 평가 및 시각화
"""

import os
import sys
import glob
import shutil
import subprocess
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations
from collections import defaultdict

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
import networkx as nx

# ==================== 사용자 설정 ====================
# ANI 도구 선택: 'fastani' 또는 'skani'
ANI_TOOL = "fastani"  # "skani"로 변경 가능

# fastANI 경로 (fastANI 사용 시 필수)
FASTANI_PATH = "/root/bin/fastANI"

# 입력/출력 경로
SOURCE_GENOMES_DIR = "data/genomes"
OUTPUT_DIR = "results"
REPRESENTATIVES_DIR = os.path.join(OUTPUT_DIR, "representative_genomes")

# 분석 설정
CPU_CORES = os.cpu_count() or 1
ANI_THRESHOLD = 99.99  # 100% ANI 판정 기준

# ==================== 공통 함수 ====================
def setup_directories():
    """필요한 디렉터리 생성"""
    os.makedirs(SOURCE_GENOMES_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(REPRESENTATIVES_DIR, exist_ok=True)

def set_korean_font():
    """matplotlib 한글 폰트 설정"""
    try:
        import matplotlib.font_manager as fm
        font_paths = [
            "C:/Windows/Fonts/malgun.ttf",
            "/usr/share/fonts/truetype/nanum/NanumGothic.ttf",
            "/System/Library/Fonts/AppleGothic.ttf"
        ]
        for fp in font_paths:
            if os.path.exists(fp):
                plt.rcParams['font.family'] = fm.FontProperties(fname=fp).get_name()
                return
    except Exception as e:
        print(f"폰트 설정 오류: {e}")

def validate_tool():
    """ANI 도구 설치 확인"""
    if ANI_TOOL == "fastani":
        if not os.path.exists(FASTANI_PATH):
            print(f"오류: fastANI를 찾을 수 없습니다: {FASTANI_PATH}")
            sys.exit(1)
        return FASTANI_PATH
    elif ANI_TOOL == "skani":
        result = subprocess.run(["which", "skani"], capture_output=True)
        if result.returncode != 0:
            print("오류: skani가 설치되어 있지 않습니다")
            print("설치: conda install -c bioconda skani")
            sys.exit(1)
        return "skani"
    else:
        print(f"오류: 지원하지 않는 도구입니다: {ANI_TOOL}")
        sys.exit(1)

# ==================== 1단계: 게놈 품질 평가 ====================
def calculate_assembly_stats(fasta_file):
    """FASTA 파일의 조립 통계 계산"""
    lengths, gc_count, total_len = [], 0, 0
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        l = len(seq)
        lengths.append(l)
        gc_count += seq.count('G') + seq.count('C')
        total_len += l
    
    if not lengths:
        return {"N50": 0, "L50": 0, "Total_Length": 0, 
                "Seq_Count": 0, "GC_Content": 0}
    
    lengths.sort(reverse=True)
    cum_sum = np.cumsum(lengths)
    half_total = total_len / 2
    
    n50_val = next((l for i, l in enumerate(lengths) 
                    if cum_sum[i] >= half_total), 0)
    l50_val = next((i + 1 for i, l in enumerate(lengths) 
                    if cum_sum[i] >= half_total), 0)
    gc_content = (gc_count / total_len * 100) if total_len else 0
    
    return {
        "N50": n50_val,
        "L50": l50_val,
        "Total_Length": total_len,
        "Seq_Count": len(lengths),
        "GC_Content": round(gc_content, 2)
    }

def assess_genome_quality(fasta_files):
    """모든 게놈의 품질 평가"""
    print("\n[1/5] 게놈 품질 평가 중...")
    stats_list = []
    
    for i, f in enumerate(fasta_files, 1):
        genome_name = os.path.splitext(os.path.basename(f))[0]
        stats = calculate_assembly_stats(f)
        stats['Genome'] = genome_name
        stats['File'] = os.path.basename(f)
        stats_list.append(stats)
        print(f"  ({i}/{len(fasta_files)}) {genome_name} | N50={stats['N50']:,}")
    
    df_quality = pd.DataFrame(stats_list)
    
    # 품질 점수 계산 (Z-score 기반)
    for col in ['N50', 'Total_Length']:
        if df_quality[col].std() > 0:
            df_quality[f'{col}_Z'] = (
                (df_quality[col] - df_quality[col].mean()) / df_quality[col].std()
            )
        else:
            df_quality[f'{col}_Z'] = 0
    
    df_quality['Quality_Score'] = (
        df_quality['N50_Z'] + 
        0.5 * df_quality['Total_Length_Z'] - 
        0.2 * df_quality['L50']
    ).rank(ascending=False)
    
    quality_file = os.path.join(OUTPUT_DIR, "genome_quality.csv")
    df_quality.to_csv(quality_file, index=False)
    print(f"✓ 품질 평가 완료: {quality_file}")
    
    return df_quality

# ==================== 2단계: ANI 분석 ====================
def run_fastani_parallel(fasta_files, max_workers=CPU_CORES):
    """병렬 fastANI 실행"""
    pairs = list(combinations(fasta_files, 2))
    total_pairs = len(pairs)
    results = []
    
    print(f"\n[2/5] fastANI 분석 중 ({total_pairs}쌍, {max_workers}코어)")
    start_time = time.time()
    
    def fastani_worker(query, ref):
        q_name = os.path.splitext(os.path.basename(query))[0]
        r_name = os.path.splitext(os.path.basename(ref))[0]
        tmp_out = os.path.join(OUTPUT_DIR, f"tmp_{q_name}_vs_{r_name}.txt")
        
        cmd = [
            FASTANI_PATH, "--query", query, "--ref", ref,
            "--threads", "1", "--output", tmp_out
        ]
        
        try:
            subprocess.run(cmd, check=True, timeout=600, 
                          capture_output=True, text=True)
            if os.path.exists(tmp_out):
                with open(tmp_out, 'r') as f:
                    line = f.readline()
                    if line:
                        parts = line.strip().split("\t")
                        if len(parts) >= 3:
                            return [q_name, r_name, float(parts[2])]
        except Exception as e:
            print(f"\n오류: {q_name} vs {r_name} - {str(e)}")
        finally:
            if os.path.exists(tmp_out):
                os.remove(tmp_out)
        return None
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(fastani_worker, q, r) 
                  for q, r in pairs}
        for i, future in enumerate(as_completed(futures), 1):
            result = future.result()
            if result:
                results.append(result)
            
            elapsed = time.time() - start_time
            pct = int(100 * i / total_pairs)
            print(f"[{pct:3d}%] ({i}/{total_pairs}) | "
                  f"경과: {int(elapsed//60)}분 {int(elapsed%60)}초", 
                  end="\r", flush=True)
    
    print("\n✓ fastANI 분석 완료")
    return pd.DataFrame(results, columns=['Genome1', 'Genome2', 'ANI'])

def run_skani(fasta_files):
    """skANI 실행"""
    print(f"\n[2/5] skANI 분석 중...")
    
    # 파일 목록 생성
    list_file = os.path.join(OUTPUT_DIR, "genome_list.txt")
    with open(list_file, 'w') as f:
        for fasta in fasta_files:
            f.write(f"{fasta}\n")
    
    # skANI 실행
    output_file = os.path.join(OUTPUT_DIR, "skani_raw_output.txt")
    cmd = [
        "skani", "triangle",
        "-l", list_file,
        "-o", output_file,
        "-t", str(CPU_CORES)
    ]
    
    try:
        subprocess.run(cmd, check=True)
        df = pd.read_csv(output_file, sep="\t")
        df.columns = ['Ref_file', 'Query_file', 'ANI', 'Aligned_fraction_ref', 
                     'Aligned_fraction_query', 'Ref_name', 'Query_name']
        
        # 파일명만 추출
        df['Genome1'] = df['Ref_file'].apply(
            lambda x: os.path.splitext(os.path.basename(x))[0]
        )
        df['Genome2'] = df['Query_file'].apply(
            lambda x: os.path.splitext(os.path.basename(x))[0]
        )
        
        print("✓ skANI 분석 완료")
        return df[['Genome1', 'Genome2', 'ANI']]
    except Exception as e:
        print(f"skANI 실행 오류: {e}")
        sys.exit(1)

# ==================== 3단계: 시각화 ====================
def create_visualizations(df_ani):
    """ANI 결과 시각화"""
    print("\n[3/5] 시각화 생성 중...")
    
    # ANI 매트릭스 생성
    all_genomes = sorted(pd.unique(
        df_ani[['Genome1', 'Genome2']].values.ravel()
    ))
    ani_matrix = pd.DataFrame(
        np.identity(len(all_genomes)) * 100,
        index=all_genomes,
        columns=all_genomes
    )
    
    for _, row in df_ani.iterrows():
        ani_matrix.loc[row['Genome1'], row['Genome2']] = row['ANI']
        ani_matrix.loc[row['Genome2'], row['Genome1']] = row['ANI']
    
    # 전체 히트맵
    plt.figure(figsize=(min(30, max(12, len(all_genomes)//2)), 
                       min(30, max(12, len(all_genomes)//2))))
    sns.heatmap(ani_matrix.astype(float), cmap='viridis', annot=False)
    plt.title("전체 게놈 ANI 유사도 히트맵", fontsize=16)
    plt.tight_layout()
    heatmap_path = os.path.join(OUTPUT_DIR, 'ani_heatmap_full.png')
    plt.savefig(heatmap_path, dpi=300)
    plt.close()
    print(f"  - 히트맵 저장: {heatmap_path}")
    
    # 구간별 히트맵
    bins = np.arange(98, 100.1, 0.5)
    labels = [f"{bins[i]:.1f}-{bins[i+1]:.1f}" 
             for i in range(len(bins)-1)]
    df_ani['ANI_Group'] = pd.cut(
        df_ani['ANI'], bins=bins, labels=labels,
        include_lowest=True, right=False
    )
    
    for grp_label in labels:
        subset = df_ani[df_ani['ANI_Group'] == grp_label]
        if subset.empty:
            continue
        
        sub_genomes = sorted(pd.unique(
            subset[['Genome1', 'Genome2']].values.ravel()
        ))
        sub_matrix = ani_matrix.loc[sub_genomes, sub_genomes]
        
        plt.figure(figsize=(max(12, len(sub_genomes)//2), 
                           max(10, len(sub_genomes)//2)))
        sns.heatmap(sub_matrix.astype(float), cmap='viridis', annot=False)
        plt.title(f'ANI {grp_label}% 구간 히트맵', fontsize=16)
        plt.tight_layout()
        fname = os.path.join(OUTPUT_DIR, f'ani_heatmap_{grp_label}.png')
        plt.savefig(fname, dpi=300)
        plt.close()
    
    print("✓ 시각화 완료")
    return ani_matrix

# ==================== 4단계: 클러스터링 ====================
def build_clusters(df_ani, threshold):
    """ANI 기반 그래프 클러스터링"""
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
            current = queue.pop(0)
            group.add(current)
            for neighbor in edges.get(current, set()):
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        
        if group:
            clusters.append(sorted(list(group)))
    
    return clusters

def select_representatives(df_quality, df_ani):
    """대표 균주 자동 선정"""
    print(f"\n[4/5] 대표 균주 선정 (ANI ≥ {ANI_THRESHOLD}%)...")
    
    # 100% ANI 쌍 필터링
    filtered_df = df_ani[df_ani['ANI'] >= ANI_THRESHOLD]
    print(f"  - {len(filtered_df)}개의 100% ANI 쌍 발견")
    
    if filtered_df.empty:
        print("  - 모든 게놈이 고유함 (대표 선정 불필요)")
        return list(df_quality['Genome'])
    
    # 네트워크 클러스터링
    G = nx.from_pandas_edgelist(filtered_df, 'Genome1', 'Genome2')
    clusters = list(nx.connected_components(G))
    
    # 품질 점수 딕셔너리
    quality_scores = dict(zip(df_quality['Genome'], 
                             df_quality['Quality_Score']))
    
    # 각 클러스터에서 최고 품질 게놈 선택
    representatives = []
    for cluster in clusters:
        rep = max(cluster, key=lambda g: quality_scores.get(g, 0))
        representatives.append(rep)
    
    print(f"✓ {len(representatives)}개 대표 균주 선정 완료")
    return representatives

# ==================== 5단계: 결과 저장 ====================
def save_results(df_ani, representatives, df_quality):
    """결과 파일 저장"""
    print("\n[5/5] 결과 저장 중...")
    
    # ANI 결과 저장
    ani_file = os.path.join(OUTPUT_DIR, "all_pairs_ani.tsv")
    df_ani.to_csv(ani_file, sep='\t', index=False)
    print(f"  - ANI 결과: {ani_file}")
    
    # 대표 균주 파일 복사
    copied = 0
    for genome in representatives:
        file_row = df_quality[df_quality['Genome'] == genome]
        if file_row.empty:
            continue
        
        source_file = file_row.iloc[0]['File']
        source_path = os.path.join(SOURCE_GENOMES_DIR, source_file)
        
        if os.path.exists(source_path):
            shutil.copy2(source_path, REPRESENTATIVES_DIR)
            copied += 1
    
    print(f"  - 대표 균주 복사: {copied}개 파일 → {REPRESENTATIVES_DIR}")
    
    # 요약 보고서
    summary_file = os.path.join(OUTPUT_DIR, "analysis_summary.txt")
    with open(summary_file, 'w', encoding='utf-8') as f:
        f.write("게놈 ANI 클러스터링 분석 요약\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"분석 도구: {ANI_TOOL.upper()}\n")
        f.write(f"전체 게놈 수: {len(df_quality)}\n")
        f.write(f"대표 균주 수: {len(representatives)}\n")
        f.write(f"ANI 임계값: {ANI_THRESHOLD}%\n\n")
        f.write("대표 균주 목록:\n")
        for rep in representatives:
            f.write(f"  - {rep}\n")
    
    print(f"  - 요약 보고서: {summary_file}")
    print("\n✓ 모든 분석 완료!")

# ==================== 메인 실행 ====================
def main():
    print("="*60)
    print("  게놈 ANI 기반 클러스터링 파이프라인")
    print("="*60)
    
    start_time = time.time()
    
    # 초기 설정
    setup_directories()
    set_korean_font()
    tool_path = validate_tool()
    
    # 입력 파일 확인
    fasta_files = (glob.glob(os.path.join(SOURCE_GENOMES_DIR, '*.fna')) +
                   glob.glob(os.path.join(SOURCE_GENOMES_DIR, '*.fasta')))
    
    if not fasta_files:
        print(f"오류: '{SOURCE_GENOMES_DIR}'에 게놈 파일이 없습니다")
        sys.exit(1)
    
    print(f"\n입력: {len(fasta_files)}개 게놈 파일")
    print(f"도구: {ANI_TOOL.upper()}")
    print(f"CPU: {CPU_CORES}코어\n")
    
    # 1. 품질 평가
    df_quality = assess_genome_quality(fasta_files)
    
    # 2. ANI 분석
    if ANI_TOOL == "fastani":
        df_ani = run_fastani_parallel(fasta_files)
    else:  # skani
        df_ani = run_skani(fasta_files)
    
    if df_ani.empty:
        print("오류: ANI 결과가 없습니다")
        sys.exit(1)
    
    # 3. 시각화
    ani_matrix = create_visualizations(df_ani)
    
    # 4. 대표 선정
    representatives = select_representatives(df_quality, df_ani)
    
    # 5. 결과 저장
    save_results(df_ani, representatives, df_quality)
    
    # 실행 시간
    elapsed = time.time() - start_time
    print(f"\n총 실행 시간: {int(elapsed//60)}분 {int(elapsed%60)}초")
    print(f"결과 폴더: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
