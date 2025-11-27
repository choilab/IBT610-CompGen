"""
L. plantarum 범유전체 분석 통합 파이프라인
- Roary 기반 Pan-genome 분석
- CompareM AAI 분석
- Surfaceome, COG, BGC 분석
- 12600K + RTX 4060 + 32GB RAM 최적화
"""

import os
import sys
import shutil
import subprocess
import time
import json
from pathlib import Path
from datetime import datetime

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from Bio import SeqIO
import psutil
import warnings
warnings.filterwarnings('ignore')

# ==================== 경로 설정 ====================
BASE_DIR = Path("/mnt/c/Users/AN/Desktop/gDrive/Study/Leacture/IBT610")
GFF_DIR = BASE_DIR / "gff"
FAA_DIR = BASE_DIR / "Prokka_faa/matched_prokka_faa"
ROARY_DIR = BASE_DIR / "Roary"
FNA_DIR = BASE_DIR / "Result_sample"
GBK_DIR = BASE_DIR / "Prokka_gbk"
OUTPUT_DIR = BASE_DIR / "output_analysis"
OUTPUT_DIR.mkdir(exist_ok=True)

# Galaxy eggNOG 결과 (선택적)
COG_ORTHOLOGS_FILE = ROARY_DIR / "core_gene_orthologs.tsv"
COG_ANNOTATIONS_FILE = ROARY_DIR / "core_gene_annotations.tsv"

# 시스템 설정
CPU_THREADS = 16
MEMORY_LIMIT = 32
GROUP_SIZE = 25

# CompareM 설정
COMPAREM_AAI_TIMEOUT = 1800  # 30분

# PyTorch GPU 감지
try:
    import torch
    GPU_AVAILABLE = torch.cuda.is_available()
    GPU_NAME = torch.cuda.get_device_name(0) if GPU_AVAILABLE else "Not detected"
except ImportError:
    GPU_AVAILABLE = False
    GPU_NAME = "PyTorch not installed"

print("="*80)
print("  L. plantarum 범유전체 분석 (12600K + RTX 4060)")
print("="*80)
print(f"✓ CPU: Intel 12600K ({CPU_THREADS}스레드)")
print(f"✓ GPU: {GPU_NAME}")
print(f"✓ RAM: 32GB (제한: {MEMORY_LIMIT}GB)")
print(f"✓ 시작: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# ==================== 공통 함수 ====================
def check_memory():
    """메모리 부족 경고"""
    mem = psutil.virtual_memory()
    if mem.available / (1024**3) < 6:
        print("⚠️ 메모리 부족! 실행을 중단합니다.")
        sys.exit(1)

def checkpoint_exists(step_name):
    """체크포인트 확인"""
    return (OUTPUT_DIR / f"checkpoint_{step_name}.flag").exists()

def mark_checkpoint(step_name):
    """체크포인트 생성"""
    with open(OUTPUT_DIR / f"checkpoint_{step_name}.flag", 'w') as f:
        f.write(f"Completed: {datetime.now()}\n")

# ==================== PART 1: Roary 범유전체 분석 ====================
def create_gff_mapping(gff_dir):
    """GFF3에서 locus_tag → gene_name 매핑"""
    print("\n[1/8] GFF3 매핑 생성 (Surfaceome용)...")
    check_memory()
    
    if checkpoint_exists("gff_mapping"):
        print("✓ 이미 완료됨 (체크포인트)")
        return {}
    
    if not gff_dir.exists():
        print(f"⚠️ GFF 디렉터리 없음: {gff_dir}")
        return {}
    
    mapping = {}
    gff_files = list(gff_dir.glob("*.gff3")) + list(gff_dir.glob("*.gff"))
    
    if not gff_files:
        print("⚠️ GFF3 파일 없음")
        return {}
    
    for idx, gff_file in enumerate(gff_files, 1):
        if idx % 50 == 0:
            print(f"  처리 중... {idx}/{len(gff_files)}")
        
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#') or '\t' not in line:
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                attrs = parts[8]
                import re
                lt_match = re.search(r'locus_tag=([^;]+)', attrs)
                gene_match = re.search(r'gene=([^;]+)', attrs)
                
                if lt_match and gene_match:
                    mapping[lt_match.group(1)] = gene_match.group(1)
    
    print(f"✓ 매핑 완료: {len(mapping):,}개")
    mark_checkpoint("gff_mapping")
    return mapping

def load_roary_data(roary_dir, output_dir):
    """Roary CSV 로드 및 CAR 분류"""
    print("\n[2/8] Roary 데이터 로드 (CAR 분류)...")
    check_memory()
    
    if checkpoint_exists("roary_load"):
        csv_file = output_dir / "roary_classified.csv"
        if csv_file.exists():
            print("✓ 이미 완료됨 (체크포인트)")
            return pd.read_csv(csv_file), None
    
    csv_file = roary_dir / "Galaxy3670-[Roary on data 3663, data 3662, and others Gene Presence Absence].csv"
    
    if not csv_file.exists():
        print(f"✗ Roary CSV 없음: {csv_file}")
        return None, None
    
    df = pd.read_csv(csv_file, low_memory=False)
    meta_cols = 14
    strain_cols = df.columns[meta_cols:]
    
    df['num_isolates'] = df[strain_cols].notna().sum(axis=1)
    df['frequency'] = (df['num_isolates'] / len(strain_cols)) * 100
    
    def categorize_gene(freq):
        if freq >= 99: return 'Core'
        elif freq >= 15: return 'Accessory'
        else: return 'Rare'
    
    df['category'] = df['frequency'].apply(categorize_gene)
    
    # CAR 분포 시각화
    car_counts = df['category'].value_counts()
    plt.figure(figsize=(8, 6))
    car_counts.plot(kind='bar', color=['green', 'orange', 'blue'])
    plt.title('CAR Classification')
    plt.xlabel('Category')
    plt.ylabel('Number of Genes')
    plt.tight_layout()
    plt.savefig(output_dir / '01_car_classification.png', dpi=300)
    plt.close()
    
    df.to_csv(output_dir / "roary_classified.csv", index=False)
    
    print(f"✓ 전체 유전자: {len(df):,}개 | 균주: {len(strain_cols)}개")
    print(f"✓ Core: {car_counts.get('Core', 0):,}개")
    print(f"✓ Accessory: {car_counts.get('Accessory', 0):,}개")
    print(f"✓ Rare: {car_counts.get('Rare', 0):,}개")
    
    mark_checkpoint("roary_load")
    return df, strain_cols

def analyze_heaps_law(df, strain_cols, output_dir):
    """Heaps' Law 계산"""
    print("\n[3/8] Heaps' Law 분석...")
    check_memory()
    
    if checkpoint_exists("heaps_law"):
        print("✓ 이미 완료됨 (체크포인트)")
        return None
    
    if len(strain_cols) < 5:
        print("⚠️ 균주 수 부족")
        return None
    
    sample_size = min(30, len(strain_cols))
    pan_sizes, core_sizes = [], []
    sample_range = range(5, sample_size + 1, 5)
    
    for n in sample_range:
        temp_pan, temp_core = [], []
        for _ in range(10):
            subset_cols = np.random.choice(strain_cols, n, replace=False)
            subset = df[subset_cols]
            temp_pan.append(subset.notna().any(axis=1).sum())
            temp_core.append(subset.notna().all(axis=1).sum())
        
        pan_sizes.append(np.mean(temp_pan))
        core_sizes.append(np.mean(temp_core))
    
    log_n = np.log(sample_range)
    log_pan = np.log(pan_sizes)
    coeffs = np.polyfit(log_n, log_pan, 1)
    lambda_val = coeffs[0]
    
    plt.figure(figsize=(10, 5))
    plt.plot(sample_range, pan_sizes, 'o-', label='Pan-genome', linewidth=2)
    plt.plot(sample_range, core_sizes, 'o-', label='Core-genome', linewidth=2)
    plt.xlabel('Number of Genomes')
    plt.ylabel('Number of Genes')
    plt.title(f'Heaps Law (λ = {lambda_val:.3f})')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.savefig(output_dir / '02_heaps_law.png', dpi=300)
    plt.close()
    
    print(f"✓ Heaps' Law λ = {lambda_val:.3f}")
    mark_checkpoint("heaps_law")
    return lambda_val

def analyze_mash(fna_dir, output_dir):
    """Mash clustering"""
    print("\n[4/8] Mash Clustering...")
    check_memory()
    
    if checkpoint_exists("mash"):
        print("✓ 이미 완료됨 (체크포인트)")
        return None
    
    if subprocess.run(["which", "mash"], capture_output=True).returncode != 0:
        print("⚠️ Mash 미설치")
        return None
    
    fna_files = list(fna_dir.glob("*.fna")) + list(fna_dir.glob("*.fasta"))
    if len(fna_files) < 2:
        print("⚠️ FASTA 파일 부족")
        return None
    
    sketch_file = output_dir / "mash_sketch.msh"
    cmd = f"mash sketch -o {sketch_file.with_suffix('')} -p {CPU_THREADS} -s 10000 " + " ".join([str(f) for f in fna_files[:50]])
    
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True)
        
        dist_file = output_dir / "mash_distances.tsv"
        cmd = f"mash dist {sketch_file} {sketch_file} > {dist_file}"
        subprocess.run(cmd, shell=True, check=True)
        
        # 시각화
        distances = []
        with open(dist_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    distances.append({
                        'genome1': Path(parts[0]).stem,
                        'genome2': Path(parts[1]).stem,
                        'distance': float(parts[2])
                    })
        
        if distances:
            dist_df = pd.DataFrame(distances)
            pivot_df = dist_df.pivot(index='genome1', columns='genome2', values='distance')
            dist_matrix = squareform(pivot_df.fillna(0).values)
            
            Z = linkage(dist_matrix, method='ward')
            
            plt.figure(figsize=(12, 8))
            dendrogram(Z, labels=pivot_df.index, leaf_rotation=90)
            plt.title('Mash-based Genome Clustering')
            plt.tight_layout()
            plt.savefig(output_dir / '03_mash_clustering.png', dpi=300)
            plt.close()
            
            print("✓ Mash clustering 완료")
            mark_checkpoint("mash")
            return dist_df
    except Exception as e:
        print(f"⚠️ Mash 실패: {e}")
    
    return None

def analyze_cog(roary_dir, output_dir):
    """Galaxy eggNOG 결과 활용 COG 분석"""
    print("\n[5/8] COG 분석 (Galaxy 결과 활용)...")
    check_memory()
    
    if checkpoint_exists("cog"):
        print("✓ 이미 완료됨 (체크포인트)")
        return None
    
    if not COG_ORTHOLOGS_FILE.exists():
        print(f"⚠️ COG orthologs 파일 없음: {COG_ORTHOLOGS_FILE}")
        print("    → COG 분석 생략 (Galaxy에서 완료 후 다시 실행)")
        return None
    
    try:
        ortho_df = pd.read_csv(COG_ORTHOLOGS_FILE, sep='\t')
        print(f"  → Orthologs: {len(ortho_df)}개")
        
        cog_col = None
        for col in ortho_df.columns:
            if 'cog' in col.lower() or 'category' in col.lower():
                cog_col = col
                break
        
        if cog_col is None:
            cog_col = ortho_df.columns[1]
        
        cog_counts = ortho_df[cog_col].value_counts().head(15)
        
        plt.figure(figsize=(14, 7))
        cog_counts.plot(kind='bar', color='steelblue')
        plt.title('COG Category Distribution')
        plt.xlabel('COG Category')
        plt.ylabel('Number of Genes')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(output_dir / '04_cog_analysis.png', dpi=300)
        plt.close()
        
        print(f"✓ COG 분석 완료: {len(cog_counts)}개 카테고리")
        mark_checkpoint("cog")
        return ortho_df
        
    except Exception as e:
        print(f"⚠️ COG 파일 파싱 실패: {e}")
        return None

def analyze_surfaceome(mapping, df, strain_cols, faa_dir, output_dir):
    """Surfaceome 분석"""
    print("\n[6/8] Surfaceome 분석...")
    check_memory()
    
    if checkpoint_exists("surfaceome"):
        csv_file = output_dir / 'surfaceome_predictions.csv'
        if csv_file.exists():
            print("✓ 이미 완료됨 (체크포인트)")
            return pd.read_csv(csv_file)
    
    if not mapping:
        print("⚠️ 매핑 테이블 없음")
        return None
    
    gene_map = dict(zip(df['Gene'], df['category']))
    faa_files = list(faa_dir.glob("*.faa")) + list(faa_dir.glob("*.fasta"))
    
    if not faa_files:
        print("⚠️ FAA 파일 없음")
        return None
    
    print(f"  {len(faa_files)}개 파일 분석 중...")
    
    class SurfaceomeAnalyzer:
        @staticmethod
        def analyze(seq):
            seq_str = str(seq).upper()
            if len(seq_str) < 30:
                return 'Intracellular'
            
            n_term = seq_str[:30]
            has_sp = sum(1 for aa in n_term if aa in 'LVIAFWM') >= 20
            
            c_term = seq_str[-60:]
            has_lpxtg = 'LP' in c_term and 'TG' in c_term
            
            hydrophobic = sum(1 for aa in seq_str if aa in 'LVIAFWM')
            is_tm = hydrophobic / len(seq_str) > 0.4
            
            if has_lpxtg and has_sp: return 'Cell_wall'
            elif has_sp: return 'Secreted'
            elif is_tm: return 'Membrane'
            else: return 'Intracellular'
    
    results = []
    for idx, faa_file in enumerate(faa_files, 1):
        if idx % 30 == 0:
            print(f"  처리 중... {idx}/{len(faa_files)}")
            check_memory()
        
        strain_name = faa_file.stem
        
        for record in SeqIO.parse(faa_file, "fasta"):
            locus_tag = record.id.split()[0]
            gene_name = mapping.get(locus_tag)
            
            if gene_name and gene_name in gene_map:
                loc = SurfaceomeAnalyzer.analyze(record.seq)
                results.append({
                    'Protein_ID': locus_tag,
                    'Strain': strain_name,
                    'Gene': gene_name,
                    'Category': gene_map[gene_name],
                    'Localization': loc
                })
    
    if not results:
        return None
    
    res_df = pd.DataFrame(results)
    res_df.to_csv(output_dir / 'surfaceome_predictions.csv', index=False)
    
    surface_df = res_df[res_df['Localization'] != 'Intracellular']
    
    if not surface_df.empty:
        ct = pd.crosstab(surface_df['Category'], surface_df['Localization'])
        ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100
        
        plt.figure(figsize=(10, 6))
        ct_pct.plot(kind='bar', stacked=True, colormap='viridis', ax=plt.gca())
        plt.title('Surfaceome Distribution by CAR Category')
        plt.xlabel('Gene Category')
        plt.ylabel('Percentage')
        plt.legend(title='Localization')
        plt.tight_layout()
        plt.savefig(output_dir / '05_surfaceome.png', dpi=300)
        plt.close()
        
        print(f"✓ 표면 단백질: {len(surface_df):,}개")
        mark_checkpoint("surfaceome")
        return res_df
    
    return None

# ==================== PART 2: CompareM AAI 분석 ====================
def find_comparem():
    """CompareM 도구 탐지"""
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if conda_prefix:
        possible_path = os.path.join(conda_prefix, 'bin', 'comparem')
        if os.path.exists(possible_path):
            return possible_path
    
    comparem_path = shutil.which('comparem')
    if comparem_path:
        return comparem_path
    
    print("⚠️ CompareM을 찾을 수 없습니다")
    print("   설치: conda install -c bioconda comparem")
    return None

def run_comparem_aai(faa_dir, output_dir):
    """CompareM AAI 분석"""
    print("\n[7/8] CompareM AAI 분석...")
    check_memory()
    
    if checkpoint_exists("comparem_aai"):
        print("✓ 이미 완료됨 (체크포인트)")
        return True
    
    comparem_path = find_comparem()
    if not comparem_path:
        return False
    
    faa_files = list(faa_dir.glob("*.fasta"))
    if not faa_files:
        print("⚠️ FAA 파일 없음")
        return False
    
    # 그룹 분할
    groups = [faa_files[i:i+GROUP_SIZE] 
             for i in range(0, len(faa_files), GROUP_SIZE)]
    
    print(f"  {len(faa_files)}개 파일 → {len(groups)}개 그룹")
    
    all_results = []
    
    for gi, files in enumerate(groups, 1):
        group_dir = output_dir / f"aai_group_{gi}"
        group_dir.mkdir(exist_ok=True)
        
        # 파일 복사
        for fname in files:
            shutil.copy2(fname, group_dir / fname.name)
        
        aai_out = group_dir / "output_aai"
        aai_out.mkdir(exist_ok=True)
        
        print(f"\n  [그룹 {gi}/{len(groups)}] 실행 중...")
        
        cmd = [
            comparem_path, "aai_wf",
            "--proteins", "--file_ext", "fasta",
            "--cpus", str(CPU_THREADS),
            str(group_dir), str(aai_out)
        ]
        
        try:
            subprocess.run(cmd, timeout=COMPAREM_AAI_TIMEOUT, check=True)
            
            result_file = aai_out / "aai" / "aai_summary.tsv"
            if result_file.exists():
                df = pd.read_csv(result_file, sep="\t")
                all_results.append(df)
                print(f"    ✓ 성공")
        except subprocess.TimeoutExpired:
            print(f"    ✗ 타임아웃")
        except Exception as e:
            print(f"    ✗ 오류: {e}")
    
    if all_results:
        merged_df = pd.concat(all_results, ignore_index=True)
        merged_df.to_csv(output_dir / "aai_summary_merged.tsv", sep='\t', index=False)
        print(f"\n✓ AAI 분석 완료: {len(merged_df)}개 비교")
        mark_checkpoint("comparem_aai")
        return True
    
    return False

# ==================== 메인 실행 ====================
def main():
    start_time = time.time()
    
    print("\n시작 분석...")
    print(f"메모리: {psutil.virtual_memory().available / (1024**3):.1f}GB 사용 가능")
    
    # 1. GFF 매핑
    mapping = create_gff_mapping(GFF_DIR)
    
    # 2. Roary 데이터
    roary_df, strain_cols = load_roary_data(ROARY_DIR, OUTPUT_DIR)
    if roary_df is None:
        sys.exit(1)
    
    # 3. Heaps' Law
    lambda_val = analyze_heaps_law(roary_df, strain_cols, OUTPUT_DIR)
    
    # 4. Mash
    mash_df = analyze_mash(FNA_DIR, OUTPUT_DIR)
    
    # 5. COG
    cog_df = analyze_cog(ROARY_DIR, OUTPUT_DIR)
    
    # 6. Surfaceome
    surface_df = analyze_surfaceome(mapping, roary_df, strain_cols, FAA_DIR, OUTPUT_DIR)
    
    # 7. CompareM AAI
    aai_success = run_comparem_aai(FAA_DIR, OUTPUT_DIR)
    
    # 결과 요약
    elapsed = time.time() - start_time
    print("\n" + "="*80)
    print("  분석 완료 요약")
    print("="*80)
    print(f"✓ Pangenome: {len(roary_df):,} 유전자")
    print(f"✓ Heaps' Law λ: {lambda_val:.3f}" if lambda_val else "✗ Heaps' Law: 생략")
    print(f"✓ Mash: {'완료' if mash_df is not None else '생략'}")
    print(f"✓ COG: {'완료' if cog_df is not None else '생략'}")
    print(f"✓ Surfaceome: {len(surface_df) if surface_df is not None else 0:,} 단백질")
    print(f"✓ AAI: {'완료' if aai_success else '생략'}")
    print(f"\n⏱️ 실행 시간: {elapsed/60:.1f}분")
    print(f"📁 결과: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
