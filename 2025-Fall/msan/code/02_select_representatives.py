import pandas as pd
import networkx as nx
import os
import shutil

# --- 1. 경로 설정 ---
# 필요한 모든 경로를 변수로 지정하여 관리합니다.
SOURCE_GENOMES_DIR = "data/genomes"
DEST_REPRESENTATIVES_DIR = "results/representative_genomes"
INPUT_ANI_FILE = "results/skani_raw_output.txt"

def main():
    """워크플로우 전체를 실행하는 메인 함수"""

    # --- 2. 100% ANI 쌍 필터링 ---
    print("Step 1: SKANI 결과에서 100% ANI 쌍을 필터링합니다...")
    try:
        df = pd.read_csv(INPUT_ANI_FILE, sep="\t")
    except FileNotFoundError:
        print(f"오류: 입력 파일 '{INPUT_ANI_FILE}'을 찾을 수 없습니다.")
        print("먼저 01_prepare_and_run_skani.sh 스크립트를 실행해주세요.")
        return

    filtered_df = df[df['ANI'] >= 99.99]
    print(f"-> 필터링 완료. {len(filtered_df)}개의 100% ANI 쌍을 찾았습니다.")

    # --- 3. 네트워크 분석으로 대표 균주 선정 ---
    print("\nStep 2: 네트워크 분석으로 대표 균주를 선정합니다...")
    G = nx.from_pandas_edgelist(filtered_df, 'Ref_file', 'Query_file')
    clusters = list(nx.connected_components(G))
    representatives = [sorted(list(cluster))[0] for cluster in clusters]
    print(f"-> 분석 완료. 총 {len(representatives)}개의 대표 균주를 선정했습니다.")

    # --- 4. 대표 균주 FASTA 파일 복사 ---
    print(f"\nStep 3: 대표 균주 파일을 '{DEST_REPRESENTATIVES_DIR}' 디렉터리로 복사합니다...")
    os.makedirs(DEST_REPRESENTATIVES_DIR, exist_ok=True)
    
    copied_count = 0
    for filename in representatives:
        source_path = os.path.join(SOURCE_GENOMES_DIR, filename)
        if os.path.exists(source_path):
            shutil.copy2(source_path, DEST_REPRESENTATIVES_DIR)
            copied_count += 1
        else:
            print(f"경고: 원본 파일 '{source_path}'를 찾을 수 없습니다.")
    
    print(f"-> 복사 완료. 총 {copied_count}개의 파일을 복사했습니다.")
    
    # --- 5. 최종 결과 확인 ---
    final_file_count = len(os.listdir(DEST_REPRESENTATIVES_DIR))
    print(f"\n최종 확인: '{DEST_REPRESENTATIVES_DIR}' 디렉터리에 {final_file_count}개의 파일이 있습니다.")

if __name__ == "__main__":
    main()
