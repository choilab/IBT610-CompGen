# MLST Analysis of 124 *S. epidermidis* Genomes

**Project Goal:** To perform Multi-Locus Sequence Typing (MLST) on the 124 representative *S. epidermidis* genomes previously selected by ANI filtering. This analysis will assign a Sequence Type (ST) to each genome, allowing for the characterization of their genetic lineages.

---

## 1. Conda Environment Setup

Create and activate a new Conda virtual environment named `mlst_env` to install the `mlst` tool.

```bash
# 1. Create the Conda environment
# (If the environment already exists, skip this step and proceed to activation)
conda create -n mlst_env -c bioconda mlst

# 2. Activate the Conda environment
conda activate mlst_env
```

## 2. Running the MLST Analysis

After activating the environment, navigate to the project directory containing your 124 genomes (e.g., `Rep_genomes/` or `dereplicated_genomes/`).

```bash
# Example: cd /mnt/c/Compgen_project
# Or: cd /mnt/c/Users/khnam_genomics/output_drep_results
```

From this location, execute the `mlst` command. This will analyze all genomes and save the results to `mlst_results.tsv`.

* **`--scheme sepidermidis`**: Explicitly specifies the standard MLST scheme for *S. epidermidis*.
* **`--legacy`**: Uses the traditional PubMLST.org definitions.
* **`--nopath`**: Omits the full file path from the results, showing only the filename.
* **`Rep_genomes/*.fna`**: Input files (adjust the folder name and extension, e.g., `*.fasta`, if needed).

```bash
# Run MLST on all 124 genomes
mlst --legacy --scheme sepidermidis --nopath Rep_genomes/*.fna > mlst_results.tsv
```

## 3. Reviewing and Interpreting Results

Check the top of the output file using `head`.

```bash
head mlst_results.tsv
```
Rename **mlst_results.tsv** to **Sepi_MLST.tsv** (File name)

### Output File Structure (`Sepi_MLST.tsv`)

The result is a Tab-Separated Value (TSV) table, structured as follows:

| FILE | SCHEME | ST | arcC | aroE | gtr | mutS | ... |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| GCF\_001.fna | sepidermidis | 184 | 16 | 1 | 2 | 1 | ... |
| GCF\_002.fna | sepidermidis | 1016 | 1 | 13 | 7 | 6 | ... |
| GCF\_003.fna | sepidermidis | 20 | 1 | 1 | 2 | 2 | ... |
| GCF\_004.fna | sepidermidis | 20 | 1 | 1 | 2 | 2 | ... |
| GCF\_005.fna | sepidermidis | - | 1 | 1 | 2 | - | ... |

* **FILE**: The input genome filename.
* **SCHEME**: The scheme used (`sepidermidis`).
* **ST**: **The most important result.** This is the **Sequence Type (ST) number**, or "genetic barcode," assigned to the strain.
* `arcC`, `aroE`...: The specific allele numbers identified for each of the 7 housekeeping genes.
* **`-` (hyphen)**: Represents a missing value. This means the gene was not found (e.g., deletion, mutation) or was truncated, preventing ST assignment.

## 4. Summarizing ST Frequencies

To understand the distribution of STs among the 124 genomes, open `Sepi_MLST.tsv` in **Excel or Google Sheets**.

1.  **Sort the data** by the **`ST` column** (Sort A-Z). This will group all genomes with the same ST together.
2.  Use the **"Pivot Table"** feature.
    * Set **`ST`** as the 'Row'.
    * Set **`ST`** as the 'Value' and ensure it is summarizing by **'Count'** (not 'Sum').

### The number of ST  (from Pivot Table)

This table summarizes how many genomes belong to each ST.

| ST | Count of ST |
| :--- | :--- |
| 2 | 18 |
| 5 | 5 |
| 20 | 3 |
| 6 | 3 |
| - | 13 |
| ... | ... |
| **Total** | **124** |

This summary shows that the **ST-2 lineage is the most prevalent** (18 genomes) in this dataset, and 13 genomes were "untypeable" (`-`).
Overall, **69 unique STs** (including the '-' untypeable group) were identified across the 124 genomes, indicating significant genetic diversity within this collection.


---

## AAI Analysis and Heatmap (124 Genomes)

**Goal:** To calculate the all-versus-all Average Amino-acid Identity (AAI) for the 124 representative genomes and visualize the similarity matrix as a heatmap.

### 1. Conda Environment Setup (`comparem`)

A dedicated Conda environment is created to install the `comparem` tool.

```bash
# 1. Create the Conda environment
# (Skip if already created)
conda create -n comparem_env comparem -c bioconda

# 2. Activate the Conda environment
conda activate comparem_env
```

### 2. Running the CompareM AAI Workflow


**2-1. Navigate to the project directory:**
The directory must contain the folder with the 124 genome files (e.g., `dereplicated_genomes/`).

```bash
# Path on the analysis machine (WSL)
cd /mnt/c/Users/khnam/khnam_genomics/output_drep_resultsj
```

**2-2. Create a temporary directory:**
`comparem` requires a pre-existing temporary directory.

```bash
mkdir comparem_temp
```

**2-3. Execute the AAI workflow using `nohup`:**
`nohup` ensures the process continues running even if the terminal is closed.

* **`dereplicated_genomes/`**: Input folder containing 124 genomes.
* **`AAI_results/`**: Output directory for results.
* **`-c 4`**: Number of CPU cores used.

```bash
nohup comparem aai_wf dereplicated_genomes/ AAI_results/ -c 4 --tmp_dir ./comparem_temp &
```

**2-4. Monitor progress (optional):**
```bash
tail -f nohup.out
```

**2-5. Key Output File:**
The analysis generates a summary file containing all pairwise AAI values:
`AAI_results/aai/aai_summary.tsv`

---

### 3. AAI Heatmap Visualization (Jupyter Notebook)

The `aai_summary.tsv` file (exported as `.xlsx` for convenience) was moved to a local notebook for visualization.

### 3-1. Environment Setup (Jupyter)

The following packages are required in the Conda environment to run the notebook.

```bash
# Ensure environment is active
conda activate comparem_env

# Install necessary libraries
conda install jupyter pandas seaborn matplotlib openpyxl
```

### 3-2. Running Jupyter Notebook

From the project directory (containing `aai_summary.xlsx`):
```bash
jupyter notebook
```

### 3-3. Python Code for Heatmap

The following code was run in a new Jupyter Notebook cell to load the data and generate the heatmap.

```python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Magic command to display plots inline in Jupyter
%matplotlib inline 

# --- 1. Load Data (from .xlsx file) ---
input_file = 'aai_summary.xlsx' 

print(f"Loading data from {input_file}...")
df = pd.read_excel(
    input_file,
    header=None, # No header in the original data rows
    skiprows=1,  # Skip the first row (which contained text headers)
    names=["Genome_A", "Genes_A", "Genome_B", "Genes_B", "Orthos", "AAI", "Std_Dev", "Orth_Frac"]
)

aai_df = df[['Genome_A', 'Genome_B', 'AAI']]
print("Data loaded successfully.")

# --- 2. Create the AAI Matrix ---
print("Pivoting data to create AAI matrix...")
aai_matrix = aai_df.pivot(
    index='Genome_A',
    columns='Genome_B',
    values='AAI'
)
# Make the matrix symmetrical
aai_matrix_symmetric = aai_matrix.combine_first(aai_matrix.T)
# Fill the diagonal (self-vs-self) with 100
aai_matrix_symmetric.fillna(100, inplace=True)
print("Matrix successfully created.")

# --- 3. Clean Genome Names (Optional) ---
clean_names = [
    str(name).replace('.fna', '').replace('.fasta', '').replace('.fa', '') 
    for name in aai_matrix_symmetric.columns
]
aai_matrix_symmetric.columns = clean_names
aai_matrix_symmetric.index = clean_names

# --- 4. Generate and Display Heatmap ---
print("Generating AAI heatmap...")

plt.figure(figsize=(18, 16)) 

sns.heatmap(
    aai_matrix_symmetric,
    cmap='viridis',     # Color map
    vmin=70,            # Minimum value for color bar
    vmax=100,           # Maximum value for color bar (fixed typo: 'vmax')
    xticklabels=False,  # Hide X-axis labels
    yticklabels=False   # Hide Y-axis labels
)

plt.title('AAI Heatmap (124 Representative Genomes)', fontsize=20)
plt.xlabel('Genomes')
plt.ylabel('Genomes')

# Save the plot (optional)
# plt.savefig('AAI_heatmap.png', dpi=300, bbox_inches='tight')

# Show the plot in the notebook
plt.show()
```

---

## 4. Final Heatmap Result

The resulting heatmap visualizes the pairwise AAI values for all 124 genomes.

<img width="1328" height="1306" alt="image" src="https://github.com/user-attachments/assets/8afff74f-8463-4fea-a6e6-8a3f7c8b5a88" />


### Interpretation

The heatmap is dominated by yellow and green colors, indicating extremely high similarity. As shown by the color bar (70-100), **almost all pairwise comparisons fall between 90% and 100% AAI.**

This demonstrates that while the genomes were filtered by 99% ANI (nucleotide) and show diverse MLST types (7 genes), their overall **proteome-level (AAI) similarity is remarkably high** and homogeneous.
