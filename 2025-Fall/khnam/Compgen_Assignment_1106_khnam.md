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
# Or: cd /mnt/c/Users/minion/Documents/Xiaoyue/Desktop/khnam_genomics/output_drep_results
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
Convert file name **mlst_results.tsv** to **Sepi_MLST.tsv**

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
