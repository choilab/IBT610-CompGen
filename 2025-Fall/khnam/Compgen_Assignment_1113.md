## 1. 📊 Roary (Pan-Genome) Summary Statistics Analysis

The `summary_statistics.txt` file from the Roary analysis provides the high-level overview of the pan-genome structure for your 124-genome collection.

---

### 1. View Summary Statistics (Terminal Command)

In notebook's Ubuntu terminal, navigate to your Roary results folder and use `cat` to display the summary file.

```bash
# (base) user@LAPTOP:~$
# Navigate to the Roary results folder (adjust the name as needed)
cd /mnt/c/Compgen_project/Roary_results_1762490311
    
# Print the contents of the summary statistics file
cat summary_statistics.txt
```

| Category | Gene Count |
| :--- | :--- |
| **Core genes** (99% <= strains <= 100%) | **1522** |
| Soft core genes (95% <= strains < 99%)| 189 |
| Shell genes (15% <= strains < 95%) | 948 |
| Cloud genes (0% <= strains < 15%) | 7355 |
| **Total genes** (0% <= strains <= 100%) | **10014** |

## 2. Functional Genome Profile

* **Contents:** This is a **10,014 x 124 matrix** (a "gene presence/absence checklist") showing which of the 10,014 pan-genes (rows) are present or absent (O/X) in each of the 124 genomes (columns).
* **Usage:**
    1.  This checklist was used to compare the presence or absence of specific genes (like `mecA`) between groups (e.g., `ST-2` vs `ST-5`).
    2.  It served as the raw data for the clustering in section 3.3.

---

## 3. Gene-Profile-Based Clustering

* **Contents:** This is a **phylogenetic tree (dendrogram)** that groups **'genomes with similar gene compositions'** together, based on the.
* **Principle:** Genomes that share a more similar set of the 8,492 accessory genes (the O/X list) are clustered onto closer branches in this tree.

<img width="1920" height="912" alt="image" src="https://github.com/user-attachments/assets/ab00dcd2-c0da-4606-bbd3-10ed721b8a2f" />

