## 1. 📊 Roary (Pan-Genome) Summary Statistics Analysis

The `summary_statistics.txt` file from the Roary analysis provides the high-level overview of the pan-genome structure for your 124-genome collection.

---

### 1. View Summary Statistics (Terminal Command)

In your notebook's Ubuntu terminal, navigate to your Roary results folder and use `cat` to display the summary file.

```bash
# (base) user@LAPTOP:~$
# Navigate to the Roary results folder (adjust the name as needed)
cd /mnt/c/Compgen_project/Roary_results_1762490311
    
# Print the contents of the summary statistics file
cat summary_statistics.txt

| Category | Gene Count |
| :--- | :--- |
| **Core genes** (99% <= strains <= 100%) | **1522** |
| Soft core genes (95% <= strains < 99%)| 189 |
| Shell genes (15% <= strains < 95%) | 948 |
| Cloud genes (0% <= strains < 15%) | 7355 |
| **Total genes** (0% <= strains <= 100%) | **10014** |
