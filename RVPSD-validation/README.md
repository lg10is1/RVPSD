# RVPSD Validation Pipeline

External validation of AlphaFold2-predicted RNA viral protein structures against experimentally determined PDB structures using Foldseek + US-align. This pipeline evaluates whether AlphaFold2 pLDDT scores reliably proxy experimental structural accuracy (TM-score / RMSD) for RNA viral proteins.

---

## Overview

| Step | Script | Purpose |
|:---|:---|:---|
| 1 | `01_foldseek_search.sh` | Search predicted PDBs against PDB100 database via Foldseek |
| 2 | `02_process_results.py` | Filter high-confidence hits, download experimental PDBs, compute TM-score/RMSD via US-align |
| 3 | `03_compute_plddt_correlation.py` | Extract mean pLDDT from B-factors, correlate with TM-score, generate stratified statistics |
| 4 | `04_plot_nature_figures.py` | Generate publication-ready Nature-style figures (hexbin scatter + stratified box plots) |

----

## Requirements

### Software
- [Foldseek](https://github.com/steineggerlab/foldseek) (conda: `bioconda::foldseek`)
- [US-align](https://zhanggroup.org/US-align/) (auto-compiled by Step 1 if `g++` available)
- Python ≥ 3.9 with: `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`

### Conda Environment
```bash
conda create -n rvpsd_val -c conda-forge -c bioconda foldseek python=3.10
conda activate rvpsd_val
pip install pandas numpy scipy matplotlib seaborn