# PLDDT Distribution Compared Across Viral Families

## Overview

This repository contains an R script for the comparative analysis of AlphaFold-predicted pLDDT (predicted Local Distance Difference Test) scores across the top 10 most abundant viral families.

The analysis includes:

- **Single-family density plots** — pLDDT distribution for each of the top-10 families
- **Global violin + boxplot** — overview across all 10 families
- **Effect-size heatmap** — rank-biserial correlation |r| for all 45 pairwise comparisons
- **Unbalanced pairwise overlays** — density comparison for all 45 family pairs (full data)
- **Balanced 1:1 pairwise overlays** — density comparison using stratified sampling to control for unequal sample sizes
- **Robustness scatter plot** — balanced vs unbalanced effect-size correlation
- **Bilingual analysis report** — comprehensive summary with biological interpretation

## File Structure

```
.
├── PLDDT_Family_Analysis.R   # Main analysis script (self-contained)
├── example_data.csv           # Synthetic example dataset (~500 rows)
├── README.md                  # This file


## Input Data Format

The script expects a CSV file named `virus_new.csv` in the working directory with three columns:

```
protein_name,family,plddt
synthetic_protein_001,Retroviridae,82.3
synthetic_protein_002,Picornaviridae,91.5
...
```

- `protein_name` — identifier for the protein (can be anonymised)
- `family` — ICTV-approved viral family name
- `plddt` — AlphaFold pLDDT confidence score (0–100)

**Note:** `example_data.csv` is a synthetic dataset included for demonstration. It contains 500 rows across the same 10 families, with pLDDT values approximated from the real distributional parameters (median, mean, SD). Replace it with your own dataset for actual analysis.

## Installation & Dependencies

The script automatically installs required packages on first run. Optionally pre-install:

```r
install.packages(c("data.table", "ggplot2", "patchwork", "dplyr", "tidyr"),
                 repos = "https://cloud.r-project.org")
```

R version ≥ 4.0 recommended.

## Usage

```bash
# Run from the directory containing virus_new.csv
cd /path/to/your/data/
Rscript PLDDT_Family_Analysis.R
```

Or in R:

```r
source("PLDDT_Family_Analysis.R")
```

## Statistical Methods

| Method | Details |
|---|---|
| Test | Wilcoxon rank-sum test (Mann-Whitney U), exact = FALSE |
| Multiple testing correction | Bonferroni (45 comparisons) |
| Effect size | Rank-biserial correlation r = (2U − n₁n₂)/(n₁n₂) |
| Effect size magnitude | Negligible (< 0.1), Small (0.1–0.3), Medium (0.3–0.5), Large (≥ 0.5) |
| Balanced 1:1 sampling | Stratified by pLDDT into 20 bins; min(n₁, n₂) drawn per bin per family |

## Output Key Files

| File | Description |
|---|---|
| `figures/overview/Top10_ViolinBox.pdf` | Main text figure — violin + boxplot overview |
| `figures/overview/Effsize_Heatmap.pdf` | Supplementary — effect-size heatmap |
| `figures/overview/Balanced_vs_Unbalanced.pdf` | Supplementary — robustness check |
| `results/ANALYSIS_BILINGUAL_REPORT.txt` | Full bilingual analysis report |
| `results/pairwise_wilcox_unbalanced.csv` | All 45-pair statistics (unbalanced) |
| `results/pairwise_balanced_all_45.csv` | All 45-pair statistics (balanced 1:1) |

## Notes

- pLDDT is a **model confidence metric**, not a direct measure of true structural quality. Interpret with caution.
- The balanced 1:1 analysis is essential when sample sizes differ by >10-fold between families, to avoid statistical over-power from inflating significance.
- All viral family names follow **ICTV italicisation conventions** in the figures.

## License

MIT License — free to use and modify with attribution.