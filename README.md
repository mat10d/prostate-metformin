# Metformin Overcomes the Consequences of NKX3.1 Loss to Suppress Prostate Cancer Progression

## Citation

Papachristodoulou A, Heidegger I, Virk RK, Di Bernardo M, Kim JY, Laplaca C, Picech F, Schäfer G, De Castro GJ, Hibshoosh H, Loda M, Klocker H, Rubin MA, Zheng T, Benson MC, McKiernan JM, Dutta A, Abate-Shen C. Metformin overcomes the consequences of NKX3.1 loss to suppress prostate cancer progression. *European Urology*. 2024;85(1):44–58. doi: [10.1016/j.eururo.2023.07.016](https://doi.org/10.1016/j.eururo.2023.07.016)

## Study Overview

This repository contains the statistical analysis code for the clinical correlative component of a study investigating the relationship between NKX3.1 expression, metformin treatment, and biochemical recurrence (BCR) in prostate cancer patients.

**NKX3.1** is a prostate-specific homeobox gene that functions in mitochondria to protect against aberrant oxidative stress. The study demonstrates that metformin, an antidiabetic drug with known anticancer effects related to its antioxidant activity, can overcome the adverse consequences of NKX3.1 loss for prostate cancer progression.

### Key Clinical Findings

- Low NKX3.1 expression is associated with worse BCR-free survival
- Metformin treatment significantly improves outcomes in patients with low NKX3.1 expression (cohort 1: log-rank p = 0.005)
- High NKX3.1 patients show no difference in survival regardless of metformin use
- Multivariate Cox regression confirms these associations are independent of EAU risk, PSA, and Gleason score

### Context

The analyses in this repository correspond to **Cohort 1** (Columbia University) in the paper. The paper also includes Cohort 2 (Medical University Innsbruck) and an active surveillance cohort, which are analyzed separately.

## Data Description

The `data/Metformin_data.xlsx` file contains de-identified patient data from radical prostatectomy cases across multiple sheets:

| Sheet | Description |
|-------|-------------|
| `NKX3.1 expression` | IHC-based NKX3.1 staining (Low/High) |
| `Metformin` | Treatment status (0 = No, 1 = Yes) |
| `EAU Risk` | European Association of Urology risk classification (Low/Intermediate/High) |
| `PSA levels` | Pre-treatment PSA values (ng/mL) |
| `Gleason Score` | Histological grade (5–9) |
| `BCR-free Estimated Survival` | Time to biochemical recurrence (months) and event status |

### Key Variables

| Variable | Type | Values |
|----------|------|--------|
| NKX3.1 | Binary | 1 = Low expression, 2 = High expression |
| Metformin | Binary | 0 = No treatment, 1 = Metformin treatment |
| EAU_risk | Ordinal | 1 = Low, 2 = Intermediate, 3 = High |
| PSA | Continuous | Pre-treatment PSA (ng/mL); dichotomized at 10 for multivariate models |
| g_score | Ordinal | Gleason score (5–9); dichotomized at 7 for multivariate models |
| Months | Continuous | Follow-up time (months) |
| BCR | Binary | Biochemical recurrence event (0 = No, 1 = Yes) |

> **Note:** `readxl` may import some numeric columns (Months, BCR, PSA, Metformin) as character type. The scripts include explicit `as.numeric()` conversion to handle this.

## Repository Structure

```
.
├── README.md                          # This file
├── data/
│   └── Metformin_data.xlsx            # Clinical data (de-identified)
└── R/
    ├── 01_survival_analysis.R         # Kaplan-Meier survival curves & Cox regression
    ├── 02_barplot_association.R       # NKX3.1 vs EAU risk distribution
    ├── 03_correlation_analysis.R      # Spearman correlation tests
    ├── 04_multivariate_forest_plots.R # Cox regression forest plots
    └── ggforest_inf.R                 # Modified ggforest for infinite CIs
```

## Script-to-Figure Mapping

Each script produces outputs corresponding to specific figures, tables, or supplementary materials in the published paper. All outputs are for **Cohort 1** (Columbia University).

| Script | Paper Figure/Table | Description |
|--------|-------------------|-------------|
| `01_survival_analysis.R` | **Fig. 5I** | KM curves: NKX3.1 × Metformin interaction (key result, p = 0.005 for low NKX3.1) |
| `01_survival_analysis.R` | Supplementary Fig. 6C | Univariate KM curves: NKX3.1, Metformin, EAU Risk, PSA, Gleason Score |
| `01_survival_analysis.R` | Supplementary Fig. 6A | KM curves: NKX3.1 × EAU Risk interaction |
| `02_barplot_association.R` | Supplementary Fig. 6A | Bar plot: NKX3.1 distribution across EAU risk groups |
| `03_correlation_analysis.R` | Table 1 (supporting) | Spearman correlations confirming covariate independence |
| `04_multivariate_forest_plots.R` | **Fig. 5K** | Forest plot: multivariate Cox regression (full cohort) |
| `04_multivariate_forest_plots.R` | Supplementary Table 4 | Forest plots: EAU risk–stratified Cox models |
| `ggforest_inf.R` | *(utility)* | Modified `survminer::ggforest()` handling infinite CIs |

### Detailed Output Map

#### `01_survival_analysis.R` — Survival Analysis

**Univariate analyses** (→ Supplementary Fig. 6C):
- `results/univariate/BCR_nkx3.1.pdf` — NKX3.1 Low vs High
- `results/univariate/BCR_metformin.pdf` — Metformin vs No Metformin
- `results/univariate/BCR_eau_no-p.pdf` — EAU risk (Low/Int/High)
- `results/univariate/BCR_psa.pdf` — PSA ≤10 vs >10
- `results/univariate/BCR_gleason.pdf` — Gleason <7 vs ≥7

**Multivariate analyses**:
- `results/multivariate/BCR_nkx3.1_metformin_no-p.pdf` — **Fig. 5I** (key result)
- `results/multivariate/BCR_nkx3.1_eau_no-p.pdf` — Supplementary Fig. 6A
- `results/multivariate/BCR_eau_metformin_no-p.pdf` — EAU × Metformin
- `results/multivariate/BCR_nkx3.1_eau_metformin_no-p.pdf` — Three-way model

#### `02_barplot_association.R` — Association Analysis

- `results/association/barplot_nkx3.1_eau.pdf` — Supplementary Fig. 6A (bar plot)
- `results/association/nkx3.1_eau_risk_counts.csv` — Patient counts by group

#### `03_correlation_analysis.R` — Correlation Analysis

- `results/correlation/all_correlations_summary.csv` — Summary of Spearman tests
- `results/correlation/EAU_risk.NKX3.1.txt` — EAU Risk vs NKX3.1 (supports Table 1)
- `results/correlation/EAU_risk.Metformin.txt` — EAU Risk vs Metformin

#### `04_multivariate_forest_plots.R` — Forest Plots

- `results/multivariate_forest/BCR_multi-full.pdf` — **Fig. 5K** (all patients)
- `results/multivariate_forest/BCR_multi-EAU_int.pdf` — Intermediate EAU stratum
- `results/multivariate_forest/BCR_multi-EAU_high.pdf` — High EAU stratum

### `ggforest_inf.R` — Utility Function

A modified version of `survminer::ggforest()` that gracefully handles Cox regression models where some covariates produce infinite confidence intervals (e.g., zero events in one level of a factor). This occurs in the High EAU risk subgroup due to small sample size. The modification replaces infinite CI bounds with the maximum/minimum finite CI values to allow rendering.

Original authors: Przemyslaw Biecek, Fabian Scheipl (`survminer` package).

## Requirements

### R Version

- R ≥ 4.0.0

### Required R Packages

```r
install.packages(c(
  "readxl",      # Reading Excel files
  "survival",    # Survival analysis (Surv, coxph, survfit)
  "survminer",   # Survival visualization (ggsurvplot, ggforest)
  "dplyr",       # Data manipulation
  "tidyr",       # Data tidying
  "cowplot",     # Plot composition
  "ggplot2",     # Plotting
  "Hmisc",       # Spearman correlation
  "corrplot",    # Correlation visualization
  "broom",       # Model tidying (used by ggforest_inf)
  "ggpubr"       # Publication-ready plots (used by ggforest_inf)
))
```

## Usage

### Run All Analyses

From the repository root directory:

```r
setwd("/path/to/prostate-metformin")

source("R/01_survival_analysis.R")
source("R/02_barplot_association.R")
source("R/03_correlation_analysis.R")
source("R/04_multivariate_forest_plots.R")
```

Or from the command line:

```bash
cd prostate-metformin
Rscript R/01_survival_analysis.R
Rscript R/02_barplot_association.R
Rscript R/03_correlation_analysis.R
Rscript R/04_multivariate_forest_plots.R
```

### Expected Runtime

| Script | Time |
|--------|------|
| `01_survival_analysis.R` | ~2–3 min |
| `02_barplot_association.R` | ~10 sec |
| `03_correlation_analysis.R` | ~10 sec |
| `04_multivariate_forest_plots.R` | ~1 min |
| **Total** | **~5 min** |

## Notes

- All p-values from pairwise survival comparisons are unadjusted
- PSA is dichotomized at 10 ng/mL and Gleason score at 7 for multivariate analyses
- Missing data are handled by complete case analysis
- The `ggforest_inf.R` helper handles infinite confidence intervals that arise in small subgroup analyses
- Some columns may be imported as character by `readxl`; scripts include explicit numeric coercion

## Authors

Alexandros Papachristodoulou, Isabel Heidegger, Renu K. Virk, Matteo Di Bernardo, Jaime Y. Kim, Caroline Laplaca, Florencia Picech, Georg Schäfer, Guarionex Joel De Castro, Hanina Hibshoosh, Massimo Loda, Helmut Klocker, Mark A. Rubin, Tian Zheng, Mitchell C. Benson, James M. McKiernan, Aditya Dutta, Cory Abate-Shen

## License

This code is provided for research purposes. If you use this code or data, please cite the original publication (see Citation above).
