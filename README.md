# Metformin Overcomes NKX3.1 Loss to Suppress Prostate Cancer Progression

## Citation

**Title:** Metformin overcomes the consequences of NKX3.1 loss to suppress prostate cancer progression

**Journal:** European Urology (2023)

**DOI:** [10.1016/j.eururo.2023.07.016](https://doi.org/10.1016/j.eururo.2023.07.016)

**Principal Investigator:** Cory Abate-Shen, PhD

## Project Overview

This repository contains the statistical analysis code for a retrospective clinical study investigating the relationship between NKX3.1 expression, metformin treatment, and biochemical recurrence (BCR) in prostate cancer patients. The analysis demonstrates that metformin intervention can overcome the negative prognostic impact of low NKX3.1 expression on disease progression.

### Key Findings

- Low NKX3.1 expression is associated with worse biochemical recurrence-free survival
- Metformin treatment significantly improves outcomes in patients with low NKX3.1 expression
- The protective effect of metformin is particularly pronounced in patients with low NKX3.1 levels
- Multivariate analysis confirms these associations are independent of EAU risk stratification, PSA levels, and Gleason score

## Repository Structure

```
.
├── README.md                          # This file
├── data/
│   └── Metformin_data.xlsx            # Clinical data (de-identified)
├── R/
│   ├── 01_survival_analysis.R         # Kaplan-Meier survival curves
│   ├── 02_barplot_association.R       # NKX3.1 vs EAU risk association
│   ├── 03_correlation_analysis.R      # Spearman correlation tests
│   ├── 04_multivariate_forest_plots.R # Cox regression forest plots
│   └── ggforest_inf.R                 # Helper function for forest plots
└── results/                           # Output directory (created by scripts)
    ├── univariate/                    # Individual survival analyses
    ├── multivariate/                  # Grouped survival analyses
    ├── multivariate_forest/           # Cox regression models
    ├── association/                   # Association analysis
    └── correlation/                   # Correlation tests
```

## Analysis Components

### 1. Survival Analysis (`01_survival_analysis.R`)

Performs Kaplan-Meier survival analyses and Cox proportional hazards regression for biochemical recurrence-free survival.

**Univariate analyses:**
- NKX3.1 expression (Low vs High)
- Metformin treatment (Yes vs No)
- EAU risk stratification (Low/Intermediate/High)
- PSA levels (≤10 vs >10)
- Gleason score (<7 vs ≥7)

**Multivariate analyses:**
- NKX3.1 × EAU Risk
- **NKX3.1 × Metformin (KEY ANALYSIS)**
- EAU Risk × Metformin
- NKX3.1 × EAU Risk × Metformin

**Outputs:**
- Kaplan-Meier survival curves (PDF)
- Cox regression statistics (TXT)
- Pairwise comparison p-values (CSV)
- Cox coefficients (CSV)

### 2. Association Analysis (`02_barplot_association.R`)

Examines the distribution of NKX3.1 expression across EAU risk categories.

**Outputs:**
- Grouped bar plot showing patient counts (PDF)
- Summary table of counts (CSV)

### 3. Correlation Analysis (`03_correlation_analysis.R`)

Tests for associations between clinical variables using Spearman correlation.

**Tests:**
- EAU Risk vs NKX3.1 expression
- EAU Risk vs Metformin treatment

**Outputs:**
- Spearman correlation statistics (TXT and CSV)
- Combined summary table (CSV)

### 4. Multivariate Forest Plots (`04_multivariate_forest_plots.R`)

Generates forest plots for multivariate Cox regression models.

**Models:**
- Full cohort (all EAU risk groups)
- Intermediate EAU risk stratum
- High EAU risk stratum

**Outputs:**
- Forest plots (PDF)
- Cox regression summaries (TXT)
- Coefficient tables (CSV)

## Prerequisites

### R Version
- R ≥ 4.0.0

### Required R Packages

```r
# Install required packages
install.packages(c(
  "readxl",      # Reading Excel files
  "survival",    # Survival analysis
  "survminer",   # Survival plot visualization
  "dplyr",       # Data manipulation
  "tidyr",       # Data tidying
  "cowplot",     # Plot composition
  "grid",        # Graphics
  "grDevices",   # Graphics devices
  "ggplot2",     # Plotting
  "Hmisc",       # Correlation analysis
  "corrplot",    # Correlation visualization
  "broom"        # Model tidying (for ggforest)
))
```

## Running the Analysis

### Option 1: Run All Analyses

From the project root directory:

```r
# Set working directory to repository root
setwd("/path/to/git_repo")

# Run all analyses in sequence
source("R/01_survival_analysis.R")
source("R/02_barplot_association.R")
source("R/03_correlation_analysis.R")
source("R/04_multivariate_forest_plots.R")
```

### Option 2: Run Individual Analyses

```r
# Set working directory
setwd("/path/to/git_repo")

# Run only survival analysis
source("R/01_survival_analysis.R")

# Or run only correlation analysis
source("R/03_correlation_analysis.R")
```

### Expected Runtime

- 01_survival_analysis.R: ~2-3 minutes
- 02_barplot_association.R: ~10 seconds
- 03_correlation_analysis.R: ~10 seconds
- 04_multivariate_forest_plots.R: ~1 minute

**Total runtime:** ~5 minutes on a standard desktop computer

## Data Description

The `Metformin_data.xlsx` file contains multiple sheets with patient data:

- **NKX3.1 expression**: Low/High NKX3.1 IHC staining
- **Metformin**: Treatment status (0=No, 1=Yes)
- **EAU Risk**: European Association of Urology risk classification
- **PSA levels**: Pre-treatment PSA values
- **Gleason Score**: Histological grade (5-9)
- **BCR-free Estimated Survival**: Time to biochemical recurrence and event status

All patient identifiers have been removed for privacy protection.

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| NKX3.1 | Binary | 1=Low expression, 2=High expression |
| Metformin | Binary | 0=No treatment, 1=Metformin treatment |
| EAU_risk | Categorical | 1=Low, 2=Intermediate, 3=High |
| PSA | Continuous | Pre-treatment PSA (ng/mL) |
| g_score | Categorical | Gleason score (5-9) |
| Months | Continuous | Follow-up time (months) |
| BCR | Binary | Biochemical recurrence event (0=No, 1=Yes) |

## Output Files

All analysis outputs are saved to the `results/` directory:

### Survival Analysis
- `results/univariate/BCR_*.pdf` - Kaplan-Meier curves
- `results/univariate/BCR_*.txt` - Cox regression statistics
- `results/multivariate/BCR_*_pairwise.csv` - Pairwise p-values
- `results/multivariate/BCR_*_coefficients.csv` - Cox coefficients

### Association Analysis
- `results/association/barplot_nkx3.1_eau.pdf` - Bar plot
- `results/association/nkx3.1_eau_risk_counts.csv` - Count table

### Correlation Analysis
- `results/correlation/EAU_risk.NKX3.1.*` - NKX3.1 correlation
- `results/correlation/EAU_risk.Metformin.*` - Metformin correlation
- `results/correlation/all_correlations_summary.csv` - Summary table

### Multivariate Forest Plots
- `results/multivariate_forest/BCR_multi-*.pdf` - Forest plots
- `results/multivariate_forest/BCR_multi-*_coefficients.csv` - Coefficients

## Notes

- All p-values from pairwise survival comparisons are unadjusted
- PSA is dichotomized at 10 ng/mL for multivariate analyses
- Gleason score is dichotomized at 7 for multivariate analyses
- Missing data are handled by complete case analysis
- The `ggforest_inf.R` helper function handles infinite confidence intervals in forest plots

## Reproducibility

This repository is designed for complete reproducibility. All analyses can be re-run from the source data to regenerate all figures and tables in the manuscript.

## Contact

For questions about the analysis or code, please contact:
- **PI:** Cory Abate-Shen, PhD
- **Institution:** Columbia University Irving Medical Center

## License

This code is provided for research purposes. If you use this code, please cite the original publication:

> [Full citation to be added upon publication]

## Acknowledgments

Analysis performed by Matteo Di Bernardo as part of the Abate-Shen laboratory.
