################################################################################
# 03_correlation_analysis.R
# Spearman Correlation Analysis
#
# Project: Metformin overcomes the consequences of NKX3.1 loss to suppress
#          prostate cancer progression
# Journal: European Urology (2024)
# DOI:     https://doi.org/10.1016/j.eururo.2023.07.016
# PI:      Cory Abate-Shen
#
# Paper Figure Mapping:
#   This script tests for confounding between EAU risk, NKX3.1 expression,
#   and metformin treatment. Results confirm independence of covariates,
#   supporting the validity of the survival analyses in Fig. 5I and 5K.
#   These correlation tests are referenced in the Methods section and
#   support the claims in Table 1.
#
# Outputs:
#   results/correlation/EAU_risk.NKX3.1.txt         - Spearman test output
#   results/correlation/EAU_risk.NKX3.1.csv         - Correlation statistics
#   results/correlation/EAU_risk.Metformin.txt       - Spearman test output
#   results/correlation/EAU_risk.Metformin.csv       - Correlation statistics
#   results/correlation/all_correlations_summary.csv - Combined summary
################################################################################

# Load required libraries
library(readxl)
library(dplyr)
library(Hmisc)

################################################################################
# Helper Functions
################################################################################

#' Read all sheets from an Excel file
#'
#' @param filename Path to the Excel file
#' @param tibble Whether to return tibbles (default: FALSE for data.frames)
#' @return Named list of data frames, one per sheet
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

################################################################################
# Data Loading and Preprocessing
################################################################################

# Read all Excel sheets
mysheets <- read_excel_allsheets("data/Metformin_data.xlsx")

# Convert Gleason Score from Y/N columns to numeric score
mysheets$`Gleason Score`[is.na(mysheets$`Gleason Score`)] <- "N"

for(i in 1:nrow(mysheets$`Gleason Score`)){
  if(mysheets$`Gleason Score`$five[i] == "Y"){
    mysheets$`Gleason Score`$g_score[i] <- 5
  }
  if(mysheets$`Gleason Score`$six[i] == "Y"){
    mysheets$`Gleason Score`$g_score[i] <- 6
  }
  if(mysheets$`Gleason Score`$seven[i] == "Y"){
    mysheets$`Gleason Score`$g_score[i] <- 7
  }
  if(mysheets$`Gleason Score`$eight[i] == "Y"){
    mysheets$`Gleason Score`$g_score[i] <- 8
  }
  if(mysheets$`Gleason Score`$nine[i] == "Y"){
    mysheets$`Gleason Score`$g_score[i] <- 9
  }
  if(mysheets$`Gleason Score`$five[i] == "N" & mysheets$`Gleason Score`$six[i] == "N" &
     mysheets$`Gleason Score`$seven[i] == "N" & mysheets$`Gleason Score`$eight[i] == "N" &
     mysheets$`Gleason Score`$nine[i] == "N"){
    mysheets$`Gleason Score`$g_score[i] <- NA
  }
}

# Convert EAU Risk from Y/N columns to numeric categories
mysheets$`EAU Risk`[is.na(mysheets$`EAU Risk`)] <- "N"

for(i in 1:nrow(mysheets$`EAU Risk`)){
  if(mysheets$`EAU Risk`$Low[i] == "Y"){
    mysheets$`EAU Risk`$EAU_risk[i] <- 1
  }
  if(mysheets$`EAU Risk`$Intermediate[i] == "Y"){
    mysheets$`EAU Risk`$EAU_risk[i] <- 2
  }
  if(mysheets$`EAU Risk`$High[i] == "Y"){
    mysheets$`EAU Risk`$EAU_risk[i] <- 3
  }
  if(mysheets$`EAU Risk`$Low[i] == "N" & mysheets$`EAU Risk`$Intermediate[i] == "N" &
     mysheets$`EAU Risk`$High[i] == "N"){
    mysheets$`EAU Risk`$EAU_risk[i] <- NA
  }
}

# Convert NKX3.1 expression from Y/N columns to numeric categories
mysheets$`NKX3.1 expression`[is.na(mysheets$`NKX3.1 expression`)] <- "N"

for(i in 1:nrow(mysheets$`NKX3.1 expression`)){
  if(mysheets$`NKX3.1 expression`$`NKX3.1 Low`[i] == "Y"){
    mysheets$`NKX3.1 expression`$NKX3.1[i] <- 1
  }
  if(mysheets$`NKX3.1 expression`$`NKX3.1 High`[i] == "Y"){
    mysheets$`NKX3.1 expression`$NKX3.1[i] <- 2
  }
  if(mysheets$`NKX3.1 expression`$`NKX3.1 Low`[i] == "N" &
     mysheets$`NKX3.1 expression`$`NKX3.1 High`[i] == "N"){
    mysheets$`NKX3.1 expression`$NKX3.1[i] <- NA
  }
}

# Create full dataframe by joining all relevant sheets
BCR_full <- mysheets$`NKX3.1 expression` %>% dplyr::select(`Case #`, NKX3.1) %>%
  dplyr::left_join(mysheets$`Metformin`, by = "Case #") %>%
  dplyr::left_join(mysheets$`EAU Risk` %>% dplyr::select(`Case #`, EAU_risk), by = "Case #") %>%
  dplyr::left_join(mysheets$`PSA levels`, by = "Case #") %>%
  dplyr::left_join(mysheets$`Gleason Score` %>% dplyr::select(`Case #`, g_score), by = "Case #") %>%
  dplyr::left_join(mysheets$`BCR-free Estimated Survival`, by = "Case #")

# NOTE: readxl may import some numeric columns as character type.
# Ensure numeric types for correlation analysis.
BCR_full$Metformin <- as.numeric(BCR_full$Metformin)

################################################################################
# Spearman Correlation Analysis
# Tests independence of covariates (supports Table 1 claims)
################################################################################

# Create results directory if it doesn't exist
dir.create("results/correlation", recursive = TRUE, showWarnings = FALSE)

## 1. EAU Risk vs NKX3.1 Expression
cat("\n=== Spearman Correlation: EAU Risk vs NKX3.1 Expression ===\n")

sink('results/correlation/EAU_risk.NKX3.1.txt')
res <- cor.test(BCR_full$EAU_risk, BCR_full$NKX3.1, method = "spearman")
print(res)
sink()

# Extract correlation statistics and save to CSV
corr_results_nkx <- data.frame(
  Variable1 = "EAU_risk",
  Variable2 = "NKX3.1",
  Rho = res$estimate,
  P_value = res$p.value,
  Method = "Spearman",
  Alternative = res$alternative
)
write.csv(corr_results_nkx, "results/correlation/EAU_risk.NKX3.1.csv", row.names = FALSE)

cat("Spearman's rho:", res$estimate, "\n")
cat("P-value:", res$p.value, "\n")

## 2. EAU Risk vs Metformin Treatment
cat("\n=== Spearman Correlation: EAU Risk vs Metformin Treatment ===\n")

sink('results/correlation/EAU_risk.Metformin.txt')
res <- cor.test(BCR_full$EAU_risk, BCR_full$Metformin, method = "spearman")
print(res)
sink()

# Extract correlation statistics and save to CSV
corr_results_met <- data.frame(
  Variable1 = "EAU_risk",
  Variable2 = "Metformin",
  Rho = res$estimate,
  P_value = res$p.value,
  Method = "Spearman",
  Alternative = res$alternative
)
write.csv(corr_results_met, "results/correlation/EAU_risk.Metformin.csv", row.names = FALSE)

cat("Spearman's rho:", res$estimate, "\n")
cat("P-value:", res$p.value, "\n")

################################################################################
# Combined Correlation Summary
################################################################################

# Combine all correlation results into a single summary table
all_correlations <- rbind(corr_results_nkx, corr_results_met)
write.csv(all_correlations, "results/correlation/all_correlations_summary.csv", row.names = FALSE)

cat("\n=== Correlation Analysis Complete ===\n")
cat("Results saved to results/correlation/\n")
cat("- EAU_risk vs NKX3.1: EAU_risk.NKX3.1.txt and .csv\n")
cat("- EAU_risk vs Metformin: EAU_risk.Metformin.txt and .csv\n")
cat("- Summary: all_correlations_summary.csv\n")
