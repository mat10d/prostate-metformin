################################################################################
# Multivariate Cox Regression with Forest Plots
# Project: Metformin overcomes the consequences of NKX3.1 loss to suppress
#          prostate cancer progression
# PI: Cory Abate-Shen
# Published: European Urology (2023)
# DOI: https://www.sciencedirect.com/science/article/pii/S0302283823030166
################################################################################

# Load required libraries
library(readxl)
library(survival)
library(survminer)
library(dplyr)
library(tidyr)
library(cowplot)
library(grid)
library(grDevices)

# Source the modified ggforest function for handling infinite CI values
source("R/ggforest_inf.R")

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

# Create full dataframe with proper variable labeling for forest plots
BCR_full <- mysheets$`NKX3.1 expression` %>% dplyr::select(`Case #`, NKX3.1) %>%
  dplyr::left_join(mysheets$`Metformin`, by = "Case #") %>%
  dplyr::left_join(mysheets$`EAU Risk` %>% dplyr::select(`Case #`, EAU_risk), by = "Case #") %>%
  dplyr::left_join(mysheets$`PSA levels`, by = "Case #") %>%
  dplyr::left_join(mysheets$`Gleason Score` %>% dplyr::select(`Case #`, g_score), by = "Case #") %>%
  dplyr::left_join(mysheets$`BCR-free Estimated Survival`, by = "Case #") %>%
  # Convert EAU risk to factor with descriptive labels
  dplyr::mutate(EAU_risk = factor(EAU_risk, levels = c(1,2,3),
                                  labels = c("Low EAU risk", "Int. EAU risk", "High EAU risk"))) %>%
  # Dichotomize PSA at 10
  dplyr::mutate(PSA = ifelse(PSA > 10, 1, 0)) %>%
  # Dichotomize Gleason score at 7
  dplyr::mutate(g_score = ifelse(g_score > 6, 1, 0)) %>%
  # Rename variables for better display in forest plots
  dplyr::rename(`High NKX3.1 expression` = NKX3.1,
                `Metformin treatment` = Metformin,
                `EAU risk group` = EAU_risk,
                `PSA > 10` = PSA,
                `GS > 6` = g_score)

################################################################################
# Multivariate Cox Regression Models
################################################################################

# Create results directory if it doesn't exist
dir.create("results/multivariate_forest", recursive = TRUE, showWarnings = FALSE)

## 1. Full Multivariate Model (All Patients)
cat("\n=== Full Multivariate Cox Regression Model ===\n")
BCR <- BCR_full %>%
  dplyr::filter(complete.cases(.))

sink('results/multivariate_forest/BCR_multi-full.txt')
res.cox <- coxph(Surv(Months,BCR) ~
                   `High NKX3.1 expression` + `Metformin treatment` +
                   `EAU risk group` + `PSA > 10` + `GS > 6`, data = BCR)
summary(res.cox)
sink()

# Save coefficients to CSV
coef_table <- as.data.frame(summary(res.cox)$coefficients)
write.csv(coef_table, "results/multivariate_forest/BCR_multi-full_coefficients.csv")

# Generate forest plot
ggforest(res.cox)
ggsave("results/multivariate_forest/BCR_multi-full.pdf", width = 9)

cat("Model fitted for", nrow(BCR), "patients\n")

## 2. Multivariate Model - Intermediate EAU Risk Stratum
cat("\n=== Multivariate Model: Intermediate EAU Risk ===\n")
BCR <- BCR_full %>%
  dplyr::filter(`EAU risk group` == "Int. EAU risk") %>%
  dplyr::filter(complete.cases(.))

sink('results/multivariate_forest/BCR_multi-EAU_int.txt')
res.cox <- coxph(Surv(Months,BCR) ~
                   `High NKX3.1 expression` + `Metformin treatment` +
                   `PSA > 10` + `GS > 6`, data = BCR)
summary(res.cox)
sink()

# Save coefficients to CSV
coef_table <- as.data.frame(summary(res.cox)$coefficients)
write.csv(coef_table, "results/multivariate_forest/BCR_multi-EAU_int_coefficients.csv")

# Generate forest plot
ggforest(res.cox)
ggsave("results/multivariate_forest/BCR_multi-EAU_int.pdf", width = 9)

cat("Model fitted for", nrow(BCR), "patients with Intermediate EAU risk\n")

## 3. Multivariate Model - High EAU Risk Stratum
cat("\n=== Multivariate Model: High EAU Risk ===\n")
BCR <- BCR_full %>%
  dplyr::filter(`EAU risk group` == "High EAU risk") %>%
  dplyr::filter(complete.cases(.))

sink('results/multivariate_forest/BCR_multi-EAU_high.txt')
res.cox <- coxph(Surv(Months,BCR) ~
                   `High NKX3.1 expression` + `Metformin treatment` +
                   `PSA > 10` + `GS > 6`, data = BCR)
summary(res.cox)
sink()

# Save coefficients to CSV
coef_table <- as.data.frame(summary(res.cox)$coefficients)
write.csv(coef_table, "results/multivariate_forest/BCR_multi-EAU_high_coefficients.csv")

# Generate forest plot using the modified ggforest_inf function to handle infinite CIs
ggforest_inf(res.cox)
ggsave("results/multivariate_forest/BCR_multi-EAU_high.pdf", width = 9)

cat("Model fitted for", nrow(BCR), "patients with High EAU risk\n")

################################################################################
# Summary
################################################################################

cat("\n=== Multivariate Forest Plot Analysis Complete ===\n")
cat("Results saved to results/multivariate_forest/\n")
cat("- Full model: BCR_multi-full.pdf and .txt\n")
cat("- Intermediate EAU risk: BCR_multi-EAU_int.pdf and .txt\n")
cat("- High EAU risk: BCR_multi-EAU_high.pdf and .txt\n")
cat("- Coefficient tables saved as CSV files\n")
