################################################################################
# 02_barplot_association.R
# Association Analysis: NKX3.1 Expression and EAU Risk Stratification
#
# Project: Metformin overcomes the consequences of NKX3.1 loss to suppress
#          prostate cancer progression
# Journal: European Urology (2024)
# DOI:     https://doi.org/10.1016/j.eururo.2023.07.016
# PI:      Cory Abate-Shen
#
# Paper Figure Mapping (Cohort 1):
#   - Grouped bar plot: NKX3.1 vs EAU risk distribution → Supplementary Fig. 6A
#     Shows patient counts by NKX3.1 expression level across EAU risk groups
#
# Outputs:
#   results/association/nkx3.1_eau_risk_counts.csv  - Summary count table
#   results/association/barplot_nkx3.1_eau.pdf      - Grouped bar plot
################################################################################

# Load required libraries
library(readxl)
library(dplyr)
library(ggplot2)

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

################################################################################
# Association Analysis: NKX3.1 vs EAU Risk → Supplementary Fig. 6A
################################################################################

cat("\n=== Association Analysis: NKX3.1 Expression and EAU Risk ===\n")

# Create dataframe with NKX3.1 and EAU risk, converting to factor labels
BCR_full <- mysheets$`NKX3.1 expression` %>% dplyr::select(`Case #`, NKX3.1) %>%
  dplyr::left_join(mysheets$`EAU Risk` %>% dplyr::select(`Case #`, EAU_risk), by = "Case #") %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::select(-`Case #`) %>%
  dplyr::mutate(EAU = as.factor(case_when(EAU_risk == 1 ~ "Low risk",
                                          EAU_risk == 2 ~ "Intermediate risk",
                                          EAU_risk == 3 ~ "High risk")),
                NKX = as.factor(case_when(NKX3.1 == 1 ~ "Low NKX3.1",
                                          NKX3.1 == 2 ~ "High NKX3.1"))) %>%
  dplyr::select(-NKX3.1, -EAU_risk) %>%
  dplyr::rename(NKX3.1 = NKX, EAU_risk = EAU) %>%
  dplyr::filter(complete.cases(.))

# Create summary table of counts
BCR_full_table <- BCR_full %>%
  group_by(NKX3.1, EAU_risk) %>%
  summarise(counts = n(), .groups = "drop")

# Display the summary table
print(BCR_full_table)

# Create results directory if it doesn't exist
dir.create("results/association", recursive = TRUE, showWarnings = FALSE)

# Save the summary table to CSV
write.csv(BCR_full_table, "results/association/nkx3.1_eau_risk_counts.csv", row.names = FALSE)

################################################################################
# Visualization: Grouped Bar Plot → Supplementary Fig. 6A
################################################################################

cat("\n=== Creating grouped bar plot ===\n")

# Create the grouped bar plot
p <- ggplot(BCR_full_table, aes(x = EAU_risk, y = counts)) +
  geom_bar(
    aes(fill = NKX3.1), color = "black",
    stat = "identity", position = position_dodge(0.8),
    width = 0.7
  ) +
  scale_color_manual(values = c("gray37", "gray70")) +
  scale_fill_manual(values = c("gray37", "gray70")) +
  geom_text(
    aes(label = counts, group = NKX3.1),
    position = position_dodge(0.8),
    vjust = -0.3, size = 3.5
  ) +
  ylim(0, 20) +
  ylab("Number of patients at risk") +
  xlab("EAU risk stratification") +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_blank())

# Save the plot
ggsave("results/association/barplot_nkx3.1_eau.pdf", plot = p, height = 5, width = 6)

cat("\n=== Association Analysis Complete ===\n")
cat("Results saved to results/association/\n")
cat("- Summary table: nkx3.1_eau_risk_counts.csv\n")
cat("- Bar plot: barplot_nkx3.1_eau.pdf\n")
