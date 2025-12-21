################################################################################
# Kaplan-Meier Survival Analysis
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

# Save the processed data
write.csv(BCR_full, "results/BCR_full_data.csv", row.names = FALSE)

################################################################################
# Individual Univariate Survival Analyses
################################################################################

# Create results directory if it doesn't exist
dir.create("results/univariate", recursive = TRUE, showWarnings = FALSE)

## 1. NKX3.1 Expression
cat("\n=== Analyzing NKX3.1 Expression ===\n")
BCR <- BCR_full %>% dplyr::select(NKX3.1, Months, BCR) %>%
  dplyr::filter(complete.cases(.))

# Cox proportional hazards model
sink('results/univariate/BCR_nkx3.1.txt')
res.cox <- coxph(Surv(Months,BCR) ~ NKX3.1, data = BCR)
summary(res.cox)
sink()

# Kaplan-Meier survival curve
fit <- survfit(Surv(Months,BCR) ~ NKX3.1, data = BCR)

BCR_nkx3.1 <- ggsurvplot(fit, data = BCR, conf.int = FALSE,
                         risk.table = TRUE, palette = c("blue", "red"),
                         legend.labs=c("Low NKX3.1", "High NKX3.1"), pval = FALSE,
                         ggtheme = theme_gray()) +
  ylab(c("BCR-free estimated survival probability")) +
  xlab(c("Time (months)"))

BCR_nkx3.1$table <- BCR_nkx3.1$table + labs(title = "Number of patients at risk") + ylab("")

BCR_nkx3.1$plot <- BCR_nkx3.1$plot +
  ggplot2::annotate("text", x = 150, y = .150,
                    label = paste("Log Rank p-value =",
                                  format(round(as.numeric(surv_pvalue(fit)[2]),digits = 4), nsmall = 4)))

pdf("results/univariate/BCR_nkx3.1.pdf", width = 7, height = 9)
print(BCR_nkx3.1, newpage = FALSE)
dev.off()

## 2. Metformin Treatment
cat("\n=== Analyzing Metformin Treatment ===\n")
BCR <- BCR_full %>% dplyr::select(Metformin, Months, BCR) %>%
  dplyr::filter(complete.cases(.))

sink('results/univariate/BCR_metformin.txt')
res.cox <- coxph(Surv(Months,BCR) ~ Metformin, data = BCR)
summary(res.cox)
sink()

fit <- survfit(Surv(Months,BCR) ~ Metformin, data = BCR)

BCR_metformin <- ggsurvplot(fit, data = BCR, conf.int = FALSE,
                            risk.table = TRUE, palette = c("blue", "red"),
                            legend.labs=c("No Metformin intervention", "Metformin intervention"),
                            pval = FALSE, ggtheme = theme_gray()) +
  ylab(c("BCR-free estimated survival probability")) +
  xlab(c("Time (months)"))

BCR_metformin$table <- BCR_metformin$table + labs(title = "Number of patients at risk") + ylab("")

BCR_metformin$plot <- BCR_metformin$plot +
  ggplot2::annotate("text", x = 150, y = .150,
                    label = paste("Log Rank p-value =",
                                  format(round(as.numeric(surv_pvalue(fit)[2]),digits = 4), nsmall = 4)))

pdf("results/univariate/BCR_metformin.pdf", width = 7, height = 9)
print(BCR_metformin, newpage = FALSE)
dev.off()

## 3. EAU Risk Stratification
cat("\n=== Analyzing EAU Risk Stratification ===\n")
BCR <- BCR_full %>% dplyr::select(EAU_risk, Months, BCR) %>%
  dplyr::filter(complete.cases(.))

sink('results/univariate/BCR_eau.txt')
res <- pairwise_survdiff(Surv(Months,BCR) ~ EAU_risk, data = BCR, p.adjust.method = "none")
res
res.cox <- coxph(Surv(Months,BCR) ~ EAU_risk, data = BCR)
summary(res.cox)
sink()

fit <- survfit(Surv(Months,BCR) ~ EAU_risk, data = BCR)

BCR_eau <- ggsurvplot(fit, data = BCR, conf.int = FALSE,
                      risk.table = TRUE,
                      palette = c("blue", "salmon","red"),
                      legend.labs=c("Low EAU risk", "Intermediate EAU risk", "High EAU risk"),
                      pval = FALSE, ggtheme = theme_gray()) +
  ylab(c("BCR-free estimated survival probability")) +
  xlab(c("Time (months)"))

BCR_eau$table <- BCR_eau$table + labs(title = "Number of patients at risk") + ylab("")

pdf("results/univariate/BCR_eau_no-p.pdf", width = 7, height = 9)
print(BCR_eau, newpage = FALSE)
dev.off()

## 4. PSA Levels (dichotomized at 10)
cat("\n=== Analyzing PSA Levels ===\n")
BCR <- BCR_full %>% dplyr::select(PSA, Months, BCR) %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::mutate(PSA = ifelse(PSA > 10, 1, 0))

sink('results/univariate/BCR_psa.txt')
res.cox <- coxph(Surv(Months,BCR) ~ PSA, data = BCR)
summary(res.cox)
sink()

fit <- survfit(Surv(Months,BCR) ~ PSA, data = BCR)

BCR_psa <- ggsurvplot(fit, data = BCR, conf.int = FALSE,
                      risk.table = TRUE,
                      palette = c("blue", "red"),
                      legend.labs=c(paste0("PSA ", "\u2264"," 10"), "PSA > 10"),
                      pval = FALSE, ggtheme = theme_gray()) +
  ylab(c("BCR-free estimated survival probability")) +
  xlab(c("Time (months)"))

BCR_psa$table <- BCR_psa$table + labs(title = "Number of patients at risk") + ylab("")

BCR_psa$plot <- BCR_psa$plot +
  ggplot2::annotate("text", x = 150, y = .150,
                    label = paste("Log Rank p-value =",
                                  format(round(as.numeric(surv_pvalue(fit)[2]),digits = 4), nsmall = 4)))

grDevices::cairo_pdf("results/univariate/BCR_psa.pdf", width = 7, height = 9)
print(BCR_psa, newpage = FALSE)
dev.off()

## 5. Gleason Score (dichotomized at 7)
cat("\n=== Analyzing Gleason Score ===\n")
BCR <- BCR_full %>% dplyr::select(g_score, Months, BCR) %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::mutate(g_score = ifelse(g_score >= 7, 1, 0))

sink('results/univariate/BCR_gleason.txt')
res.cox <- coxph(Surv(Months,BCR) ~ g_score, data = BCR)
summary(res.cox)
sink()

fit <- survfit(Surv(Months,BCR) ~ g_score, data = BCR)

BCR_gleason <- ggsurvplot(fit, data = BCR, conf.int = FALSE,
                          risk.table = TRUE,
                          palette = c("blue", "red"),
                          legend.labs=c("GS < 7", paste0("GS ", "\u2265"," 7")),
                          pval = FALSE, ggtheme = theme_gray()) +
  ylab(c("BCR-free estimated survival probability")) +
  xlab(c("Time (months)"))

BCR_gleason$table <- BCR_gleason$table + labs(title = "Number of patients at risk") + ylab("")

BCR_gleason$plot <- BCR_gleason$plot +
  ggplot2::annotate("text", x = 150, y = .150,
                    label = paste("Log Rank p-value =",
                                  format(round(as.numeric(surv_pvalue(fit)[2]),digits = 4), nsmall = 4)))

grDevices::cairo_pdf("results/univariate/BCR_gleason.pdf", width = 7, height = 9)
print(BCR_gleason, newpage = FALSE)
dev.off()

################################################################################
# Grouped Multivariate Survival Analyses
################################################################################

dir.create("results/multivariate", recursive = TRUE, showWarnings = FALSE)

## 1. NKX3.1 Expression & EAU Risk
cat("\n=== Analyzing NKX3.1 Expression & EAU Risk ===\n")
BCR <- BCR_full %>% dplyr::select(NKX3.1, EAU_risk, Months, BCR) %>%
  dplyr::filter(complete.cases(.))

sink('results/multivariate/BCR_nkx3.1_eau.txt')
res <- pairwise_survdiff(Surv(Months,BCR) ~ EAU_risk + NKX3.1, data = BCR, p.adjust.method = "none")
res
res.cox <- coxph(Surv(Months,BCR) ~ EAU_risk + NKX3.1, data = BCR)
summary(res.cox)
sink()

fit <- survfit(Surv(Months,BCR) ~ EAU_risk + NKX3.1, data = BCR)

BCR_nkx3.1_eau <- ggsurvplot(fit, data = BCR, conf.int = FALSE,
                             risk.table = TRUE,
                             palette = c("blue","blue1","salmon","salmon1","red", "red1"),
                             legend.labs=c("Low NKX3.1, Low EAU", "High NKX3.1, Low EAU",
                                           "Low NKX3.1, Int. EAU", "High NKX3.1, Int. EAU",
                                           "Low NKX3.1, High EAU","High NKX3.1, High EAU"),
                             linetype = c(1,4,1,4,1,4), pval = FALSE,
                             ggtheme = theme_gray()) +
  ylab(c("BCR-free estimated survival probability")) +
  xlab(c("Time (months)"))

BCR_nkx3.1_eau$table <- BCR_nkx3.1_eau$table + labs(title = "Number of patients at risk") + ylab("")

pdf("results/multivariate/BCR_nkx3.1_eau_no-p.pdf", width = 7, height = 9)
print(BCR_nkx3.1_eau, newpage = FALSE)
dev.off()

## 2. NKX3.1 Expression & Metformin Intervention (KEY ANALYSIS)
cat("\n=== Analyzing NKX3.1 Expression & Metformin Intervention ===\n")
BCR <- BCR_full %>% dplyr::select(NKX3.1, Metformin, Months, BCR) %>%
  dplyr::filter(complete.cases(.))

sink('results/multivariate/BCR_nkx3.1_metformin.txt')
res <- pairwise_survdiff(Surv(Months,BCR) ~ NKX3.1 + Metformin, data = BCR, p.adjust.method = "none")
res
res.cox <- coxph(Surv(Months,BCR) ~ NKX3.1 + Metformin, data = BCR)
summary(res.cox)
sink()

# Save pairwise comparison results to CSV
pairwise_results <- as.data.frame(res$p.value)
write.csv(pairwise_results, "results/multivariate/BCR_nkx3.1_metformin_pairwise.csv")

fit <- survfit(Surv(Months,BCR) ~ NKX3.1 + Metformin, data = BCR)

BCR_nkx3.1_metformin <- ggsurvplot(fit, data = BCR, conf.int = FALSE,
                                   risk.table = TRUE,
                                   palette = c("blue","blue1","red", "red1"),
                                   legend.labs=c("Low NKX3.1, No trt", "Low NKX3.1, Metformin",
                                                 "High NKX3.1, No trt", "High NKX3.1, Metformin"),
                                   linetype = c(1,4,1,4), pval = FALSE,
                                   ggtheme = theme_gray()) +
  ylab(c("BCR-free estimated survival probability")) +
  xlab(c("Time (months)"))

BCR_nkx3.1_metformin$table <- BCR_nkx3.1_metformin$table +
  labs(title = "Number of patients at risk") + ylab("")

pdf("results/multivariate/BCR_nkx3.1_metformin_no-p.pdf", width = 7, height = 9)
print(BCR_nkx3.1_metformin, newpage = FALSE)
dev.off()

## 3. EAU Risk & Metformin Intervention
cat("\n=== Analyzing EAU Risk & Metformin Intervention ===\n")
BCR <- BCR_full %>% dplyr::select(EAU_risk, Metformin, Months, BCR) %>%
  dplyr::filter(complete.cases(.))

sink('results/multivariate/BCR_eau_metformin.txt')
res <- pairwise_survdiff(Surv(Months,BCR) ~ EAU_risk + Metformin, data = BCR, p.adjust.method = "none")
res
res.cox <- coxph(Surv(Months,BCR) ~ EAU_risk + Metformin, data = BCR)
summary(res.cox)
sink()

fit <- survfit(Surv(Months,BCR) ~ EAU_risk + Metformin, data = BCR)

BCR_eau_metformin <- ggsurvplot(fit, data = BCR, conf.int = FALSE,
                                risk.table = TRUE,
                                palette = c("blue","blue1","salmon","salmon1","red", "red1"),
                                legend.labs=c("No trt, Low EAU", "Metformin, Low EAU",
                                              "No trt, Int. EAU", "Metformin, Int. EAU",
                                              "No trt, High EAU","Metformin, High EAU"),
                                linetype = c(4,1,4,1,4,1), pval = FALSE,
                                ggtheme = theme_gray()) +
  ylab(c("BCR-free estimated survival probability")) +
  xlab(c("Time (months)"))

BCR_eau_metformin$table <- BCR_eau_metformin$table + labs(title = "Number of patients at risk") + ylab("")

pdf("results/multivariate/BCR_eau_metformin_no-p.pdf", width = 7, height = 9)
print(BCR_eau_metformin, newpage = FALSE)
dev.off()

## 4. Full three-way model: NKX3.1 & EAU Risk & Metformin
cat("\n=== Analyzing NKX3.1 Expression & EAU Risk & Metformin Intervention ===\n")
BCR <- BCR_full %>% dplyr::select(NKX3.1, EAU_risk, Metformin, Months, BCR) %>%
  dplyr::filter(complete.cases(.))

sink('results/multivariate/BCR_nkx3.1_eau_metformin.txt')
res.cox <- coxph(Surv(Months,BCR) ~ NKX3.1 + EAU_risk + Metformin, data = BCR)
summary(res.cox)
sink()

# Save Cox model coefficients to CSV
coef_table <- as.data.frame(summary(res.cox)$coefficients)
write.csv(coef_table, "results/multivariate/BCR_nkx3.1_eau_metformin_coefficients.csv")

fit <- survfit(Surv(Months,BCR) ~ NKX3.1 + EAU_risk + Metformin, data = BCR)

BCR_nkx3.1_eau_metformin <- ggsurvplot(fit, data = BCR, conf.int = FALSE,
                                       risk.table = TRUE, pval = FALSE,
                                       ggtheme = theme_gray()) +
  ylab(c("BCR-free estimated survival probability")) +
  xlab(c("Time (months)"))

BCR_nkx3.1_eau_metformin$table <- BCR_nkx3.1_eau_metformin$table +
  labs(title = "Number of patients at risk") + ylab("")

pdf("results/multivariate/BCR_nkx3.1_eau_metformin_no-p.pdf", width = 7, height = 9)
print(BCR_nkx3.1_eau_metformin, newpage = FALSE)
dev.off()

cat("\n=== Survival Analysis Complete ===\n")
cat("Results saved to results/univariate/ and results/multivariate/\n")
