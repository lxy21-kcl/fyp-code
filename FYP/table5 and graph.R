library(ggplot2)
library(boot)

# Read in the dataset and set encoding
ds <- read.csv("datainfo_transpose_1.20.csv", stringsAsFactors = FALSE, fileEncoding = "latin1")

# Subset the dataset to Phase 3
ds <- ds[ds$Study.Phase == "3", ]

# (Note: The following line overwrites the previous subset to only keep Phase 2;
#  it might be intentional or an error, but is left as-is.)
ds <- ds[ds$Study.Phase == "2", ]

# Create ORds dataframe
ORds <- as.data.frame(cbind(
  ds$ï..PMID, ds$Cancer.Type, ds$Blinded.open.label, ds$Sample.size, 
  ds$BICR.ORR.experimental, ds$X95..CI.BICR.ORR.L.experimental, ds$X95..CI.BICR.ORR.U.experimental,
  ds$BICR.ORR.control, ds$X95..CI.BICR.ORR.L.control, ds$X95..CI.BICR.ORR.U.control,
  ds$INV.ORR.experimental, ds$X95..CI.INV.ORR.L.experimental, ds$X95..CI.INV.ORR.U.experimental,
  ds$INV.ORR.control, ds$X95..CI.INV.ORR.L.control, ds$X95..CI.INV.ORR.U.control
))

# Convert to numeric/character
num_cols <- 4:16  # These columns are numeric
char_cols <- c(1, 2, 3)  # These columns are character
for (i in num_cols) {
  ORds[, i] <- as.numeric(ORds[, i])
}
for (i in char_cols) {
  ORds[, i] <- as.character(ORds[, i])
}

# Rename columns
colnames(ORds) <- c(
  "PMID", "Cancer_type", "Label", "Sample_size", 
  "OR_BICR_exp", "OR_BICR_L_exp", "OR_BICR_U_exp",
  "OR_BICR_c", "OR_BICR_L_c", "OR_BICR_U_c",
  "OR_INV_exp", "OR_INV_L_exp", "OR_INV_U_exp",
  "OR_INV_c", "OR_INV_L_c", "OR_INV_U_c"
)

# Calculate OR for BICR and INV
ORds$OR_BICR <- (ORds$OR_BICR_exp * (1 - ORds$OR_BICR_exp)) / (ORds$OR_BICR_c * (1 - ORds$OR_BICR_c))
ORds$OR_INV <- (ORds$OR_INV_exp * (1 - ORds$OR_INV_exp)) / (ORds$OR_INV_c * (1 - ORds$OR_INV_c))

# Remove NA values
ORds <- ORds[!is.na(ORds$OR_BICR) & !is.na(ORds$OR_INV), ]

# Calculate Pearson correlation coefficient (r) and 95% CI
cor_test <- cor.test(log(ORds$OR_BICR), log(ORds$OR_INV), method = "pearson")
r_value <- cor_test$estimate        # correlation coefficient r
r_ci <- cor_test$conf.int           # 95% confidence interval

r_value
r_ci
cat("Pearson r:", round(r_value, 3), "95% CI:", round(r_ci, 3), "\n")

# Calculate R² and 95% confidence intervals (using weighted regression)
lm_model <- lm(log(ORds$OR_INV) ~ log(ORds$OR_BICR), weights = ORds$Sample_size)
R2_value <- summary(lm_model)$r.squared  # R²

# Bootstrap calculation of the R² confidence interval
rsq_function <- function(data, indices) {
  d <- data[indices, ]
  model <- lm(log(OR_INV) ~ log(OR_BICR), data = d, weights = d$Sample_size)
  return(summary(model)$r.squared)
}

boot_result <- boot(ORds, rsq_function, R = 1000)
R2_ci <- boot.ci(boot_result, type = "perc")$percent[4:5]  # 95% CI

boot_result
R2_ci
cat("R²:", round(R2_value, 3), "95% CI:", round(R2_ci, 3), "\n")

# Calculate OddsRR and 95% CI
ORds$OddsRR <- log(ORds$OR_INV) / log(ORds$OR_BICR)
mean_OR <- mean(ORds$OddsRR, na.rm = TRUE)
ci_OR <- t.test(ORds$OddsRR)$conf.int

cat("OddsRR:", round(mean_OR, 3), "95% CI:", round(ci_OR, 3), "\n")

# Plot scatter + weighted regression line
# The 95% CI of the OR regression line is shown with shading.
library(ggplot2)
library(dplyr)

# Merge "blinded" and "double-blind" into "blinded"
ORds <- ORds %>%
  mutate(Label = case_when(
    Label %in% c("blinded", "double-blind") ~ "blinded",
    TRUE ~ Label  # Keep other labels unchanged
  ))

ggplot(ORds, aes(x = log(OR_INV), y = log(OR_BICR), size = Sample_size, color = Label)) +
  geom_point(alpha = 0.3) +
  geom_smooth(
    method = "lm", se = TRUE, color = "black", size = 0.5,
    inherit.aes = FALSE, aes(x = log(OR_INV), y = log(OR_BICR))
  ) +
  scale_size_continuous(name = "Sample Size", range = c(2, 10)) +
  scale_color_manual(values = c("open label" = "red", "blinded" = "blue")) +
  labs(
    title = "Correlation between OR_BICR and OR_INV",
    subtitle = paste0(
      "r = ", round(r_value, 3), 
      " (95% CI: ", round(r_ci[1], 3), "-", round(r_ci[2], 3), ")",
      "\nR² = ", round(R2_value, 3), 
      " (95% CI: ", round(R2_ci[1], 3), ")"
    ),
    x = "Log LE OR for Response",
    y = "Log BICR OR for Response",
    color = "Study Type"
  ) +
  theme_minimal() +
  guides(size = "none")  # Hide size legend

# The 95% CI of the HR regression line is shown with shading.
HRds <- as.data.frame(cbind(
  ds$ï..PMID, ds$Sample.size, ds$HR.for.BICR.PFS, ds$X95..CI.BICR.HR.L, ds$X95..CI.BICR.HR.U,
  ds$HR.for.INV.PFS, ds$X95..CI.INV.L, ds$X95..CI.INV.U, ds$follow.up.time, ds$Cancer.Type, ds$Blinded.open.label
))

print(ncol(HRds))

num_cols <- 2:9   # These columns are numeric
char_cols <- c(1, 10, 11)  # These columns are character
for (i in num_cols) {
  HRds[, i] <- as.numeric(HRds[, i])
}
for (i in char_cols) {
  HRds[, i] <- as.character(HRds[, i])
}

colnames(HRds) <- c(
  "PMID", "Sample_size", "HR_BICR", "HR_BICR_L", "HR_BICR_U",
  "HR_INV", "HR_INV_L", "HR_INV_U", "follow-up time", "Cancer_type", "Label"
)

HRds <- HRds[!is.na(HRds$HR_BICR) & !is.na(HRds$HR_INV), ]

# Calculate Pearson correlation coefficient (r) and 95% CI
cor_test <- cor.test(log(HRds$HR_INV), log(HRds$HR_BICR), method = "pearson")
r_value <- cor_test$estimate   # correlation coefficient r
r_ci <- cor_test$conf.int      # 95% confidence interval

r_value
r_ci
cat("Pearson r:", round(r_value, 3), "95% CI:", round(r_ci, 3), "\n")

# Calculate R² and 95% confidence intervals (using weighted regression)
lm_model <- lm(log(HRds$HR_INV) ~ log(HRds$HR_BICR), weights = HRds$Sample_size)
R2_value <- summary(lm_model)$r.squared  # R²

# Bootstrap calculation of the R² confidence interval
rsq_function <- function(data, indices) {
  d <- data[indices, ]
  model <- lm(log(HR_INV) ~ log(HR_BICR), data = d, weights = d$Sample_size)
  return(summary(model)$r.squared)
}

boot_result <- boot(HRds, rsq_function, R = 1000)
R2_ci <- boot.ci(boot_result, type = "perc")$percent[4:5]  # 95% CI

HRds <- HRds %>%
  mutate(Label = case_when(
    Label %in% c("blinded", "double-blind") ~ "blinded",
    TRUE ~ Label
  ))

ggplot(HRds, aes(x = log(HR_INV), y = log(HR_BICR), size = Sample_size, color = Label)) +
  geom_point(alpha = 0.3) +
  geom_smooth(
    method = "lm", se = TRUE, color = "black", size = 0.5,
    inherit.aes = FALSE, aes(x = log(HR_INV), y = log(HR_BICR))
  ) +
  scale_size_continuous(name = "Sample Size", range = c(2, 10)) +
  scale_color_manual(values = c("open label" = "red", "blinded" = "blue")) +
  labs(
    title = "Correlation between HR_BICR and HR_INV",
    subtitle = paste0(
      "r = ", round(r_value, 3), 
      " (95% CI: ", round(r_ci[1], 3), "-", round(r_ci[2], 3), ")",
      "\nR² = ", round(R2_value, 3), 
      " (95% CI: ", round(R2_ci[1], 3), "-", round(R2_ci[2], 3), ")"
    ),
    x = "Log LE HR for Response",
    y = "Log BICR HR for Response",
    color = "Study Type"
  ) +
  theme_minimal() +
  guides(size = "none")

# Calculate the number of categories (There are 22 obs in OR, divided into 3 categories by label)

# Filter out data that are both 'blinded' and Phase 3
ds_blind_phase3 <- ds[ds$Blinded.open.label == "blinded" & ds$Study.Phase == "3", ]

# Then proceed with subsequent data processing
ORds_blind <- as.data.frame(cbind(
  ds_blind_phase3$PMID, ds_blind_phase3$Sample.size, 
  ds_blind_phase3$BICR.ORR.experimental, ds_blind_phase3$X95..CI.BICR.ORR.L.experimental, ds_blind_phase3$X95..CI.BICR.ORR.U.experimental,
  ds_blind_phase3$BICR.ORR.control, ds_blind_phase3$X95..CI.BICR.ORR.L.control, ds_blind_phase3$X95..CI.BICR.ORR.U.control,
  ds_blind_phase3$INV.ORR.experimental, ds_blind_phase3$X95..CI.INV.ORR.L.experimental, ds_blind_phase3$X95..CI.INV.ORR.U.experimental,
  ds_blind_phase3$INV.ORR.control, ds_blind_phase3$X95..CI.INV.ORR.L.control, ds_blind_phase3$X95..CI.INV.ORR.U.control
))

# Subsequent operations remain the same
for (i in 2:ncol(ORds_blind)) {
  ORds_blind[, i] <- as.numeric(ORds_blind[, i])
}

colnames(ORds_blind) <- c(
  "PMID", "Sample_size", 
  "OR_BICR_exp", "OR_BICR_L_exp", "OR_BICR_U_exp",
  "OR_BICR_c", "OR_BICR_L_c", "OR_BICR_U_c",
  "OR_INV_exp", "OR_INV_L_exp", "OR_INV_U_exp",
  "OR_INV_c", "OR_INV_L_c", "OR_INV_U_c"
)

ORds_blind$OR_BICR <- (ORds_blind$OR_BICR_exp * (1 - ORds_blind$OR_BICR_exp)) / 
  (ORds_blind$OR_BICR_c * (1 - ORds_blind$OR_BICR_c))
ORds_blind$OR_INV <- (ORds_blind$OR_INV_exp * (1 - ORds_blind$OR_INV_exp)) / 
  (ORds_blind$OR_INV_c * (1 - ORds_blind$OR_INV_c))

ORds_blind <- ORds_blind[!is.na(ORds_blind$OR_BICR) & !is.na(ORds_blind$OR_INV), ]

# Correlation test
cor_test_blind <- cor.test(log(ORds_blind$OR_BICR), log(ORds_blind$OR_INV), method = "pearson")
r_value_blind <- cor_test_blind$estimate    # correlation coefficient r
r_ci_blind <- cor_test_blind$conf.int       # 95% confidence interval

cat("Pearson r (blind phase 3 studies):", round(r_value_blind, 3), "95% CI:", round(r_ci_blind, 3), "\n")

# Linear model and R² calculation
lm_model_blind <- lm(log(ORds_blind$OR_INV) ~ log(ORds_blind$OR_BICR), weights = ORds_blind$Sample_size)
R2_value_blind <- summary(lm_model_blind)$r.squared

boot_result_blind <- boot(ORds_blind, rsq_function, R = 1000)
R2_ci_blind <- boot.ci(boot_result_blind, type = "perc")$percent[4:5]

cat("R² (blind phase 3 studies):", round(R2_value_blind, 3), "95% CI:", round(R2_ci_blind, 3), "\n")

# OddsRR calculation
ORds_blind$OddsRR <- log(ORds_blind$OR_BICR) / log(ORds_blind$OR_INV)
mean_OR_blind <- mean(ORds_blind$OddsRR, na.rm = TRUE)
ci_OR_blind <- t.test(ORds_blind$OddsRR)$conf.int

cat("OddsRR (blind phase 3 studies):", round(mean_OR_blind, 3), "95% CI:", round(ci_OR_blind, 3), "\n")

##################################################################
# Filter out data that are both 'double-blind' and Phase 3
ds_double_blind_phase3 <- ds[ds$Blinded.open.label == "double-blind" & ds$Study.Phase == "3", ]

# Then proceed with subsequent data processing
ORds_double_blind <- as.data.frame(cbind(
  ds_double_blind_phase3$PMID, ds_double_blind_phase3$Sample.size, 
  ds_double_blind_phase3$BICR.ORR.experimental, ds_double_blind_phase3$X95..CI.BICR.ORR.L.experimental, ds_double_blind_phase3$X95..CI.BICR.ORR.U.experimental,
  ds_double_blind_phase3$BICR.ORR.control, ds_double_blind_phase3$X95..CI.BICR.ORR.L.control, ds_double_blind_phase3$X95..CI.BICR.ORR.U.control,
  ds_double_blind_phase3$INV.ORR.experimental, ds_double_blind_phase3$X95..CI.INV.ORR.L.experimental, ds_double_blind_phase3$X95..CI.INV.ORR.U.experimental,
  ds_double_blind_phase3$INV.ORR.control, ds_double_blind_phase3$X95..CI.INV.ORR.L.control, ds_double_blind_phase3$X95..CI.INV.ORR.U.control
))

for (i in 2:ncol(ORds_double_blind)) {
  ORds_double_blind[, i] <- as.numeric(ORds_double_blind[, i])
}

colnames(ORds_double_blind) <- c(
  "PMID", "Sample_size", 
  "OR_BICR_exp", "OR_BICR_L_exp", "OR_BICR_U_exp",
  "OR_BICR_c", "OR_BICR_L_c", "OR_BICR_U_c",
  "OR_INV_exp", "OR_INV_L_exp", "OR_INV_U_exp",
  "OR_INV_c", "OR_INV_L_c", "OR_INV_U_c"
)

ORds_double_blind$OR_BICR <- (ORds_double_blind$OR_BICR_exp * (1 - ORds_double_blind$OR_BICR_exp)) /
  (ORds_double_blind$OR_BICR_c * (1 - ORds_double_blind$OR_BICR_c))
ORds_double_blind$OR_INV <- (ORds_double_blind$OR_INV_exp * (1 - ORds_double_blind$OR_INV_exp)) /
  (ORds_double_blind$OR_INV_c * (1 - ORds_double_blind$OR_INV_c))

ORds_double_blind <- ORds_double_blind[!is.na(ORds_double_blind$OR_BICR) & !is.na(ORds_double_blind$OR_INV), ]

# Correlation test
cor_test_double_blind <- cor.test(log(ORds_double_blind$OR_BICR), log(ORds_double_blind$OR_INV), method = "pearson")
r_value_double_blind <- cor_test_double_blind$estimate
r_ci_double_blind <- cor_test_double_blind$conf.int

cat("Pearson r (double-blind phase 3 studies):", round(r_value_double_blind, 3), "95% CI:", round(r_ci_double_blind, 3), "\n")

# Linear model and R² calculation
lm_model_double_blind <- lm(log(ORds_double_blind$OR_INV) ~ log(ORds_double_blind$OR_BICR),
                            weights = ORds_double_blind$Sample_size)
R2_value_double_blind <- summary(lm_model_double_blind)$r.squared

boot_result_double_blind <- boot(ORds_double_blind, rsq_function, R = 1000)
R2_ci_double_blind <- boot.ci(boot_result_double_blind, type = "perc")$percent[4:5]

cat("R² (double-blind phase 3 studies):", round(R2_value_double_blind, 3), "95% CI:", round(R2_ci_double_blind, 3), "\n")

# OddsRR calculation
ORds_double_blind$OddsRR <- log(ORds_double_blind$OR_BICR) / log(ORds_double_blind$OR_INV)
mean_OR_double_blind <- mean(ORds_double_blind$OddsRR, na.rm = TRUE)
ci_OR_double_blind <- t.test(ORds_double_blind$OddsRR)$conf.int

cat("OddsRR (double-blind phase 3 studies):", round(mean_OR_double_blind, 3), 
    "95% CI:", round(ci_OR_double_blind, 3), "\n")


##################################################################
# Filter out data that are both 'open label' and Phase 3
ds_open_label_phase3 <- ds[ds$Blinded.open.label == "open label" & ds$Study.Phase == "3", ]

# Then proceed with subsequent data processing
ORds_open_label <- as.data.frame(cbind(
  ds_open_label_phase3$PMID, ds_open_label_phase3$Sample.size, 
  ds_open_label_phase3$BICR.ORR.experimental, ds_open_label_phase3$X95..CI.BICR.ORR.L.experimental, ds_open_label_phase3$X95..CI.BICR.ORR.U.experimental,
  ds_open_label_phase3$BICR.ORR.control, ds_open_label_phase3$X95..CI.BICR.ORR.L.control, ds_open_label_phase3$X95..CI.BICR.ORR.U.control,
  ds_open_label_phase3$INV.ORR.experimental, ds_open_label_phase3$X95..CI.INV.ORR.L.experimental, ds_open_label_phase3$X95..CI.INV.ORR.U.experimental,
  ds_open_label_phase3$INV.ORR.control, ds_open_label_phase3$X95..CI.INV.ORR.L.control, ds_open_label_phase3$X95..CI.INV.ORR.U.control
))

for (i in 2:ncol(ORds_open_label)) {
  ORds_open_label[, i] <- as.numeric(ORds_open_label[, i])
}

colnames(ORds_open_label) <- c(
  "PMID", "Sample_size", 
  "OR_BICR_exp", "OR_BICR_L_exp", "OR_BICR_U_exp",
  "OR_BICR_c", "OR_BICR_L_c", "OR_BICR_U_c",
  "OR_INV_exp", "OR_INV_L_exp", "OR_INV_U_exp",
  "OR_INV_c", "OR_INV_L_c", "OR_INV_U_c"
)

ORds_open_label$OR_BICR <- (ORds_open_label$OR_BICR_exp * (1 - ORds_open_label$OR_BICR_exp)) /
  (ORds_open_label$OR_BICR_c * (1 - ORds_open_label$OR_BICR_c))
ORds_open_label$OR_INV <- (ORds_open_label$OR_INV_exp * (1 - ORds_open_label$OR_INV_exp)) /
  (ORds_open_label$OR_INV_c * (1 - ORds_open_label$OR_INV_c))

ORds_open_label <- ORds_open_label[!is.na(ORds_open_label$OR_BICR) & !is.na(ORds_open_label$OR_INV), ]

# Correlation test
cor_test_open_label <- cor.test(log(ORds_open_label$OR_BICR), log(ORds_open_label$OR_INV), method = "pearson")
r_value_open_label <- cor_test_open_label$estimate
r_ci_open_label <- cor_test_open_label$conf.int

cat("Pearson r (open-label phase 3 studies):", round(r_value_open_label, 3), 
    "95% CI:", round(r_ci_open_label, 3), "\n")

# Linear model and R² calculation
lm_model_open_label <- lm(log(ORds_open_label$OR_INV) ~ log(ORds_open_label$OR_BICR),
                          weights = ORds_open_label$Sample_size)
R2_value_open_label <- summary(lm_model_open_label)$r.squared

boot_result_open_label <- boot(ORds_open_label, rsq_function, R = 1000)
R2_ci_open_label <- boot.ci(boot_result_open_label, type = "perc")$percent[4:5]

cat("R² (open-label phase 3 studies):", round(R2_value_open_label, 3), 
    "95% CI:", round(R2_ci_open_label, 3), "\n")

# OddsRR calculation
ORds_open_label$OddsRR <- log(ORds_open_label$OR_BICR) / log(ORds_open_label$OR_INV)
mean_OR_open_label <- mean(ORds_open_label$OddsRR, na.rm = TRUE)
ci_OR_open_label <- t.test(ORds_open_label$OddsRR)$conf.int

cat("OddsRR (open-label phase 3 studies):", round(mean_OR_open_label, 3), 
    "95% CI:", round(ci_OR_open_label, 3), "\n")
