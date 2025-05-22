# Read the CSV file, specifying encoding
ds <- read.csv("datainfo_transpose_1.20.csv", stringsAsFactors = FALSE, fileEncoding = "latin1")

# Subset for Phase 3 OR Phase 2 (Note: The second line overwrites the first)
ds <- ds[ds$Study.Phase == "3", ]  # or
ds <- ds[ds$Study.Phase == "2", ]

# Create a data frame for HR-related data
HRds <- as.data.frame(cbind(
  ds$ï..PMID, ds$Sample.size, ds$HR.for.BICR.PFS, ds$X95..CI.BICR.HR.L, ds$X95..CI.BICR.HR.U,
  ds$HR.for.INV.PFS, ds$X95..CI.INV.L, ds$X95..CI.INV.U, ds$follow.up.time, 
  ds$Cancer.Type, ds$Blinded.open.label, ds$Study.Phase
))

print(ncol(HRds))

# Identify which columns are numeric and which are character
num_cols <- 2:9   # These are numeric columns
char_cols <- c(1, 10, 11)  # These are character columns

# Convert numeric columns
for (i in num_cols) {
  HRds[, i] <- as.numeric(HRds[, i])
}

# Convert character columns
for (i in char_cols) {
  HRds[, i] <- as.character(HRds[, i])
}

# Set column names
colnames(HRds) <- c(
  "PMID", "Sample_size", "HR_BICR", "HR_BICR_L", "HR_BICR_U",
  "HR_INV", "HR_INV_L", "HR_INV_U", "follow-up time", "Cancer_type",
  "Label", "Phase"
)

# Remove rows where both HR_BICR and HR_INV are NA
HRds <- HRds[!is.na(HRds$HR_BICR) & !is.na(HRds$HR_INV), ]

# Count the number of rows where Label is "open label", "blinded", and "double-blind"
open_label_count <- sum(HRds$Label == "open label", na.rm = TRUE)
open_label_count <- nrow(subset(HRds, Label == "open label"))
print(open_label_count)

blinded_count <- sum(HRds$Label == "blinded", na.rm = TRUE)
print(blinded_count)
blinded_count <- sum(HRds$Label == "double-blind", na.rm = TRUE)
print(blinded_count)

# Count the sample size categories
count_lt_100 <- sum(HRds$Sample_size < 100, na.rm = TRUE)                  # Sample_size < 100
count_100_500 <- sum(HRds$Sample_size >= 100 & HRds$Sample_size <= 500, na.rm = TRUE)  # 100 <= Sample_size <= 500
count_gt_500 <- sum(HRds$Sample_size > 500, na.rm = TRUE)                  # Sample_size > 500

print(paste("Sample size < 100:", count_lt_100))
print(paste("Sample size 100-500:", count_100_500))
print(paste("Sample size > 500:", count_gt_500))

# Count cancer types
count_LYM <- sum(HRds$Cancer_type == "LYM", na.rm = TRUE)
count_LEU <- sum(HRds$Cancer_type == "LEU", na.rm = TRUE)
count_MM <- sum(HRds$Cancer_type == "MM", na.rm = TRUE)

print(paste("LYM count:", count_LYM))
print(paste("LEU count:", count_LEU))
print(paste("MM count:", count_MM))

# Count the follow-up duration
# Create a new column that categorizes the follow-up time into Q1, Q2, Q3
HRds$followup_category <- cut(
  HRds$`follow-up time`,
  breaks = c(-Inf, 24, 36, Inf),
  labels = c("<24", "24-36", ">36"),
  right = TRUE
)
# or
HRds$followup_category <- cut(
  HRds$`follow-up time`,
  breaks = c(-Inf, 30, Inf),
  labels = c("<30", ">30"),
  right = TRUE
)

# Count how many entries fall into each category
followup_counts <- table(HRds$followup_category)

# Compute percentages
followup_percentages <- prop.table(followup_counts) * 100

print(followup_counts)
print(followup_percentages)

# Get the total number of rows
total_rows <- nrow(HRds)

# Count each Label
open_label_count <- sum(HRds$Label == "open label", na.rm = TRUE)
blinded_count <- sum(HRds$Label == "blinded", na.rm = TRUE)
double_blind_count <- sum(HRds$Label == "double-blind", na.rm = TRUE)

# Calculate percentages
open_label_percent <- (open_label_count / total_rows) * 100
blinded_percent <- (blinded_count / total_rows) * 100
double_blind_percent <- (double_blind_count / total_rows) * 100

# Print results
print(paste("Open label count:", open_label_count, "(", round(open_label_percent, 2), "% )"))
print(paste("Blinded count:", blinded_count, "(", round(blinded_percent, 2), "% )"))
print(paste("Double-blind count:", double_blind_count, "(", round(double_blind_percent, 2), "% )"))

# Sample size counts
count_lt_100 <- sum(HRds$Sample_size < 100, na.rm = TRUE)
count_100_500 <- sum(HRds$Sample_size >= 100 & HRds$Sample_size <= 500, na.rm = TRUE)
count_gt_500 <- sum(HRds$Sample_size > 500, na.rm = TRUE)

# Calculate percentages
percent_lt_100 <- (count_lt_100 / total_rows) * 100
percent_100_500 <- (count_100_500 / total_rows) * 100
percent_gt_500 <- (count_gt_500 / total_rows) * 100

# Print results
print(paste("Sample size < 100:", count_lt_100, "(", round(percent_lt_100, 2), "% )"))
print(paste("Sample size 100-500:", count_100_500, "(", round(percent_100_500, 2), "% )"))
print(paste("Sample size > 500:", count_gt_500, "(", round(percent_gt_500, 2), "% )"))

# Cancer type counts
count_LYM <- sum(HRds$Cancer_type == "LYM", na.rm = TRUE)
count_LEU <- sum(HRds$Cancer_type == "LEU", na.rm = TRUE)
count_MM <- sum(HRds$Cancer_type == "MM", na.rm = TRUE)

# Calculate percentages
percent_LYM <- (count_LYM / total_rows) * 100
percent_LEU <- (count_LEU / total_rows) * 100
percent_MM <- (count_MM / total_rows) * 100

# Print results
print(paste("LYM count:", count_LYM, "(", round(percent_LYM, 2), "% )"))
print(paste("LEU count:", count_LEU, "(", round(percent_LEU, 2), "% )"))
print(paste("MM count:", count_MM, "(", round(percent_MM, 2), "% )"))

########################################################################
# ORR
ds <- read.csv("datainfo_transpose_1.20.csv", stringsAsFactors = FALSE, fileEncoding = "latin1")
ds <- ds[ds$Study.Phase == "2", ]

ORds <- as.data.frame(cbind(
  ds$ï..PMID, ds$Sample.size, ds$BICR.ORR.experimental, ds$X95..CI.BICR.ORR.L.experimental,
  ds$X95..CI.BICR.ORR.U.experimental, ds$BICR.ORR.control, ds$X95..CI.BICR.ORR.L.control, ds$X95..CI.BICR.ORR.U.control,
  ds$INV.ORR.experimental, ds$X95..CI.INV.ORR.L.experimental, ds$X95..CI.INV.ORR.U.experimental,
  ds$INV.ORR.control, ds$X95..CI.INV.ORR.L.control, ds$X95..CI.INV.ORR.U.control, ds$follow.up.time
))

# Add columns for Label, Cancer.Type, and Phase
ORds$Blinded.open.label <- ds$Blinded.open.label
ORds$Cancer.Type <- ds$Cancer.Type
ORds$Phase <- ds$Study.Phase

# Identify numeric and character columns
num_cols <- 2:15  # These are numeric columns
char_cols <- c(1, 16, 17)  # These are character columns

# Convert numeric columns
for (i in num_cols) {
  ORds[, i] <- as.numeric(ORds[, i])
}

# Convert character columns
for (i in char_cols) {
  ORds[, i] <- as.character(ORds[, i])
}

# Rename columns
colnames(ORds) <- c(
  "PMID", "Sample_size", "OR_BICR_exp", "OR_BICR_L_exp", "OR_BICR_U_exp",
  "OR_BICR_c", "OR_BICR_L_c", "OR_BICR_U_c",
  "OR_INV_exp", "OR_INV_L_exp", "OR_INV_U_exp",
  "OR_INV_c", "OR_INV_L_c", "OR_INV_U_c", "follow-up time",
  "Label", "Cancer_type", "Phase"
)

# Calculate BICR OR and INV OR
ORds$OR_BICR <- (ORds$OR_BICR_exp * (1 - ORds$OR_BICR_exp)) / (ORds$OR_BICR_c * (1 - ORds$OR_BICR_c))
ORds$OR_INV <- (ORds$OR_INV_exp * (1 - ORds$OR_INV_exp)) / (ORds$OR_INV_c * (1 - ORds$OR_INV_c))

# Remove rows with NA for OR_BICR or OR_INV (note: there might be an issue with data cleaning for Phase 2)
ORds <- ORds[!is.na(ORds$OR_BICR) & !is.na(ORds$OR_INV), ]

total_rows <- nrow(ORds)

# Count the number of rows where Label is "open label", "blinded", and "double-blind"
open_label_count <- sum(ORds$Label == "open label", na.rm = TRUE)
blinded_count <- sum(ORds$Label == "blinded", na.rm = TRUE)
double_blind_count <- sum(ORds$Label == "double-blind", na.rm = TRUE)

# Calculate percentages
open_label_percent <- (open_label_count / total_rows) * 100
blinded_percent <- (blinded_count / total_rows) * 100
double_blind_percent <- (double_blind_count / total_rows) * 100

# Print results
print(paste("Open label count:", open_label_count, "(", round(open_label_percent, 2), "% )"))
print(paste("Blinded count:", blinded_count, "(", round(blinded_percent, 2), "% )"))
print(paste("Double-blind count:", double_blind_count, "(", round(double_blind_percent, 2), "% )"))

# Sample size counts
count_lt_100 <- sum(ORds$Sample_size < 100, na.rm = TRUE)
count_100_500 <- sum(ORds$Sample_size >= 100 & ORds$Sample_size <= 500, na.rm = TRUE)
count_gt_500 <- sum(ORds$Sample_size > 500, na.rm = TRUE)

# Calculate percentages
percent_lt_100 <- (count_lt_100 / total_rows) * 100
percent_100_500 <- (count_100_500 / total_rows) * 100
percent_gt_500 <- (count_gt_500 / total_rows) * 100

# Print results
print(paste("Sample size < 100:", count_lt_100, "(", round(percent_lt_100, 2), "% )"))
print(paste("Sample size 100-500:", count_100_500, "(", round(percent_100_500, 2), "% )"))
print(paste("Sample size > 500:", count_gt_500, "(", round(percent_gt_500, 2), "% )"))

# Cancer type counts
count_LYM <- sum(ORds$Cancer_type == "LYM", na.rm = TRUE)
count_LEU <- sum(ORds$Cancer_type == "LEU", na.rm = TRUE)
count_MM <- sum(ORds$Cancer_type == "MM", na.rm = TRUE)

# Calculate percentages
percent_LYM <- (count_LYM / total_rows) * 100
percent_LEU <- (count_LEU / total_rows) * 100
percent_MM <- (count_MM / total_rows) * 100

# Print results
print(paste("LYM count:", count_LYM, "(", round(percent_LYM, 2), "% )"))
print(paste("LEU count:", count_LEU, "(", round(percent_LEU, 2), "% )"))
print(paste("MM count:", count_MM, "(", round(percent_MM, 2), "% )"))

# Count the follow-up duration
# Create a new column that categorizes the follow-up time
ORds$followup_category <- cut(
  ORds$`follow-up time`,
  breaks = c(-Inf, 24, 36, Inf),
  labels = c("<24", "24-36", ">36"),
  right = TRUE
)
# or
ORds$followup_category <- cut(
  ORds$`follow-up time`,
  breaks = c(-Inf, 30, Inf),
  labels = c("<30", ">30"),
  right = TRUE
)

# Count how many entries fall into each category
followup_counts <- table(ORds$followup_category)

# Compute percentages
followup_percentages <- prop.table(followup_counts) * 100

print(followup_counts)
print(followup_percentages)

###########################################################################
# Calculate HRR
library(lme4)

# Ensure numeric columns are indeed numeric
HRds$HR_BICR <- as.numeric(HRds$HR_BICR)
HRds$HR_INV <- as.numeric(HRds$HR_INV)
HRds$Sample_size <- as.numeric(HRds$Sample_size)

# Define HRR only if both HR_BICR and HR_INV are non-NA and HR_INV != 0
HRds$HRR <- ifelse(
  !is.na(HRds$HR_BICR) & !is.na(HRds$HR_INV) & HRds$HR_INV != 0,
  HRds$HR_INV / HRds$HR_BICR, 
  NA
)

# Define the HRR intervals
HRds$HRR_bin <- cut(
  HRds$HRR,
  breaks = c(-Inf, 0.85, 1, 1.15, Inf),
  labels = c("≤0.85", "(0.85, 1]", "(1, 1.15]", ">1.15"),
  include.lowest = TRUE
)

# Count each interval
HRR_counts <- table(HRds$HRR_bin)

# Compute percentages
HRR_percentages <- round((HRR_counts / sum(HRR_counts)) * 100, 2)

# Print results
print(HRR_counts)
print(HRR_percentages)

# Define intervals for sample size
sample_size_bins <- cut(
  HRds$Sample_size,
  breaks = c(-Inf, 100, 500, Inf),
  labels = c("<100", "100-500", ">500"),
  include.lowest = TRUE
)

# Compute how many HRR_bin fall into each sample size interval
HRR_sample_counts <- table(sample_size_bins, HRds$HRR_bin)
HRR_sample_percentages <- round(prop.table(HRR_sample_counts, margin = 1) * 100, 2)
print(HRR_sample_counts)
print(HRR_sample_percentages)

# Define masking (label)
HRR_label_sample_counts <- table(HRds$Label, HRds$HRR_bin)
print(HRR_label_sample_counts)
HRR_label_percentages <- round(prop.table(HRR_label_sample_counts, margin = 1) * 100, 2)
print(HRR_label_percentages)

# Define cancer type
HRR_cancer_type_counts <- table(HRds$Cancer_type, HRds$HRR_bin)
print(HRR_cancer_type_counts)
HRR_cancer_percentages <- round(prop.table(HRR_cancer_type_counts, margin = 1) * 100, 2)
print(HRR_cancer_percentages)

# Define follow-up time
HRR_followup_counts <- table(HRds$followup_category, HRds$HRR_bin)
HRR_followup_percent <- prop.table(HRR_followup_counts, margin = 1) * 100
print(round(HRR_followup_counts, 2))
print(round(HRR_followup_percent, 2))

########################################################################
# Calculate the Odds Ratio ratio (ORR)
ORds$OR_BICR <- as.numeric(ORds$OR_BICR)
ORds$OR_INV <- as.numeric(ORds$OR_INV)
ORds$Sample_size <- as.numeric(ORds$Sample_size)

# Define ORR only if both OR_BICR and OR_INV are non-NA and OR_INV != 0
ORds$ORR <- ifelse(
  !is.na(ORds$OR_BICR) & !is.na(ORds$OR_INV) & ORds$OR_INV != 0, 
  ORds$OR_INV / ORds$OR_BICR, 
  NA
)

# Define intervals for ORR
ORds$ORR_bin <- cut(
  ORds$ORR,
  breaks = c(-Inf, 0.85, 1, 1.15, Inf),
  labels = c("≤0.85", "(0.85, 1]", "(1, 1.15]", ">1.15"),
  include.lowest = TRUE
)

# Count each interval
ORR_counts <- table(ORds$ORR_bin)

# Compute percentages
ORR_percentages <- round((ORR_counts / sum(ORR_counts)) * 100, 2)

# Print results
print(ORR_counts)
print(ORR_percentages)

# Define intervals for sample size
sample_size_bins <- cut(
  ORds$Sample_size,
  breaks = c(-Inf, 100, 500, Inf),
  labels = c("<100", "100-500", ">500"),
  include.lowest = TRUE
)

# Compute how many ORR_bin fall into each sample size interval
ORR_sample_counts <- table(sample_size_bins, ORds$ORR_bin)
ORR_sample_percentages <- round(prop.table(ORR_sample_counts, margin = 1) * 100, 2)
print(ORR_sample_counts)
print(ORR_sample_percentages)

# Define masking (label)
ORR_label_sample_counts <- table(ORds$Label, ORds$ORR_bin)
print(ORR_label_sample_counts)
ORR_label_percentages <- round(prop.table(ORR_label_sample_counts, margin = 1) * 100, 2)
print(ORR_label_percentages)

# Define cancer type
ORR_cancer_type_counts <- table(ORds$Cancer_type, ORds$ORR_bin)
print(ORR_cancer_type_counts)
ORR_cancer_percentages <- round(prop.table(ORR_cancer_type_counts, margin = 1) * 100, 2)
print(ORR_cancer_percentages)

# Define follow-up time
ORR_followup_counts <- table(ORds$followup_category, ORds$ORR_bin)
ORR_followup_percent <- prop.table(ORR_followup_counts, margin = 1) * 100
print(round(ORR_followup_counts, 2))
print(round(ORR_followup_percent, 2))
