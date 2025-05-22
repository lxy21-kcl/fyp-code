library(readxl)
ds<- read_excel("data（RCT 2 and 3）.xlsx")


HRds <- as.data.frame(cbind(ds$PMID, ds$`Sample size`, ds$`HR for BICR-PFS`, ds$`95% CI BICR-HR L`, ds$`95% CI BICR-HR U`,
                            ds$`HR for INV-PFS`, ds$`95% CI INV L`, ds$`95% CI INV U`,ds$`follow up time`, ds$`Cancer Type`, ds$`Blinded/open label`,
                            ds$`Study Phase`))

print(ncol(HRds))

# 选择哪些列是数值列，哪些是字符列
num_cols <- 2:9  # 这些是数值列
char_cols <- c(1, 10, 11)  # 这些是字符列

# 转换数值列
for (i in num_cols) {
  HRds[, i] <- as.numeric(HRds[, i])
}

# 转换字符列
for (i in char_cols) {
  HRds[, i] <- as.character(HRds[, i])
}

# 设置列名
colnames(HRds) <- c("PMID", "Sample_size", "HR_BICR", "HR_BICR_L", "HR_BICR_U",
                    "HR_INV", "HR_INV_L", "HR_INV_U", "follow-up time","Cancer_type", "Label","Phase")

# 删除 HR_BICR 和 HR_INV 都是 NA 的行
HRds <- HRds[!is.na(HRds$HR_BICR) & !is.na(HRds$HR_INV),]

# 统计 Label 为 "open label" "blinded" 和"double-blind"的行数
open_label_count <- sum(HRds$Label == "open label", na.rm = TRUE)
open_label_count <- nrow(subset(HRds, Label == "open label"))
print(open_label_count)

blinded_count <- sum(HRds$Label == "blinded", na.rm = TRUE)
print(blinded_count)
blinded_count <- sum(HRds$Label == "double-blind", na.rm = TRUE)
print(blinded_count)

# 统计sample size
q1 <- quantile(HRds$Sample_size, probs = 0.25, na.rm = TRUE)
q3 <- quantile(HRds$Sample_size, probs = 0.75, na.rm = TRUE)

count_lt_Q1 <- sum(HRds$Sample_size < q1, na.rm = TRUE)
count_Q1_Q3 <- sum(HRds$Sample_size >= q1 & HRds$Sample_size <= q3, na.rm = TRUE)
count_gt_Q3 <- sum(HRds$Sample_size > q3, na.rm = TRUE)

print(paste("样本量 Q1（25%分位数）:", round(q1, 1)))
print(paste("样本量 Q3（75%分位数）:", round(q3, 1)))

print(paste("Sample size < Q1:", count_lt_Q1))
print(paste("Sample size Q1~Q3:", count_Q1_Q3))
print(paste("Sample size > Q3:", count_gt_Q3))

#or
count_lt_250 <- sum(HRds$Sample_size < 250, na.rm = TRUE)        
count_250_450 <- sum(HRds$Sample_size >= 250 & HRds$Sample_size <= 450, na.rm = TRUE)  
count_gt_450 <- sum(HRds$Sample_size > 450, na.rm = TRUE)         

print(paste("Sample size < 250:", count_lt_250))
print(paste("Sample size 250~450:", count_250_450))
print(paste("Sample size > 450:", count_gt_450))



#统计cancer type
count_LYM <- sum(HRds$Cancer_type == "LYM", na.rm = TRUE)
count_LEU <- sum(HRds$Cancer_type == "LEU", na.rm = TRUE)
count_MM <- sum(HRds$Cancer_type == "MM", na.rm = TRUE)

print(paste("LYM count:", count_LYM))
print(paste("LEU count:", count_LEU))
print(paste("MM count:", count_MM))

#统计follow-up duration
# 创建一个新列，将 follow-up time 分类为 Q1, Q2, Q3
HRds$followup_category <- cut(HRds$`follow-up time`, 
                              breaks = c(-Inf, 24, 36, Inf), 
                              labels = c("<24", "24-36", ">36"),
                              right = TRUE)
#or
HRds$followup_category <- cut(HRds$`follow-up time`, 
                              breaks = c(-Inf, 30, Inf), 
                              labels = c("<30", ">30"),
                              right = TRUE)
# 统计每个类别的数量
followup_counts <- table(HRds$followup_category)

# 计算百分比
followup_percentages <- prop.table(followup_counts) * 100

print(followup_counts)
print(followup_percentages)

# Get total number of rows
total_rows <- nrow(HRds)

# Count labels
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
total_valid <- sum(!is.na(HRds$Sample_size))
count_lt_Q1 <- sum(HRds$Sample_size < q1, na.rm = TRUE)
count_Q1_Q3 <- sum(HRds$Sample_size >= q1 & HRds$Sample_size <= q3, na.rm = TRUE)
count_gt_Q3 <- sum(HRds$Sample_size > q3, na.rm = TRUE)

percent_lt_Q1 <- (count_lt_Q1 / total_valid) * 100
percent_Q1_Q3 <- (count_Q1_Q3 / total_valid) * 100
percent_gt_Q3 <- (count_gt_Q3 / total_valid) * 100

# 输出统计结果（带百分比）
print(paste("Sample size < Q1:", count_lt_Q1, 
            "(", round(percent_lt_Q1, 1), "% )"))
print(paste("Sample size Q1~Q3:", count_Q1_Q3, 
            "(", round(percent_Q1_Q3, 1), "% )"))
print(paste("Sample size > Q3:", count_gt_Q3, 
            "(", round(percent_gt_Q3, 1), "% )"))

# Calculate percentages
percent_lt_250 <- (count_lt_250 / total_rows) * 100
percent_250_450 <- (count_250_450 / total_rows) * 100
percent_gt_450 <- (count_gt_450 / total_rows) * 100

# Print results
print(paste("Sample size < 250:", count_lt_250, "(", round(percent_lt_250, 2), "% )"))
print(paste("Sample size 250~450:", count_250_450, "(", round(percent_250_450, 2), "% )"))
print(paste("Sample size > 450:", count_gt_450, "(", round(percent_gt_450, 2), "% )"))



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
#ORR
library(readxl)
ds<- read_excel("data（RCT 2 and 3）.xlsx")

ORds<-as.data.frame(cbind(ds$PMID,ds$`Sample size`,ds$`BICR-ORR-experimental`,ds$`95% CI BICR ORR L-experimental`,
                          ds$`95% CI BICR ORR U-experimental`,ds$`BICR-ORR-control`,ds$`95% CI BICR ORR L-control`,ds$`95% CI BICR ORR U-control`,
                          ds$`INV-ORR-experimental`,ds$`95% CI INV ORR L-experimental`,ds$`95% CI INV ORR U-experimental`,
                          ds$`INV-ORR-control`,ds$`95% CI INV ORR L-control`,
                          ds$`95% CI INV ORR U-control`, ds$`follow up time`))

ORds$Blinded.open.label <- ds$`Blinded/open label`
ORds$Cancer.Type <- ds$`Cancer Type`
ORds$Phase <- ds$`Study Phase`
# 选择哪些列是数值列，哪些是字符列
num_cols <- 2:15  # 这些是数值列
char_cols <- c(1, 16, 17)  # 这些是字符列

# 转换数值列
for (i in num_cols) {
  ORds[, i] <- as.numeric(ORds[, i])
}

# 转换字符列
for (i in char_cols) {
  ORds[, i] <- as.character(ORds[, i])
}

colnames(ORds)<-c("PMID","Sample_size","OR_BICR_exp","OR_BICR_L_exp","OR_BICR_U_exp",
                  "OR_BICR_c","OR_BICR_L_c","OR_BICR_U_c",
                  "OR_INV_exp","OR_INV_L_exp","OR_INV_U_exp",
                  "OR_INV_c","OR_INV_L_c","OR_INV_U_c","follow-up time","Label","Cancer_type","Phase")

ORds$OR_BICR<-(ORds$OR_BICR_exp*(1-ORds$OR_BICR_exp))/(ORds$OR_BICR_c*(1-ORds$OR_BICR_c))
ORds$OR_INV<-(ORds$OR_INV_exp*(1-ORds$OR_INV_exp))/(ORds$OR_INV_c*(1-ORds$OR_INV_c))

ORds<-ORds[!is.na(ORds$OR_BICR)&!is.na(ORds$OR_INV),]        ###计算phase 2时清洗数据有问题

total_rows <- nrow(ORds)
# 统计 Label 为 "open label" "blinded" 和"double-blind"的行数
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


# 统计Sample size counts
q1 <- quantile(ORds$Sample_size, probs = 0.25, na.rm = TRUE)
q3 <- quantile(ORds$Sample_size, probs = 0.75, na.rm = TRUE)

total_valid <- sum(!is.na(ORds$Sample_size))

# 统计区间数量
count_lt_Q1 <- sum(ORds$Sample_size < q1, na.rm = TRUE)
count_Q1_Q3 <- sum(ORds$Sample_size >= q1 & ORds$Sample_size <= q3, na.rm = TRUE)
count_gt_Q3 <- sum(ORds$Sample_size > q3, na.rm = TRUE)

# 计算百分比
percent_lt_Q1 <- (count_lt_Q1 / total_valid) * 100
percent_Q1_Q3 <- (count_Q1_Q3 / total_valid) * 100
percent_gt_Q3 <- (count_gt_Q3 / total_valid) * 100

# 输出统计结果（带百分比）
print(paste("Sample size < Q1:", count_lt_Q1, 
            "(", round(percent_lt_Q1, 1), "% )"))
print(paste("Sample size Q1~Q3:", count_Q1_Q3, 
            "(", round(percent_Q1_Q3, 1), "% )"))
print(paste("Sample size > Q3:", count_gt_Q3, 
            "(", round(percent_gt_Q3, 1), "% )"))

# or
count_lt_250 <- sum(ORds$Sample_size < 250, na.rm = TRUE)
count_250_450 <- sum(ORds$Sample_size >= 250 & ORds$Sample_size <= 450, na.rm = TRUE)
count_gt_450 <- sum(ORds$Sample_size > 450, na.rm = TRUE)

# Calculate percentages
percent_lt_250 <- (count_lt_250 / total_rows) * 100
percent_250_450 <- (count_250_450 / total_rows) * 100
percent_gt_450 <- (count_gt_450 / total_rows) * 100

# Print results
print(paste("Sample size < 250:", count_lt_250, "(", round(percent_lt_250, 2), "% )"))
print(paste("Sample size 250~450:", count_250_450, "(", round(percent_250_450, 2), "% )"))
print(paste("Sample size > 450:", count_gt_450, "(", round(percent_gt_450, 2), "% )"))



# 统计Cancer type counts
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



#统计follow-up duration
# 创建一个新列，将 follow-up time 分类为 Q1, Q2, Q3
ORds$followup_category <- cut(ORds$`follow-up time`, 
                              breaks = c(-Inf, 24, 36, Inf), 
                              labels = c("<24", "24-36", ">36"),
                              right = TRUE)
#or
ORds$followup_category <- cut(ORds$`follow-up time`, 
                              breaks = c(-Inf, 30, Inf), 
                              labels = c("<30", ">30"),
                              right = TRUE)

# 统计每个类别的数量
followup_counts <- table(ORds$followup_category)

# 计算百分比
followup_percentages <- prop.table(followup_counts) * 100

print(followup_counts)
print(followup_percentages)


###########################################################################
# 计算HRR
library(lme4)

# 确保数值列转换为数值型
HRds$HR_BICR <- as.numeric(HRds$HR_BICR)
HRds$HR_INV <- as.numeric(HRds$HR_INV)
HRds$Sample_size <- as.numeric(HRds$Sample_size)

HRds$HRR <- ifelse(!is.na(HRds$HR_BICR) & !is.na(HRds$HR_INV) & HRds$HR_INV != 0, 
                   HRds$HR_INV / HRds$HR_BICR, NA)

# 定义 HRR 的区间
HRds$HRR_bin <- cut(HRds$HRR, 
                    breaks = c(-Inf, 0.85, 1, 1.15, Inf), 
                    labels = c("≤0.85", "(0.85, 1]", "(1, 1.15]", ">1.15"), 
                    include.lowest = TRUE)
# 统计每个区间的数量
HRR_counts <- table(HRds$HRR_bin)

# 计算百分比
HRR_percentages <- round((HRR_counts / sum(HRR_counts)) * 100, 2)

# 打印结果
print(HRR_counts)
print(HRR_percentages)

# 定义 Sample_size 的区间
q1 <- quantile(HRds$Sample_size, probs = 0.25, na.rm = TRUE)
q3 <- quantile(HRds$Sample_size, probs = 0.75, na.rm = TRUE)

# 定义基于分位数的区间
sample_size_bins <- cut(HRds$Sample_size,
                        breaks = c(-Inf, q1, q3, Inf),
                        labels = c(
                          paste0("<Q1(", round(q1,1), ")"), 
                          paste0("Q1-Q3(", round(q1,1), "-", round(q3,1), ")"),
                          paste0(">Q3(", round(q3,1), ")")
                        ),
                        include.lowest = TRUE)

# 计算交叉统计
HRR_sample_counts <- table(sample_size_bins, HRds$HRR_bin)
HRR_sample_percentages <- round(prop.table(HRR_sample_counts, margin = 1) * 100, 2)

# 打印统计结果
cat("样本量区间 vs HRR分布 - 频数表:\n")
print(HRR_sample_counts)

cat("\n样本量区间 vs HRR分布 - 行百分比:\n")
print(HRR_sample_percentages)

#or
sample_size_bins <- cut(HRds$Sample_size, 
                        breaks = c(-Inf, 250, 450, Inf), 
                        labels = c("<250", "250-450", ">450"), 
                        include.lowest = TRUE)

# 计算每个 Sample_size 区间内的 HRR_bin 统计数
HRR_sample_counts <- table(sample_size_bins, HRds$HRR_bin)
HRR_sample_percentages <- round(prop.table(HRR_sample_counts, margin = 1) * 100, 2)
print(HRR_sample_counts)
print(HRR_sample_percentages)


#定义masking
HRR_label_sample_counts <- table(HRds$Label, HRds$HRR_bin)
print(HRR_label_sample_counts)
HRR_label_percentages <- round(prop.table(HRR_label_sample_counts, margin = 1) * 100, 2)
print(HRR_label_percentages)
#定义masking
HRR_label_sample_counts <- table(HRds$Label, HRds$HRR_bin)
print(HRR_label_sample_counts)
HRR_label_percentages <- round(prop.table(HRR_label_sample_counts, margin = 1) * 100, 2)
print(HRR_label_percentages)


#定义Cancer type
HRR_cancer_type_counts <- table(HRds$Cancer_type,HRds$HRR_bin)
print(HRR_cancer_type_counts)
HRR_cancer_percentages <- round(prop.table(HRR_cancer_type_counts, margin = 1) * 100, 2)
print(HRR_cancer_percentages)


#定义follow-up time
HRR_followup_counts <- table(HRds$followup_category, HRds$HRR_bin)
HRR_followup_percent <- prop.table(HRR_followup_counts, margin = 1) * 100
print(round(HRR_followup_counts, 2))
print(round(HRR_followup_percent, 2))



########################################################################
# Odds Ratio ratio 计算
ORds$OR_BICR <- as.numeric(ORds$OR_BICR)
ORds$OR_INV <- as.numeric(ORds$OR_INV)
ORds$Sample_size <- as.numeric(ORds$Sample_size)

ORds$ORR <- ifelse(!is.na(ORds$OR_BICR) & !is.na(ORds$OR_INV) & ORds$OR_INV != 0, 
                   ORds$OR_INV / ORds$OR_BICR, NA)

ORds$ORR_bin <- cut(ORds$ORR, 
                    breaks = c(-Inf, 0.85, 1, 1.15, Inf), 
                    labels = c("≤0.85", "(0.85, 1]", "(1, 1.15]", ">1.15"), 
                    include.lowest = TRUE)

# 统计每个区间的数量
ORR_counts <- table(ORds$ORR_bin)

# 计算百分比
ORR_percentages <- round((ORR_counts / sum(ORR_counts)) * 100, 2)

# 打印结果
print(ORR_counts)
print(ORR_percentages)

# 定义 Sample_size 的区间
sample_size_bins <- cut(ORds$Sample_size, 
                        breaks = c(-Inf, 100, 500, Inf), 
                        labels = c("<100", "100-500", ">500"), 
                        include.lowest = TRUE)

# 计算每个 Sample_size 区间内的 ORR_bin 统计数
ORR_sample_counts <- table(sample_size_bins, ORds$ORR_bin)
ORR_sample_percentages <- round(prop.table(ORR_sample_counts, margin = 1) * 100, 2)
print(ORR_sample_counts)
print(ORR_sample_percentages)

# 定义 masking
ORR_label_sample_counts <- table(ORds$Label, ORds$ORR_bin)
print(ORR_label_sample_counts)
ORR_label_percentages <- round(prop.table(ORR_label_sample_counts, margin = 1) * 100, 2)
print(ORR_label_percentages)

# 定义 Cancer type
ORR_cancer_type_counts <- table(ORds$Cancer_type, ORds$ORR_bin)
print(ORR_cancer_type_counts)
ORR_cancer_percentages <- round(prop.table(ORR_cancer_type_counts, margin = 1) * 100, 2)
print(ORR_cancer_percentages)

# 定义 follow-up time
ORR_followup_counts <- table(ORds$followup_category, ORds$ORR_bin)
ORR_followup_percent <- prop.table(ORR_followup_counts, margin = 1) * 100
print(round(ORR_followup_counts, 2))
print(round(ORR_followup_percent, 2))
