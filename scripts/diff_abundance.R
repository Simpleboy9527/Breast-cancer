# 数据标准化 (Data Normalization)  
normalized_meta_data <- meta_data_WGS_tumor %>%   
  mutate(age = scale(age))  

# 线性混合效应模型分析 (Linear Mixed Effects Model Analysis)  
# 设置 PAM50 水平 (Set PAM50 Levels)  
PAM50_levels <- c("Basal", "Her2", "LumA", "LumB", "Normal")  
# 初始化结果数据框 (Initialize Result Data Frame)  
linear_model_summary <- data.frame()  

micro_data_tumor_tss <- norm_tss(tse_tumor)@otu_table %>% as.matrix() %>% t() %>% as.data.frame()

# 循环遍历每个特征 (Loop Through Each Feature)  
for (feature in colnames(micro_data_tumor_tss)) {  
  # 计算流行率 (Calculate Prevalence)  
  prevalence <- sum(micro_data_tumor_tss[, feature] > 0) / nrow(micro_data_tumor_tss)  
  
  # 仅在流行率大于 0.1 时继续 (Continue Only if Prevalence > 0.1)  
  if (prevalence <= 0.1) {  
    next  
  }  
  
  # 创建临时数据框 (Create Temporary Data Frame)  
  temp_data <- normalized_meta_data %>%  
    mutate(response_variable = micro_data_tumor_tss[match(normalized_meta_data$sample, rownames(micro_data_tumor_tss)), feature])  
  
  # 确保 PAM50 是一个因子 (Ensure PAM50 is a Factor)  
  temp_data$PAM50 <- factor(temp_data$PAM50, levels = PAM50_levels)  
  
  # 检查 temp_data 是否为空 (Check if temp_data is Empty)  
  if (nrow(temp_data) == 0) {  
    next  
  }  
  
  # 拟合混合线性模型 (Fit Mixed Linear Model)  
  model <- lmerTest::lmer(response_variable ~ age + PAM50 + age:PAM50 + (1 | batch), data = temp_data)  
  
  # 获取模型的摘要 (Get Model Summary)  
  summary_model <- summary(model)  
  
  # 提取 p 值 (Extract p-values)  
  p_values <- summary_model$coefficients[, "Pr(>|t|)"]  
  
  # 将结果添加到 linear_model_summary 数据框 (Add Results to linear_model_summary Data Frame)  
  linear_model_summary <- rbind(linear_model_summary, c(taxa = feature, p_values))  
}  

# 设置列名 (Set Column Names)  
colnames(linear_model_summary) <- c("taxa", names(summary_model$coefficients[, "Pr(>|t|)"]))   

# 查看结果 (View Results)  
print(linear_model_summary)  

write.csv(linear_model_summary, "results/linear_model_summary_tumor.csv")  

# 对这些列进行 p 值调整 (Adjust p-values for these Columns)  
linear_model_summary <- linear_model_summary %>%  
  mutate(across(-c(1, 2), ~ p.adjust(., method = "BH"), .names = "q_{col}"))   

# 根据 q 值进行筛选，获取 q < 0.05 的行 (Filter by q-value to Get Rows with q < 0.05)  
q_value_columns <- grep("^q_", colnames(linear_model_summary), value = TRUE)  
filtered_linear_model_summary <- linear_model_summary %>%  
  filter(if_any(all_of(q_value_columns), ~ . < 0.05))  

# 热图 - 所有差异 (Heatmap - All Differences)  
filtered_micro_data <- micro_data_tumor_tss[, colnames(micro_data_tumor_tss) %in% filtered_linear_model_summary$taxa] %>%   
  t() %>% as.data.frame()  

annotation_columns <- meta_data_WGS_tumor[order(meta_data_WGS_tumor$age), c("age", "PAM50", "batch", "age_group_55", "age_group_45_60","sample","LNM.0.No.1.Yes.")]  
annotation_columns <- annotation_columns[!is.na(annotation_columns$age),]  
rownames(annotation_columns) <- annotation_columns$sample  
annotation_columns <- annotation_columns[, c("age", "age_group_55", "age_group_45_60","PAM50", "batch","LNM.0.No.1.Yes.")]  
filtered_micro_data <- filtered_micro_data[, rownames(annotation_columns)]  
rownames(filtered_micro_data) <- gsub(".*\\|g__", "g__", rownames(filtered_micro_data))  
write.csv(filtered_micro_data, "results/diff_micro_TT.csv")  

all_heatmap <- pheatmap(filtered_micro_data, cluster_rows = TRUE, cluster_cols = FALSE,   
                        annotation_col = annotation_columns, show_colnames = FALSE, color = my_color_palette(100))  
ggsave("results/figures/all_diff.png", plot = all_heatmap, width = 8, height = 5, units = "in", dpi = 300)  

# 热图 - 在年龄变量中 (Heatmap - By Age Variable)  
filtered_micro_data_age <- micro_data_tumor_tss[, colnames(micro_data_tumor_tss) %in% filtered_linear_model_summary$taxa[filtered_linear_model_summary$q_age <= 0.05]] %>%   
  t() %>% as.data.frame()  

filtered_micro_data_age <- filtered_micro_data_age[, rownames(annotation_columns)]   
rownames(filtered_micro_data_age) <- gsub(".*\\|g__", "g__", rownames(filtered_micro_data_age))  

age_heatmap <- pheatmap(filtered_micro_data_age, cluster_rows = TRUE, cluster_cols = FALSE,   
                        annotation_col = annotation_columns, show_colnames = FALSE, color = my_color_palette(100))  
ggsave("results/figures/age_diff.png", plot = age_heatmap, width = 8, height = 5, units = "in", dpi = 300)  

# 热图 - 在 PAM50 变量中 (Heatmap - By PAM50 Variable)  
filtered_micro_data_pam50 <- micro_data_tumor_tss[, colnames(micro_data_tumor_tss) %in% filtered_linear_model_summary$taxa[!filtered_linear_model_summary$q_age <= 0.05]] %>%   
  t() %>% as.data.frame()  
annotation_columns <- annotation_columns[order(annotation_columns$PAM50),]  

filtered_micro_data_pam50 <- filtered_micro_data_pam50[, rownames(annotation_columns)]  
rownames(filtered_micro_data_pam50) <- gsub(".*\\|g__", "g__", rownames(filtered_micro_data_pam50))  

pam_heatmap <- pheatmap(filtered_micro_data_pam50, cluster_rows = TRUE, cluster_cols = FALSE,   
                        annotation_col = annotation_columns, show_colnames = FALSE, color = my_color_palette(100))  
ggsave("results/figures/pam_diff.png", plot = pam_heatmap, width = 8, height = 5, units = "in", dpi = 300)  

# log10转换 (Log10 Transformation)  
filtered_micro_data_pam50_log10 <- log10(filtered_micro_data_pam50 + min(filtered_micro_data_pam50[filtered_micro_data_pam50 != 0]))  
pam50_heatmap_log10 <- pheatmap(filtered_micro_data_pam50_log10, cluster_rows = TRUE, cluster_cols = FALSE,   
                                annotation_col = annotation_columns, show_colnames = FALSE, color = my_color_palette(100))  

ggsave("results/figures/pam_diff_log10.png", plot = pam50_heatmap_log10, width = 8, height = 5, units = "in", dpi = 300)  

# ANCOMBC2 分析 (cut_age: 55) (ANCOMBC2 Analysis - Cut Age: 55)  
ancombc2_results_55 <- ancombc2(tse_tumor, group = "PAM50",   
                                fix_formula = "PAM50 + age_group_55",   
                                rand_formula = "(1|batch)",   
                                pairwise = TRUE,  
                                global = TRUE,  
                                struc_zero = TRUE,   
                                prv_cut = 0.10,  
                                p_adj_method = "fdr",  
                                pseudo_sens = TRUE,  
                                n_cl = 8)  

ancombc2_res_55 <- ancombc2_results_55$res  
taxa_55 <- ancombc2_res_55$taxon[ancombc2_res_55$`diff_age_group_55>=55` == "TRUE"]  


# ANCOMBC2 分析 (cut_age: 45-60) (ANCOMBC2 Analysis - Cut Age: 45-60)  
ancombc2_results_45_60 <- ancombc2(tse_tumor, group = "PAM50",   
                                   fix_formula = "PAM50 + age_group_45_60",   
                                   rand_formula = "(1|batch)",   
                                   pairwise = TRUE,  
                                   global = TRUE,  
                                   struc_zero = TRUE,   
                                   prv_cut = 0.10,  
                                   p_adj_method = "fdr",  
                                   pseudo_sens = TRUE,  
                                   n_cl = 8)  

ancombc2_res_45_60 <- ancombc2_results_45_60$res  
taxa_45_60 <- ancombc2_res_45_60$taxon[ancombc2_res_45_60$`diff_age_group_45_6042~60` == "TRUE" | ancombc2_res_45_60$`diff_age_group_45_60>60` == "TRUE"]  
filtered_micro_data_age_45_60 <- micro_data_tumor_tss[, colnames(micro_data_tumor_tss) %in% taxa_45_60, drop = FALSE] %>% t() %>% as.data.frame()  

annotation_columns_age_45_60 <- meta_data_WGS_tumor[order(meta_data_WGS_tumor$age), c("age_group_45_60", "PAM50", "batch", "age_group_55", "age","sample")]   
rownames(annotation_columns_age_45_60) <- annotation_columns_age_45_60$sample  
annotation_columns_age_45_60 <- annotation_columns_age_45_60[, c("age_group_45_60", "age_group_55", "age", "PAM50", "batch")]  

filtered_micro_data_age_45_60 <- filtered_micro_data_age_45_60[, rownames(annotation_columns_age_45_60)]  
rownames(filtered_micro_data_age_45_60) <- gsub(".*\\|s__", "s__", rownames(filtered_micro_data_age_45_60))  

# 年龄组 45-60 的热图 (Heatmap for Age Group 45-60)   
age_45_60_heatmap <- pheatmap(filtered_micro_data_age_45_60, cluster_rows = TRUE, cluster_cols = FALSE,  
                              annotation_col = annotation_columns_age_45_60, show_colnames = FALSE,color = my_color_palette(100))  
ggsave("results/figures/ancombc2_heatmap_age_45_60.png", plot = age_45_60_heatmap, width = 8, height = 5, units = "in", dpi = 300)  

# log10 转换 (Log10 Transformation)  
filtered_micro_data_age_45_60_log10 <- log10(filtered_micro_data_age_45_60 + min(filtered_micro_data_age_45_60[filtered_micro_data_age_45_60 != 0]))  
pam_heatmap_age_45_60_log10 <- pheatmap(filtered_micro_data_age_45_60_log10, cluster_rows = TRUE, cluster_cols = FALSE,   
                                        annotation_col = annotation_columns_age_45_60, show_colnames = FALSE, color = my_color_palette(100))  
ggsave("results/figures/ancombc2_heatmap_age_45_60_log10.png", plot = pam_heatmap_age_45_60_log10, width = 8, height = 5, units = "in", dpi = 300)  

# ANCOMBC2 分析 (年龄) (ANCOMBC2 Analysis - By Age)  
ancombc2_results_age <- ancombc2(tse_tumor, group = "PAM50",   
                                 fix_formula = "PAM50 + age",   
                                 rand_formula = "(1|batch)",   
                                 pairwise = TRUE,  
                                 global = TRUE,  
                                 struc_zero = TRUE,   
                                 prv_cut = 0.10,  
                                 p_adj_method = "fdr",  
                                 pseudo_sens = TRUE,  
                                 n_cl = 16)  
ancombc2_res_age <- ancombc2_results_age$res   

diff_columns_age <- c("diff_(Intercept)", "diff_PAM50Her2", "diff_PAM50LumA", "diff_PAM50LumB", "diff_PAM50Normal", "diff_age")  

# 提取以 diff_ 开头的列中存在 TRUE 的行 (Extract Rows with TRUE in Columns Starting with diff_)  
filtered_ancombc2_res_age <- ancombc2_res_age %>%  
  filter(rowSums(.[, diff_columns_age] == TRUE) > 0)  

# 提取差异物种 (Extract Differential Taxa)  
taxa_age <- filtered_ancombc2_res_age$taxon  
filtered_micro_data_age <- micro_data_tumor_tss[, colnames(micro_data_tumor_tss) %in% taxa_age, drop = FALSE] %>% t() %>% as.data.frame()  

annotation_columns_age <- meta_data_WGS_tumor[order(meta_data_WGS_tumor$age), c("age", "PAM50", "batch", "age_group_55", "age_group_45_60","sample","LNM.0.No.1.Yes.")]  
rownames(annotation_columns_age) <- annotation_columns_age$sample  
annotation_columns_age <- annotation_columns_age[, c("age", "age_group_55", "age_group_45_60","PAM50", "batch","LNM.0.No.1.Yes.")]  

filtered_micro_data_age <- filtered_micro_data_age[, rownames(annotation_columns_age)]  

rownames(filtered_micro_data_age) <- gsub(".*\\|s__", "s__", rownames(filtered_micro_data_age))  

# 年龄相关热图 (Age-Related Heatmap)  
age_heatmap_final <- pheatmap(filtered_micro_data_age, cluster_rows = TRUE, cluster_cols = FALSE,   
                              annotation_col = annotation_columns_age, show_colnames = FALSE, color = my_color_palette(100))  
ggsave("results/figures/ancombc2_heatmap_age_final.png", plot = age_heatmap_final, width = 8, height = 5, units = "in", dpi = 300)

annotation_columns_age <- annotation_columns_age[, c("PAM50","age", "age_group_55", "age_group_45_60", "batch","LNM.0.No.1.Yes.")]
annotation_columns_age <- annotation_columns_age[order(annotation_columns_age$PAM50),]
filtered_micro_data_age <- filtered_micro_data_age[, rownames(annotation_columns_age)]

age_heatmap_PAM50_final <- pheatmap(filtered_micro_data_age, cluster_rows = TRUE, cluster_cols = FALSE,   
                              annotation_col = annotation_columns_age, show_colnames = FALSE, color = my_color_palette(100))
ggsave("results/figures/ancombc2_heatmap_age_PAM50.png", plot = age_heatmap_PAM50_final, width = 8, height = 5, units = "in", dpi = 300)
