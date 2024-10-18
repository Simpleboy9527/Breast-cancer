# 设置随机种子  
set.seed(123)  

# 正规化 TSS 处理并转换为数据框  
micro_data_clean <- norm_tss(tse_all@otu_table) %>% as.data.frame()  # seed 123 , depth 947  

# Alpha 多样性计算  
alpha_diversity <- micro_data_clean %>%  
  mutate_all(as.numeric) %>%  
  as.matrix() %>%  
  microbiome::alpha(index = "all")  
write.csv(alpha_diversity, "results/alpha_diversity.csv", row.names = TRUE)  

# Alpha 多样性绘图  
alpha_plot <- list()  
for (i in "Group") {  
  alpha_plot[[1]] <- alpha_plot_fun(alpha_diversity, meta_data_metagenome, group_var = i)  
}  
combine_alpha_plot <- alpha_plot[[1]]  
ggsave("results/figures/alpha_plot_all.png", plot = combine_alpha_plot, width = 6, height = 4.5, units = "in", dpi = 300)  

# Beta 多样性计算  
micro_data_clean <- micro_data_clean[, meta_data_metagenome$sample]  
beta_diversity <- beta_diversity_fun(micro_data_clean, meta_data_metagenome)  

# Beta 多样性绘图  
beta_plot_bray <- beta_plot_fun(beta_diversity[[1]], meta_data_metagenome, variables_to_plot = c("Group", "PAM50", "age_group_45_60", "age", "stage", "LNM.0.No.1.Yes."))  
beta_plot_jaccard <- beta_plot_fun(beta_diversity[[2]], meta_data_metagenome, variables_to_plot = c("Group", "PAM50", "age_group_45_60", "age", "stage", "LNM.0.No.1.Yes."))  
combine_beta_plot <- beta_plot_bray / beta_plot_jaccard  
ggsave("results/figures/combine_beta_plot.png", plot = combine_beta_plot, width = 25, height = 6, units = "in", dpi = 300)  

# 肿瘤组织 Beta 多样性分析  
beta_diversity_tumor <- beta_diversity_fun(micro_data_tumor, meta_data_WGS_tumor)  
beta_plot_bray_tumor <- beta_plot_fun(beta_diversity_tumor[[1]], meta_data_WGS_tumor, variables_to_plot = c("age", "PAM50", "stage", "LNM.0.No.1.Yes.", "Grade", "HER2..clean."))  
beta_plot_jaccard_tumor <- beta_plot_fun(beta_diversity_tumor[[2]], meta_data_WGS_tumor, variables_to_plot = c("age", "PAM50", "stage", "LNM.0.No.1.Yes.", "Grade", "HER2..clean."))  
combine_beta_plot_tumor <- beta_plot_bray_tumor / beta_plot_jaccard_tumor  
ggsave("results/figures/combine_beta_plot_tumor.png", plot = combine_beta_plot_tumor, width = 20, height = 6, units = "in", dpi = 300)  

# Adonis 检验  
select_vars <- c("Group", "stage", "PAM50", "age", "ER.clean.1..", "Grade", "PR..clean..1..", "HER2..clean.", "Ki67..clean.", "type_of_organization", "batch", "LNM.0.No.1.Yes.")  
adonis_summary <- adonis_fun(beta_diversity[[1]], meta_data_metagenome, selected_variables = select_vars)  

# 绘制 Adonis 检验图  
adonis_summary <- adonis_summary %>%  
  arrange(R2) %>%  
  mutate(Variable = factor(Variable, levels = Variable))  

adonis_summary_plot <- ggplot(adonis_summary, aes(x = Variable, y = R2)) +  
  geom_bar(stat = "identity", fill = "skyblue") +  
  theme_minimal() +  
  geom_text(aes(label = ifelse(PrF < 0.001, "***", ifelse(PrF < 0.01, "**", ifelse(PrF < 0.05, "*", "")))),   
            hjust = -0.5, size = 4) +  
  coord_flip() +  
  xlab(NULL) +  
  ylab("R²") +  
  ggtitle("Adonis2 R² for each variable") +  
  theme(  
    plot.title = element_text(hjust = 0.5),  
    axis.text.x = element_text(size = 12, family = "Arial"),  
    axis.text.y = element_text(size = 12, family = "Arial", hjust = 1),  
    axis.title = element_text(size = 8, family = "Arial"),  
    legend.text = element_text(size = 8, family = "Arial"),  
    legend.title = element_text(size = 8, family = "Arial")  
  )  

ggsave("results/figures/adonis_summary_plot.png", plot = adonis_summary_plot, width = 6, height = 4.5, units = "in", dpi = 300)  

# 肿瘤组织 Adonis 检验  
select_vars_tumor <- c("stage", "PAM50", "age", "ER.clean.1..", "Grade", "PR..clean..1..", "HER2..clean.", "Ki67..clean.", "type_of_organization", "batch","OS")  
adonis_summary_tumor <- adonis_fun(beta_diversity_tumor[[1]], meta_data_WGS_tumor, selected_variables = select_vars_tumor)  

adonis_summary_tumor <- adonis_summary_tumor %>%  
  arrange(R2) %>%  
  mutate(Variable = factor(Variable, levels = Variable))  

adonis_summary_tumor_plot <- ggplot(adonis_summary_tumor, aes(x = Variable, y = R2)) +  
  geom_bar(stat = "identity", fill = "skyblue") +  
  theme_minimal() +  
  geom_text(aes(label = ifelse(PrF < 0.001, "***", ifelse(PrF < 0.01, "**", ifelse(PrF < 0.05, "*", "")))),   
            hjust = -0.5, size = 4) +  
  coord_flip() +  
  xlab(NULL) +  
  ylab("R²") +  
  ggtitle("Adonis2 R² for each variable in tumor tissue") +  
  theme(  
    plot.title = element_text(hjust = 0.5),  
    axis.text.x = element_text(size = 12, family = "Arial"),  
    axis.text.y = element_text(size = 12, family = "Arial", hjust = 1),  
    axis.title = element_text(size = 8, family = "Arial"),  
    legend.text = element_text(size = 8, family = "Arial"),  
    legend.title = element_text(size = 8, family = "Arial")  
  )  

ggsave("results/figures/adonis_summary_tumor_plot.png", plot = adonis_summary_tumor_plot, width = 7, height = 4.5, units = "in", dpi = 300)  

# 血液样本的 Adonis 检验  
meta_data_blood <- meta_data_metagenome %>% filter(Group == "Blood")  
adonis_summary_blood <- adonis_fun(beta_diversity[[1]], meta_data_blood, selected_variables = select_vars)  

adonis_summary_blood <- adonis_summary_blood %>%  
  arrange(R2) %>%  
  mutate(Variable = factor(Variable, levels = Variable))  

adonis_summary_blood_plot <- ggplot(adonis_summary_blood, aes(x = Variable, y = R2)) +  
  geom_bar(stat = "identity", fill = "skyblue") +  
  theme_minimal() +  
  geom_text(aes(label = ifelse(PrF < 0.001, "***", ifelse(PrF < 0.01, "**", ifelse(PrF < 0.05, "*", "")))),   
            hjust = -0.5, size = 4) +  
  coord_flip() +  
  xlab(NULL) +  
  ylab("R²") +  
  ggtitle("Adonis2 R² for each variable in blood samples") +  
  theme(  
    plot.title = element_text(hjust = 0.5),  
    axis.text.x = element_text(size = 12, family = "Arial"),  
    axis.text.y = element_text(size = 12, family = "Arial", hjust = 1),  
    axis.title = element_text(size = 8, family = "Arial"),  
    legend.text = element_text(size = 8, family = "Arial"),  
    legend.title = element_text(size = 8, family = "Arial")  
  )  

ggsave("results/figures/adonis_summary_blood_plot.png", plot = adonis_summary_blood_plot, width = 7, height = 4.5, units = "in", dpi = 300)  

# 合并 Alpha 多样性、Meta 数据和 Micro 数据  
merge_data <- merge(meta_data_metagenome, alpha_diversity, by.x = "sample", by.y = "row.names")  
micro_data_t <- t(micro_data) %>% as.data.frame()  
merge_data <- merge(merge_data, micro_data_t, by.x = "sample", by.y = "row.names")  
write.csv(merge_data, "results/merge_data.csv", row.names = TRUE)  

# 异常样本分析  
# 选择样本  
bray_distance <- beta_diversity[[2]] %>% as.matrix()  
temp_data <- cmdscale(bray_distance, k = 3, eig = TRUE)  
dune_pcoa_points <- as.data.frame(temp_data$points)  

temp_plot <- ggplot(dune_pcoa_points, aes(x = V1, y = V2)) +  
  geom_point(aes(color = meta_data_metagenome$Group)) +  
  theme_minimal() +  
  xlab("PCoA1") +  
  ylab("PCoA2") +  
  ggtitle("PCoA plot") +  
  theme(  
    plot.title = element_text(hjust = 0.5),  
    axis.text.x = element_text(size = 12, family = "Arial"),  
    axis.text.y = element_text(size = 12, family = "Arial"),  
    axis.title = element_text(size = 8, family = "Arial"),  
    legend.text = element_text(size = 8, family = "Arial"),  
    legend.title = element_text(size = 8, family = "Arial")  
  )  
print(temp_plot)  

# 找到异常样本  
abnormal_samples <- rownames(dune_pcoa_points)[dune_pcoa_points$V1 < 0 & dune_pcoa_points$V2 < -0.2]  
meta_data_abnormal <- meta_data_metagenome %>% filter(sample %in% abnormal_samples)  
write.csv(meta_data_abnormal, "results/meta_data_abnormal.csv")  
micro_data_abnormal <- micro_data_clean[, meta_data_abnormal$sample]  

# PAM50 统计  
PAM50_abnormal_plot <- ggplot(meta_data_abnormal, aes(x = PAM50)) +  
  geom_bar(fill = "skyblue") +  
  theme_minimal() +  
  xlab("PAM50") +  
  ylab("Sample number") +  
  ggtitle("PAM50 distribution in abnormal samples") +  
  theme(  
    plot.title = element_text(hjust = 0.5),  
    axis.text.x = element_text(size = 12, family = "Arial"),  
    axis.text.y = element_text(size = 12, family = "Arial"),  
    axis.title = element_text(size = 8, family = "Arial"),  
    legend.text = element_text(size = 8, family = "Arial"),  
    legend.title = element_text(size = 8, family = "Arial")  
  )  

# 在异常样本中进行 Adonis 检验  
beta_diversity_abnormal <- beta_diversity_fun(micro_data_abnormal, meta_data_abnormal)  
adonis_summary_abnormal <- adonis_fun(beta_diversity_abnormal[[1]], meta_data_abnormal, selected_variables = select_vars)  

adonis_summary_abnormal <- adonis_summary_abnormal %>%  
  arrange(R2) %>%  
  mutate(Variable = factor(Variable, levels = Variable))  

adonis_summary_abnormal_plot <- ggplot(adonis_summary_abnormal, aes(x = Variable, y = R2)) +  
  geom_bar(stat = "identity", fill = "skyblue") +  
  theme_minimal() +  
  geom_text(aes(label = ifelse(PrF < 0.001, "***", ifelse(PrF < 0.01, "**", ifelse(PrF < 0.05, "*", "")))),   
            hjust = -0.5, size = 4) +  
  coord_flip() +  
  xlab(NULL) +  
  ylab("R²") +  
  ggtitle("Adonis2 R² for each variable in abnormal samples") +  
  theme(  
    plot.title = element_text(hjust = 0.5),  
    axis.text.x = element_text(size = 12, family = "Arial"),  
    axis.text.y = element_text(size = 12, family = "Arial", hjust = 1),  
    axis.title = element_text(size = 8, family = "Arial"),  
    legend.text = element_text(size = 8, family = "Arial"),  
    legend.title = element_text(size = 8, family = "Arial")  
  )  

ggsave("results/figures/adonis_summary_abnormal_plot.png", plot = adonis_summary_abnormal_plot, width = 7, height = 4.5, units = "in", dpi = 300)  

# 处理正常样本数据  
meta_data_WGS_temp <- meta_data_metagenome %>%  
  mutate(beta_Group = case_when(  
    sample %in% abnormal_samples ~ "Abnormal",  
    TRUE ~ "Normal"  
  )) %>%  
  filter(Group %in% unique(meta_data_abnormal$Group))  

micro_data_clean_temp <- micro_data_clean[, meta_data_WGS_temp$sample]  
tse_temp <- phyloseq_creat(micro_data_clean, meta_data_WGS_temp)

# 物种组成绘图  
# taxa_stack_plot 函数用于绘制物种组成图  
for (i in c("Genus", "Species")) {  
  tse_temp_plot <- taxa_stack_plot(tse_temp, group_var = "beta_Group", taxa_level = i, top_n = 15)  
  ggsave(paste0("results/figures/tse_temp_plot_", i, ".png"), plot = tse_temp_plot, width = 12, height = 6, units = "in", dpi = 300)  
}  

# Alpha 多样性分析  
alpha_plot_abnormal <- alpha_plot_fun(alpha_diversity, meta_data_WGS_temp, group_var = "beta_Group")  
ggsave("results/figures/alpha_plot_abnormal.png", plot = alpha_plot_abnormal, width = 5, height = 4.5, units = "in", dpi = 300)  

# 物种组成分析 - 所有组  
for (i in c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {  
  tse_all_plot <- taxa_stack_plot(tse_all, group_var = "Group", taxa_level = i, top_n = 15, faced_group_ordered = c("Tumor_tissue", "Adjacent_tissue", "Blood"))  
  ggsave(paste0("results/figures/tse_all_", i, ".png"), plot = tse_all_plot, width = 18, height = 6, units = "in", dpi = 300)  
}  

# 肿瘤组织物种组成分析 - 按 PAM50 分组  
for (i in c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {  
  tse_tumor_plot <- taxa_stack_plot(tse_tumor, group_var = "PAM50", taxa_level = i, top_n = 15, faced_group_ordered = c("LumA", "LumB", "Her2", "Basal", "Normal"))  
  ggsave(paste0("results/figures/tse_tumor_", i, ".png"), plot = tse_tumor_plot, width = 15, height = 6, units = "in", dpi = 300)  
}  

# 按年龄分组的肿瘤组织物种组成分析  
for (i in c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {  
  tse_tumor_plot <- taxa_stack_plot(tse_tumor, group_var = "age_group_55", taxa_level = i, top_n = 15, faced_group_ordered = c("<55", ">=55"))  
  ggsave(paste0("results/figures/tse_tumor_age_group_55_", i, ".png"), plot = tse_tumor_plot, width = 12, height = 6, units = "in", dpi = 300)  
}  

# 按年龄组 45-60 的肿瘤组织物种组成分析  
for (i in c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {  
  tse_tumor_plot <- taxa_stack_plot(tse_tumor, group_var = "age_group_45_60", taxa_level = i, top_n = 15, faced_group_ordered = c("0~41", "42~60", ">60"))  
  ggsave(paste0("results/figures/tse_tumor_age_group_45_60_", i, ".png"), plot = tse_tumor_plot, width = 12, height = 6, units = "in", dpi = 300)  
}  

# 肿瘤组织的 Alpha 多样性分析  
alpha_plot_tumor_PAM50 <- alpha_plot_fun(alpha_diversity, meta_data_WGS_tumor, group_var = "PAM50")  
ggsave("results/figures/alpha_plot_tumor_PAM50.png", plot = alpha_plot_tumor_PAM50, width = 5, height = 4.5, units = "in", dpi = 300)  

alpha_plot_tumor_age_group_55 <- alpha_plot_fun(alpha_diversity, meta_data_WGS_tumor, group_var = "age_group_55")  
ggsave("results/figures/alpha_plot_tumor_age_group_55.png", plot = alpha_plot_tumor_age_group_55, width = 5, height = 4.5, units = "in", dpi = 300)  

alpha_plot_tumor_age_group_45_60 <- alpha_plot_fun(alpha_diversity, meta_data_WGS_tumor, group_var = "age_group_45_60")  
ggsave("results/figures/alpha_plot_tumor_age_group_45_60.png", plot = alpha_plot_tumor_age_group_45_60, width = 5, height = 4.5, units = "in", dpi = 300)  

# 年龄与 Alpha 多样性的点图  
temp_data <- meta_data_WGS_tumor %>%  
  mutate(observed = alpha_diversity$observed[match(sample, rownames(alpha_diversity))],  
         evenness_pielou = alpha_diversity$evenness_pielou[match(sample, rownames(alpha_diversity))],  
         shannon = alpha_diversity$diversity_shannon[match(sample, rownames(alpha_diversity))])  

temp_plot <- list()  
for (i in c("observed", "shannon", "evenness_pielou")) {  
  cor_test <- cor.test(temp_data$age, temp_data[[i]], method = "spearman")  
  cor_value <- round(cor_test$estimate, 2)  # 相关系数  
  p_value <- round(cor_test$p.value, 4)     # p值  
  
  temp_plot[[i]] <- ggplot(temp_data, aes(x = age, y = !!sym(i))) +  
    geom_point(aes(color = PAM50)) +  
    geom_smooth(method = "lm", se = TRUE) +  
    theme_minimal() +  
    xlab("Age") +  
    ylab(i) +  
    theme(  
      plot.title = element_text(hjust = 0.5),  
      axis.text.x = element_text(size = 12, family = "Arial"),  
      axis.text.y = element_text(size = 12, family = "Arial"),  
      axis.title = element_text(size = 8, family = "Arial"),  
      legend.text = element_text(size = 8, family = "Arial"),  
      legend.title = element_text(size = 8, family = "Arial")  
    ) +  
    annotate("text", x = Inf, y = Inf,  
             label = paste("Spearman r =", cor_value, "\np =", p_value),  
             hjust = 1.1, vjust = 1.5, size = 4, color = "black",  
             fontface = "italic")  # 添加Spearman相关系数和p值  
}  

# 合并相关性图  
cor_plot <- temp_plot[[1]] + temp_plot[[2]] + temp_plot[[3]] + plot_layout(guides = "collect")  
ggsave("results/figures/cor_plot.png", plot = cor_plot, width = 12, height = 4, units = "in", dpi = 300)  

# 保存结果  
save(combine_alpha_plot, combine_beta_plot, combine_beta_plot_tumor, 
     adonis_summary_plot, adonis_summary_tumor_plot,adonis_summary_blood_plot,adonis_summary_abnormal_plot,
     alpha_plot_tumor_PAM50, alpha_plot_tumor_age_group_55, alpha_plot_tumor_age_group_45_60,
     cor_plot,file = "results/alpha_beta_diversity.RData")
remove_non_functions() 