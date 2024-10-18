#meta_data statistics
#临床特征统计
#pam50
pam50_stats_plot <- ggplot(meta_data, aes(x = PAM50)) +  
  geom_bar() +  
  ylab("Sample number") +  
  xlab(NULL) +  
  theme_minimal() +  
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5)  

pam50_stats_plot
ggsave("results/figures/data_statistic/pam50_stats_plot.png", plot = pam50_stats_plot, width = 4, height = 4, units = "in", dpi = 300)

#TNM
TNM_stats_plot <- ggplot(meta_data,aes(x= TNM))+
  geom_bar()+
  ylab("Sample number")+
  xlab(NULL)+
  theme_minimal()+
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5) 
TNM_stats_plot
ggsave("results/figures/data_statistic/TNM_stats_plot.png", plot = TNM_stats_plot, width = 12, height = 4, units = "in", dpi = 300)

#Grade
Grade_stats_plot <- ggplot(meta_data,aes(x= Grade))+
  geom_bar()+
  ylab("Sample number")+
  xlab(NULL)+
  theme_minimal()+
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5) 
Grade_stats_plot
ggsave("results/figures/data_statistic/Grade_stats_plot.png", plot = Grade_stats_plot, width = 4, height = 4, units = "in", dpi = 300)

#age
# 假设你想将年龄分成以下区间  
breaks <- c(0, 20, 30, 40, 50, 60, 70, 80, 90, 100)  
labels <- c("0-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100")  

# 创建一个新的变量 age_group，处理 NA 值  
meta_data$age_group <- cut(meta_data$age, breaks = breaks, labels = labels, right = FALSE)  

# 绘制直方图  
age_histogram_plot <- ggplot(meta_data, aes(x = age_group)) +  
  geom_bar() +  
  ylab("Sample number") +  
  xlab("Age Group") +  
  theme_minimal() +  
  geom_text(stat = 'count', aes(label = ..count..),   
            vjust = -0.5)  # 添加计数标签  
age_histogram_plot
ggsave("results/figures/data_statistic/age_histogram_plot.png", plot = age_histogram_plot, width = 5, height = 4, units = "in", dpi = 300)


#stage
stage_stats_plot <- ggplot(meta_data,aes(x= stage))+
  geom_bar()+
  ylab("Sample number")+
  xlab(NULL)+
  theme_minimal()+
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5) 
stage_stats_plot
ggsave("results/figures/data_statistic/stage_stats_plot.png", plot = stage_stats_plot, width = 4, height = 4, units = "in", dpi = 300)


#map_reads statistics
# Identify indices for Bacteria, Archaea, and Viruses
indices <- which(map_reads$`#Classification` %in% c("d__Bacteria", "d__Archaea", "d__Viruses"))

# Filter, transpose, and set column names
map_reads_filter <- map_reads[indices, ] %>% t()
colnames(map_reads_filter) <- map_reads_filter[1, ]
map_reads_filter <- map_reads_filter[-1, ] %>% as.data.frame()
map_reads_filter$sample <- rownames(map_reads_filter)

merged_data <- merge(map_reads_filter, meta_data_metagenome, by = "sample")
save(merged_data, file = "results/merged_data.RData")
remove(merged_data)
map_reads_filter <- map_reads_filter %>% dplyr::select(-sample)

# Add Total_reads and unclassfied_reads columns
map_reads_filter <- map_reads_filter %>%
  mutate(Total_reads = total_reads$Clean_Reads[match(rownames(map_reads_filter), total_reads$Sample)],
         unclassfied_reads = unclassfied_reads$reads[match(rownames(map_reads_filter), unclassfied_reads$Sample)])

# Convert all columns to numeric and calculate human_reads
map_reads_filter <- map_reads_filter %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(human_reads = Total_reads - d__Bacteria - d__Archaea - d__Viruses - unclassfied_reads)


# Normalize the data by dividing each value by the total reads
map_reads_filter_prop <- map_reads_filter %>%
  mutate(across(everything(), ~ ./ Total_reads))

# Add a 'group' column based on the sample names
map_reads_filter_prop <- map_reads_filter_prop %>%
  mutate(
    group = ifelse(grepl("T$", rownames(map_reads_filter_prop)), "Tumor_tissue",
                   ifelse(grepl("N$", rownames(map_reads_filter_prop)), "Adjacent_tissue", "Blood"))
  )


# Convert 'group' column to a factor with specified levels
map_reads_filter_prop <- map_reads_filter_prop[rownames(map_reads_filter_prop) %in% meta_data_metagenome$sample,]
map_reads_filter_prop$PAM50 <- meta_data_metagenome$PAM50[match(rownames(map_reads_filter_prop), meta_data_metagenome$sample)]
map_reads_filter_prop$stage <- meta_data_metagenome$stage[match(rownames(map_reads_filter_prop), meta_data_metagenome$sample)]
map_reads_filter_prop$group <- factor(map_reads_filter_prop$group, levels = c("Tumor_tissue", "Adjacent_tissue", "Blood"))


# Create a boxplot with beeswarm overlay for Bacteria abundance
fig0_0 <- ggplot(map_reads_filter_prop, aes(x = group, y = log10(d__Bacteria))) +
  geom_boxplot() +
  geom_beeswarm(alpha = 0.5) +
  xlab(NULL) +
  theme_minimal() +
  ylab("log10(Bacteria_Abundance)") +
  ggtitle("not removing contamination Bacteria") +
  theme(plot.title = element_text(hjust = 0.5))

# Save the plot as a PNG file
ggsave("results/figures/fig0_0.png", plot = fig0_0, width = 4, height = 3, units = "in", dpi = 300)

fig0_1 <- ggplot(map_reads_filter_prop, aes(x = PAM50, y = log10(d__Bacteria))) +
  geom_boxplot() +
  geom_beeswarm(alpha = 0.5) +
  xlab(NULL) +
  theme_minimal() +
  ylab("log10(Bacteria_Abundance)") +
  ggtitle("not removing contamination Bacteria") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("results/figures/fig0_1.png", plot = fig0_1, width = 4, height = 3, units = "in", dpi = 300)

fig0_2 <- ggplot(map_reads_filter_prop, aes(x = stage, y = log10(d__Bacteria))) +
  geom_boxplot() +
  geom_beeswarm(alpha = 0.5) +
  xlab(NULL) +
  theme_minimal() +
  ylab("log10(Bacteria_Abundance)") +
  ggtitle("not removing contamination Bacteria") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("results/figures/fig0_2.png", plot = fig0_2, width = 4, height = 3, units = "in", dpi = 300)

save(age_histogram_plot, pam50_stats_plot, TNM_stats_plot, Grade_stats_plot, stage_stats_plot, fig0_0, fig0_1, fig0_2, file = "results/data_statistic_plot.RData")

remove(age_histogram_plot, pam50_stats_plot, TNM_stats_plot, Grade_stats_plot, stage_stats_plot, fig0_0, fig0_1, fig0_2, map_reads_filter_prop, map_reads_filter, map_reads, total_reads, unclassfied_reads)


#meta_data statistics

#热图
temp_data <- meta_data %>% mutate(
  WGS_T_tissue = ifelse(is.na(WGS_T_tissue), 0, 1),
  WGS_N_tissue = ifelse(is.na(WGS_N_tissue), 0, 1),
  WGS_N_blood = ifelse(is.na(WGS_N_blood), 0, 1),
  RNA_NR_Tissue = ifelse(is.na(RNA_NR_Tissue), 0, 1),
  RNA_TR_Tissue = ifelse(is.na(RNA_TR_Tissue), 0, 1),
  WGBS_T_tissue = ifelse(is.na(WGBS_T_tissue), 0, 1),
  WGBS_N_tissue = ifelse(is.na(WGBS_N_tissue), 0, 1),
  WGBS_N_blood = ifelse(is.na(WGBS_N_blood), 0, 1)
) %>% dplyr::select(WGBS_N_blood, WGBS_N_tissue, WGBS_T_tissue, RNA_TR_Tissue, 
                    RNA_NR_Tissue, WGS_N_blood, WGS_N_tissue, WGS_T_tissue,
                    age, PAM50, stage, Grade, TNM,
                    Ki67,LNM.0.No.1.Yes.,PR..clean..1..,ER.clean.1..,HER2..clean.,patient_record_number)

#热图
annotation_col <- temp_data %>% dplyr::select(age,PAM50, stage, Grade, TNM, Ki67, LNM.0.No.1.Yes., PR..clean..1.., ER.clean.1.., HER2..clean.,patient_record_number)
rownames(annotation_col) <- temp_data$patient_record_number
annotation_col <- annotation_col %>% dplyr::select(-patient_record_number)

# Create a heatmap of the metadata
rownames(temp_data) <- temp_data$patient_record_number
temp_data <- temp_data %>% dplyr::select(WGBS_N_blood, WGBS_N_tissue, WGBS_T_tissue, RNA_TR_Tissue, 
                                         RNA_NR_Tissue, WGS_N_blood, WGS_N_tissue, WGS_T_tissue)
temp_data <- temp_data %>% t() %>% as.data.frame()

long_data <- temp_data %>%  
  rownames_to_column(var = "Sample") %>%  
  pivot_longer(-Sample, names_to = "Feature", values_to = "Value") %>%  
  mutate(Value = factor(Value, levels = c(0, 1)))  # 将Value转换为因子



# 绘制热图  
long_data <- long_data %>%
  mutate(PAM50 = annotation_col$PAM50[match(Feature, rownames(annotation_col))],
         Stage = annotation_col$stage[match(Feature, rownames(annotation_col))],
         Grade = annotation_col$Grade[match(Feature, rownames(annotation_col))],
         TNM = annotation_col$TNM[match(Feature, rownames(annotation_col))],
         LNM.0.No.1.Yes. = annotation_col$LNM.0.No.1.Yes.[match(Feature, rownames(annotation_col))],
         PR_stage = annotation_col$PR..clean..1..[match(Feature, rownames(annotation_col))],
         ER_stage = annotation_col$ER.clean.1..[match(Feature, rownames(annotation_col))],
         HER2_stage = annotation_col$HER2..clean.[match(Feature, rownames(annotation_col))],
         Age = annotation_col$age[match(Feature, rownames(annotation_col))])

long_data <- long_data %>%  
  mutate(across(where(is.character), ~ na_if(., "")))  # 仅将字符列中的空字符串转换为NA 

long_data <- long_data[order(long_data$Age,decreasing = F),]
long_data$Feature <- factor(long_data$Feature,levels = rev(unique(long_data$Feature)),ordered = T)
long_data$Grade <- factor(long_data$Grade)
long_data$LNM.0.No.1.Yes. <- factor(long_data$LNM.0.No.1.Yes.)

temp_pheatmap <- ggplot(long_data, aes(x = Feature, y = Sample)) +  
  geom_tile(aes(fill = Value)) +  # 使用因子Value  
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue")) +  # 手动设置颜色  
  theme_minimal() +  
  labs(x = "Sample", y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",         
        panel.grid.major = element_blank(),  # 去除主要网格线           
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
temp_pheatmap

temp_age <- ggplot(long_data, aes(x = Feature, y = "Age")) +  
  geom_tile(aes(fill = Age)) +  # 使用因子Value  
  theme_minimal() +  
  labs(x = NULL, y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",         
        panel.grid.major = element_blank(),  # 去除主要网格线           
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
temp_age

temp_PAM50 <- ggplot(long_data, aes(x = Feature, y = "PAM50")) +  
  geom_tile(aes(fill = PAM50)) +  # 使用因子Value  
  theme_minimal() +  
  scale_fill_manual(values = my_pal) +  # 手动设置颜色
  labs(x = NULL, y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",         
        panel.grid.major = element_blank(),  # 去除主要网格线           
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
temp_PAM50

temp_Stage <- ggplot(long_data, aes(x = Feature, y = "Stage")) +  
  geom_tile(aes(fill = Stage)) +  # 使用因子Value  
  theme_minimal() +  
  scale_fill_manual(values = my_pal) +  # 手动设置颜色
  labs(x = NULL, y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",         
        panel.grid.major = element_blank(),  # 去除主要网格线           
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
temp_Stage

temp_Grade <- ggplot(long_data, aes(x = Feature, y = "Grade")) +  
  geom_tile(aes(fill = Grade)) +  # 使用因子Value  
  theme_minimal() +  
  scale_fill_manual(values = my_pal) +  # 手动设置颜色
  labs(x = NULL, y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",         
        panel.grid.major = element_blank(),  # 去除主要网格线           
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
temp_Grade

temp_LNM <- ggplot(long_data, aes(x = Feature, y = "LNM.0.No.1.Yes.")) +  
  geom_tile(aes(fill = LNM.0.No.1.Yes.)) +  # 使用因子Value  
  theme_minimal() +  
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue")) +  # 手动设置颜色
  labs(x = NULL, y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",         
        panel.grid.major = element_blank(),  # 去除主要网格线           
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
temp_LNM

temp_PR <- ggplot(long_data, aes(x = Feature, y = "PR_stage")) +  
  geom_tile(aes(fill = PR_stage)) +  # 使用因子Value  
  theme_minimal() +  
  scale_fill_manual(values = my_pal) +  # 手动设置颜色
  labs(x = NULL, y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",
        panel.grid.major = element_blank(),  # 去除主要网格线           
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
temp_PR

temp_ER <- ggplot(long_data, aes(x = Feature, y = "ER_stage")) +  
  geom_tile(aes(fill = ER_stage)) +  # 使用因子Value  
  theme_minimal() +  
  scale_fill_manual(values = my_pal) +  # 手动设置颜色
  labs(x = NULL, y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",
        panel.grid.major = element_blank(),  # 去除主要网格线  
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
temp_ER

temp_HER2 <- ggplot(long_data, aes(x = Feature, y = "HER2_stage")) +  
  geom_tile(aes(fill = HER2_stage)) +  # 使用因子Value  
  theme_minimal() +  
  scale_fill_manual(values = my_pal) +  # 手动设置颜色
  labs(x = NULL, y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",
        panel.grid.major = element_blank(),  # 去除主要网格线  
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))  # 去除次要网格线
temp_HER2

meta_data_heatmap <- temp_pheatmap %>%  
  aplot::insert_top(temp_age, height = 0.1) %>%  
  aplot::insert_top(temp_PAM50, height = 0.1) %>%  
  aplot::insert_top(temp_Stage, height = 0.1) %>%  
  aplot::insert_top(temp_Grade, height = 0.1) %>%  
  aplot::insert_top(temp_LNM, height = 0.1) %>%  
  aplot::insert_top(temp_PR, height = 0.1) %>%  
  aplot::insert_top(temp_ER, height = 0.1) %>%  
  aplot::insert_top(temp_HER2, height = 0.1) 

# 显示最终的热图  
print(meta_data_heatmap)  
ggsave("results/figures/data_statistic/meta_data_heatmap.png", plot = meta_data_heatmap, width = 12, height = 10, units = "in", dpi = 300)


#meta_data_raw
meta_data_raw <- read.csv("data/metadata.csv")
temp_data <- meta_data_raw %>% mutate(
  WGS_T_tissue = ifelse(is.na(WGS_T_tissue), 0, 1),
  WGS_N_tissue = ifelse(is.na(WGS_N_tissue), 0, 1),
  WGS_N_blood = ifelse(is.na(WGS_N_blood), 0, 1),
  RNA_NR_Tissue = ifelse(is.na(RNA_NR_Tissue), 0, 1),
  RNA_TR_Tissue = ifelse(is.na(RNA_TR_Tissue), 0, 1),
  WGBS_T_tissue = ifelse(is.na(WGBS_T_tissue), 0, 1),
  WGBS_N_tissue = ifelse(is.na(WGBS_N_tissue), 0, 1),
  WGBS_N_blood = ifelse(is.na(WGBS_N_blood), 0, 1)
) %>% dplyr::select(WGBS_N_blood, WGBS_N_tissue, WGBS_T_tissue, RNA_TR_Tissue, 
                    RNA_NR_Tissue, WGS_N_blood, WGS_N_tissue, WGS_T_tissue,
                    age, PAM50, stage, Grade, TNM,
                    Ki67,LNM.0.No.1.Yes.,PR..clean..1..,ER.clean.1..,HER2..clean.,patient_record_number)

#热图
annotation_col <- temp_data %>% dplyr::select(age,PAM50, stage, Grade, TNM, Ki67, LNM.0.No.1.Yes., PR..clean..1.., ER.clean.1.., HER2..clean.,patient_record_number)
rownames(annotation_col) <- temp_data$patient_record_number
annotation_col <- annotation_col %>% dplyr::select(-patient_record_number)

# Create a heatmap of the metadata
rownames(temp_data) <- temp_data$patient_record_number
temp_data <- temp_data %>% dplyr::select(WGBS_N_blood, WGBS_N_tissue, WGBS_T_tissue, RNA_TR_Tissue, 
                                         RNA_NR_Tissue, WGS_N_blood, WGS_N_tissue, WGS_T_tissue)
temp_data <- temp_data %>% t() %>% as.data.frame()

long_data <- temp_data %>%  
  rownames_to_column(var = "Sample") %>%  
  pivot_longer(-Sample, names_to = "Feature", values_to = "Value") %>%  
  mutate(Value = factor(Value, levels = c(0, 1)))  # 将Value转换为因子



# 绘制热图  
long_data <- long_data %>%
  mutate(PAM50 = annotation_col$PAM50[match(Feature, rownames(annotation_col))],
         Stage = annotation_col$stage[match(Feature, rownames(annotation_col))],
         Grade = annotation_col$Grade[match(Feature, rownames(annotation_col))],
         TNM = annotation_col$TNM[match(Feature, rownames(annotation_col))],
         LNM.0.No.1.Yes. = annotation_col$LNM.0.No.1.Yes.[match(Feature, rownames(annotation_col))],
         PR_stage = annotation_col$PR..clean..1..[match(Feature, rownames(annotation_col))],
         ER_stage = annotation_col$ER.clean.1..[match(Feature, rownames(annotation_col))],
         HER2_stage = annotation_col$HER2..clean.[match(Feature, rownames(annotation_col))],
         Age = annotation_col$age[match(Feature, rownames(annotation_col))])

long_data <- long_data %>%  
  mutate(across(where(is.character), ~ na_if(., "")))  # 仅将字符列中的空字符串转换为NA 

long_data <- long_data[order(long_data$Age,decreasing = F),]
long_data$Feature <- factor(long_data$Feature,levels = rev(unique(long_data$Feature)),ordered = T)
long_data$Grade <- factor(long_data$Grade)
long_data$LNM.0.No.1.Yes. <- factor(long_data$LNM.0.No.1.Yes.)

temp_pheatmap <- ggplot(long_data, aes(x = Feature, y = Sample)) +  
  geom_tile(aes(fill = Value)) +  # 使用因子Value  
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue")) +  # 手动设置颜色  
  theme_minimal() +  
  labs(x = "Sample", y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",         
        panel.grid.major = element_blank(),  # 去除主要网格线           
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
temp_pheatmap

temp_age <- ggplot(long_data, aes(x = Feature, y = "Age")) +  
  geom_tile(aes(fill = Age)) +  # 使用因子Value  
  theme_minimal() +  
  labs(x = NULL, y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",         
        panel.grid.major = element_blank(),  # 去除主要网格线           
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
temp_age

temp_PAM50 <- ggplot(long_data, aes(x = Feature, y = "PAM50")) +  
  geom_tile(aes(fill = PAM50)) +  # 使用因子Value  
  theme_minimal() +  
  scale_fill_manual(values = my_pal) +  # 手动设置颜色
  labs(x = NULL, y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",         
        panel.grid.major = element_blank(),  # 去除主要网格线           
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
temp_PAM50

temp_Stage <- ggplot(long_data, aes(x = Feature, y = "Stage")) +  
  geom_tile(aes(fill = Stage)) +  # 使用因子Value  
  theme_minimal() +  
  scale_fill_manual(values = my_pal) +  # 手动设置颜色
  labs(x = NULL, y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",         
        panel.grid.major = element_blank(),  # 去除主要网格线           
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
temp_Stage

temp_Grade <- ggplot(long_data, aes(x = Feature, y = "Grade")) +  
  geom_tile(aes(fill = Grade)) +  # 使用因子Value  
  theme_minimal() +  
  scale_fill_manual(values = my_pal) +  # 手动设置颜色
  labs(x = NULL, y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",         
        panel.grid.major = element_blank(),  # 去除主要网格线           
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
temp_Grade

temp_LNM <- ggplot(long_data, aes(x = Feature, y = "LNM.0.No.1.Yes.")) +  
  geom_tile(aes(fill = LNM.0.No.1.Yes.)) +  # 使用因子Value  
  theme_minimal() +  
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue")) +  # 手动设置颜色
  labs(x = NULL, y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",         
        panel.grid.major = element_blank(),  # 去除主要网格线           
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
temp_LNM

temp_PR <- ggplot(long_data, aes(x = Feature, y = "PR_stage")) +  
  geom_tile(aes(fill = PR_stage)) +  # 使用因子Value  
  theme_minimal() +  
  scale_fill_manual(values = my_pal) +  # 手动设置颜色
  labs(x = NULL, y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",
        panel.grid.major = element_blank(),  # 去除主要网格线           
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
temp_PR

temp_ER <- ggplot(long_data, aes(x = Feature, y = "ER_stage")) +  
  geom_tile(aes(fill = ER_stage)) +  # 使用因子Value  
  theme_minimal() +  
  scale_fill_manual(values = my_pal) +  # 手动设置颜色
  labs(x = NULL, y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",
        panel.grid.major = element_blank(),  # 去除主要网格线  
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
temp_ER

temp_HER2 <- ggplot(long_data, aes(x = Feature, y = "HER2_stage")) +  
  geom_tile(aes(fill = HER2_stage)) +  # 使用因子Value  
  theme_minimal() +  
  scale_fill_manual(values = my_pal) +  # 手动设置颜色
  labs(x = NULL, y =NULL) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_text(size = 8),  
        legend.position = "right",
        panel.grid.major = element_blank(),  # 去除主要网格线  
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))  # 去除次要网格线
temp_HER2

meta_data_heatmap <- temp_pheatmap %>%  
  aplot::insert_top(temp_age, height = 0.1) %>%  
  aplot::insert_top(temp_PAM50, height = 0.1) %>%  
  aplot::insert_top(temp_Stage, height = 0.1) %>%  
  aplot::insert_top(temp_Grade, height = 0.1) %>%  
  aplot::insert_top(temp_LNM, height = 0.1) %>%  
  aplot::insert_top(temp_PR, height = 0.1) %>%  
  aplot::insert_top(temp_ER, height = 0.1) %>%  
  aplot::insert_top(temp_HER2, height = 0.1) 
ggsave("results/figures/data_statistic/meta_data_raw_heatmap.png", plot = meta_data_heatmap, width = 12, height = 10, units = "in", dpi = 300)

remove(temp_data, annotation_col, long_data, temp_pheatmap, temp_age, temp_PAM50, temp_Stage, temp_Grade, temp_LNM, temp_PR, temp_ER, temp_HER2, meta_data_heatmap, meta_data_raw)


#稀疏曲线
# 计算每个样本的丰度  
temp_data <- colSums(tse_all@otu_table) %>% as.data.frame()   
colnames(temp_data) <- "Abundance"  

# 转换丰度为 numeric  
temp_data$Abundance <- as.numeric(as.character(temp_data$Abundance))  

# 定义特定深度  
specific_depths <- c(0, 500, 700, 900,1200,1500,3000,5000,7000)   

# 计算大于每个特定深度的样本数  
sample_counts <- sapply(specific_depths, function(depth) {  
  sum(temp_data$Abundance > depth)  
})  

# 创建数据框  
depth_counts <- data.frame(Depth = specific_depths, SampleCount = sample_counts)  

# 绘制图形  
temp_1 <- ggplot(depth_counts, aes(x = Depth, y = SampleCount)) +  
  geom_line() +  # 使用折线图  
  geom_point() +  # 添加点  
  labs(x = "Sequencing Depth (Abundance)", y = "Number of Samples ") +  
  theme_minimal()  
ggsave("results/figures/data_statistic/depth_counts.png", plot = temp_1, width = 4, height = 3, units = "in", dpi = 300)
#计算shannon稀疏曲线
# 计算每个样本的 Shannon 指数  
shannon_indices <- diversity(tse_all@otu_table, index = "shannon",MARGIN = 2)  

# 创建数据框  
# 初始化一个数据框来存储 Shannon 指数  
shannon_results <- colnames(tse_all@otu_table) %>% as.data.frame()
colnames(shannon_results) <- "Sample"
# 对于每个采样深度，计算每个样本的 Shannon 指数  
for (depth in specific_depths[-1]) {  
  # 随机抽样  
  sampled_data <- norm_rarefy(tse_all,size = depth)@otu_table %>% as.matrix() %>% as.data.frame()
  shannon_index <- diversity(sampled_data, index = "shannon", MARGIN = 2)  %>% as.data.frame()
  shannon_index$Sample <- rownames(shannon_index)
  colnames(shannon_index) <- c(depth, "Sample")
  # 将结果添加到总结果数据框
  shannon_results <- merge(shannon_results, shannon_index, by = "Sample", all = T)
}  

# 绘制 Shannon 稀疏曲线  
shannon_results <- shannon_results %>%  
  pivot_longer(cols = -Sample, names_to = "Depth", values_to = "ShannonIndex") %>%  
  mutate(Depth = as.numeric(as.character(Depth)))
# 将 Depth 列转换为因子，并设置水平顺序  
shannon_results$Depth <- factor(shannon_results$Depth, levels = specific_depths)  

# 绘制 Shannon 稀疏曲线  
ggplot(shannon_results, aes(x = Depth, y = ShannonIndex, group = Sample, color = Sample)) +  
  geom_line() +  # 使用折线图  
  geom_point() +  # 添加点  
  labs(x = "Sampling Depth", y = "Shannon Diversity Index") +  
  theme_minimal() +  
  theme(legend.position = "none")  # 如果样本数量较多，可以隐藏图例  
temp_2 <- ggplot(shannon_results, aes(x = Depth, y = ShannonIndex, group = Sample, color = Sample)) +  
  geom_line() +  # 使用折线图  
  geom_point() +  # 添加点  
  labs(x = "Sampling Depth", y = "Shannon Diversity Index") +  
  theme_minimal() +  
  theme(legend.position = "none")  # 如果样本数量较多，可以隐藏图例 
ggsave("results/figures/data_statistic/rarefaction_curve.png", plot = temp_2, width = 4, height = 23, units = "in", dpi = 300)

shannon_results$Group <- meta_data_metagenome$Group[match(shannon_results$Sample, meta_data_metagenome$sample)]
shannon_results$PAM50 <- meta_data_metagenome$PAM50[match(shannon_results$Sample, meta_data_metagenome$sample)]
shannon_results$stage <- meta_data_metagenome$stage[match(shannon_results$Sample, meta_data_metagenome$sample)]

# 绘制 Shannon 稀疏曲线 Group \ PAM50 \ stage
shannon_means <- shannon_results %>%  
  group_by(Group, Depth) %>%  
  summarise(Mean_Shannon = mean(ShannonIndex, na.rm = TRUE))  
Group_shannon <- ggplot(shannon_means, aes(x = Depth, y = Mean_Shannon,group = Group, color = Group)) +  
  geom_line() +  
  geom_point() +  
  labs(x = "Sampling Depth", y = "Mean Shannon Diversity Index") +  
  theme_minimal() +  
  theme(legend.position = "bottom")
Group_shannon
ggsave("results/figures/data_statistic/Group_shannon.png", plot = Group_shannon, width = 4, height = 3, units = "in", dpi = 300)

shannon_means <- shannon_results %>%  
  group_by(PAM50, Depth) %>%  
  summarise(Mean_Shannon = mean(ShannonIndex, na.rm = TRUE))
PAM50_shannon <- ggplot(shannon_means, aes(x = Depth, y = Mean_Shannon,group = PAM50, color = PAM50)) +
  geom_line() +
  geom_point() +
  labs(x = "Sampling Depth", y = "Mean Shannon Diversity Index") +
  theme_minimal() +
  theme(legend.position = "bottom")
PAM50_shannon
ggsave("results/figures/data_statistic/PAM50_shannon.png", plot = PAM50_shannon, width = 4, height = 3, units = "in", dpi = 300)

shannon_means <- shannon_results %>%  
  group_by(stage, Depth) %>%  
  summarise(Mean_Shannon = mean(ShannonIndex, na.rm = TRUE))
stage_shannon <- ggplot(shannon_means, aes(x = Depth, y = Mean_Shannon,group = stage, color = stage)) +
  geom_line() +
  geom_point() +
  labs(x = "Sampling Depth", y = "Mean Shannon Diversity Index") +
  theme_minimal() +
  theme(legend.position = "bottom")
stage_shannon
ggsave("results/figures/data_statistic/stage_shannon.png", plot = stage_shannon, width = 4, height = 3, units = "in", dpi = 300)
