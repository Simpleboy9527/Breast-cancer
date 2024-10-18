micro_data$species <- gsub(".*s__", "", rownames(micro_data))
relative_abundance$species <- gsub(".*s__", "", rownames(relative_abundance))
RNA_seq_micro <- read.csv("data/merged_output_RNA_seq.csv", row.names = 1)
rownames(RNA_seq_micro) <- gsub("\\ ", "_", rownames(RNA_seq_micro))
micro_data_clean <- micro_data[micro_data$species %in% rownames(RNA_seq_micro),] #filter species in RNA-seq data
micro_data_clean <- micro_data_clean[, -ncol(micro_data_clean)]  # 去除最后一列的species列
#统计去除的微生物情况
removed_species <- relative_abundance[!relative_abundance$species %in% rownames(RNA_seq_micro),]  # 去除RNA-seq数据中的物种 
removed_species <- removed_species[, -ncol(removed_species)]  # 去除最后一列的species列
rownames(removed_species) <- gsub(".*s__", "", rownames(removed_species))  # 提取物种名称
# 统计去除的微生物数量  
removed_abundance_sum <- rowSums(removed_species, na.rm = TRUE)  

# 计算丰度平均值  
removed_abundance_mean <- rowMeans(removed_species, na.rm = TRUE)  

# 将结果整理到数据框中  
removed_summary <- data.frame(  
  Species = gsub(".*s__", "", rownames(removed_species)),  # 提取物种名称  
  Total_Abundance = removed_abundance_sum,  
  Average_Abundance = removed_abundance_mean  
)   
#plot
# 选择丰度最高的前10个物种  
top_n_species <- removed_summary[order(-removed_summary$Total_Abundance), ][1:10, ]  

top_n_species_data <- removed_species[rownames(removed_species) %in% top_n_species$Species, ]  
top_n_species_data$species <- rownames(top_n_species_data)
# 将数据转换为长格式  
top_n_species_long <- melt(top_n_species_data, varnames = c("Species", "Sample"), value.name = "Abundance")

# 创建箱型图  
removed_plot <- ggplot(top_n_species_long, aes(x = Abundance, y = species)) +  
  geom_boxplot(fill = "steelblue") +  
  labs(title = NULL, x = "Abundance", y = NULL) +  
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) +  # y轴标签不旋转  
  theme_minimal()  
ggsave("results/figures/removed_plot.png", removed_plot, width = 5, height = 5, dpi = 300)

