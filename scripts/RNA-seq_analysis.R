#load RNA-seq data (TPM)
TPM_data <- read.table("data/all.sample.TPM.all.txt", header = T, row.names = 1)
RNA_count_data <- read.table("data/all.sample.rawcount.all.txt", header = T, row.names = 1)
colnames(TPM_data) <- gsub("X", "", colnames(TPM_data))
colnames(RNA_count_data) <- gsub("X", "", colnames(RNA_count_data))
#select tumor samples 
TPM_data_tumor <- TPM_data[, colnames(TPM_data) %in% meta_data_WGS_inter$RNA_TR_Tissue]
TPM_data_tumor$GeneID <- TPM_data$GeneID[match(rownames(TPM_data_tumor), rownames(TPM_data))]
RNA_count_data_tumor <- RNA_count_data[, colnames(RNA_count_data) %in% meta_data_WGS_inter$RNA_TR_Tissue]


#PCA analysis
# 进行 PCA 分析  
temp_data <- TPM_data_tumor %>% t() %>% as.data.frame()
temp_data <- temp_data[-nrow(temp_data),]
temp_data <- temp_data[, apply(temp_data, 2, function(col) var(col) != 0)]
pca_result <- prcomp(temp_data, center = TRUE, scale. = TRUE)  

# 查看 PCA 结果  
summary(pca_result)  

# 提取主成分  
pca_data <- as.data.frame(pca_result$x)  
pca_data$sample <- gsub("R","",rownames(pca_data))
pca_data <- merge(pca_data, meta_data_WGS_inter, by = "sample", all = F)
# 绘制 PCA 图  
ggplot(pca_data, aes(x = PC1, y = PC2, color = age, shape = PAM50)) +  
  geom_point(size = 3) +  
  xlab("PC1") +  
  ylab("PC2") +  
  theme_minimal() +
  theme(  
    axis.title = element_text(size = 8, face = "bold"),  
    axis.text = element_text(size = 8, face = "bold"),  
    legend.text = element_text(size = 8, face = "bold"),  
    legend.title = element_text(size = 8, face = "bold"),  
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid") 
  )  


#DA analysis for age_group
RNA_seq_meta_data <- meta_data_WGS_inter[paste(meta_data_WGS_inter$sample,"R",sep = "") %in% colnames(TPM_data_tumor),]
rownames(RNA_seq_meta_data) <- paste(RNA_seq_meta_data$sample, "R", sep = "")
RNA_seq_meta_data <- RNA_seq_meta_data[colnames(RNA_count_data_tumor),]
RNA_seq_meta_data <- RNA_seq_meta_data[!is.na(RNA_seq_meta_data$age_group_45_60),]
RNA_count_data_tumor <- RNA_count_data_tumor[,rownames(RNA_seq_meta_data)]
RNA_count_data_tumor <- round(RNA_count_data_tumor,0)


coldata <- RNA_seq_meta_data[,c("age_group_45_60","PAM50","Grade")]
coldata$age_group_45_60 <- as.factor(coldata$age_group_45_60)
coldata$PAM50 <- as.factor(coldata$PAM50)
coldata$Grade <- as.factor(coldata$Grade)


dss <- DESeqDataSetFromMatrix(countData = RNA_count_data_tumor, colData = coldata, design = ~ Grade + PAM50 + age_group_45_60)
dss <- DESeq(dss)
res <- results(dss,contrast = c('age_group_45_60', '42~60','0~41')) 
res_2 <- results(dss,contrast = c('age_group_45_60', '>60','0~41'))
res <- res[order(res$padj),]
write.csv(res, "data/DESeq2_age_group_45_60.csv")

#diff gene correlation analysis with microbiome
res1 <- res %>% as.data.frame() %>% na.omit()
res1[which(res$log2FoldChange >= 1 & res$padj < 0.05),'sig'] <- 'up'
res1[which(res$log2FoldChange <= -1 & res$padj < 0.05),'sig'] <- 'down'
res1[which(abs(res$log2FoldChange) <= 1 | res$padj >= 0.01),'sig'] <- 'none'
write.csv(res1, "data/DESeq2_age_group_45_60_sig.csv")
#火山图
p_volcano <- ggplot(res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +  
  geom_point(size = 2) +  
  scale_color_manual(values = c("red", "black", "blue")) +  
  theme_classic() +  
  labs(x = "log2FoldChange", y = "-log10(padj)") +  
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "red") +  
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +  
  theme(legend.position = "none")

p_volcano
ggsave("results/figures/volcano_age_group_45_60.png", p_volcano, width = 3, height = 3, units = "in")

#筛选差异基因
sig_gene <- rownames(res1)[which(res1$sig != 'none')]
TPM_data_tumor_selected <- TPM_data_tumor[rownames(TPM_data_tumor) %in% sig_gene,] %>% t() %>% as.data.frame()
TPM_data_tumor_selected$sample <- gsub("R", "", rownames(TPM_data_tumor_selected))
#差异基因与微生物组成的相关性分析
#合并RNA-seq数据和meta数据、差异分析结果
temp_data <- filtered_micro_data_age %>% t() %>% as.data.frame()
temp_data$sample <- rownames(temp_data)
merge_data_temp <- merge(temp_data, TPM_data_tumor_selected, by = "sample",all = F)


# 检查并转换非数值列
merge_data_temp[, 2:ncol(merge_data_temp)] <- lapply(merge_data_temp[, 2:ncol(merge_data_temp)], function(x) {
  if (!is.numeric(x)) as.numeric(as.character(x)) else x
})

cor_matrix_all <- cor(merge_data_temp[, 25:ncol(merge_data_temp)], merge_data_temp[, 2:24], method = "spearman", use = "pairwise.complete.obs")  

# 创建一个数据框用于存储结果  
results <- data.frame(  
  Variable1 = rep(colnames(cor_matrix_all), each = nrow(cor_matrix_all)),  
  Variable2 = rep(rownames(cor_matrix_all), times = ncol(cor_matrix_all)),  
  Correlation = as.vector(cor_matrix_all),  
  P_Value = NA,  
  R2 = NA  
)  

# 循环遍历所有列  
for (i in 1:ncol(cor_matrix_all)) {  
  # 循环遍历所有行  
  for (j in 1:nrow(cor_matrix_all)) {  
    # 获取merged_data中与cor_all_test的列名相匹配的列索引  
    a <- which(colnames(merge_data_temp) == colnames(cor_matrix_all)[i])  
    # 获取merged_data中与cor_all_test的行名相匹配的列索引  
    b <- which(colnames(merge_data_temp) == rownames(cor_matrix_all)[j])  
    
    # 使用spearman方法计算merged_data中第a列和第b列的相关性,并提取p值  
    test_result <- cor.test(merge_data_temp[, a], merge_data_temp[, b], method = "spearman")  
    results$P_Value[(i - 1) * nrow(cor_matrix_all) + j] <- test_result$p.value  
    results$R2[(i - 1) * nrow(cor_matrix_all) + j] <- test_result$estimate^2  
  }  
}  

results$padj <- p.adjust(results$P_Value, method = "BH")
results <- results %>%  
  mutate(text = case_when(  
    padj < 0.001 ~ "***",  
    padj < 0.01 ~ "**",  
    padj < 0.05 ~ "*",  
    TRUE ~ ""  
  ))
 

dist_matrix <- vegdist(t(cor_matrix_all), method = "euclidean")  # 确保 cor_matrix_all 是有效的相关性矩阵  
hclust_col <- hclust(dist_matrix, method = "average")
dist_matrix_col <- vegdist(cor_matrix_all, method = "euclidean")  # 确保 cor_matrix_all 是有效的相关性矩阵
hclust_row <- hclust(dist_matrix_col, method = "average")  

# 重新排序相关性矩阵  
cor_matrix_ordered <- cor_matrix_all[hclust_row$order,hclust_col$order]  


# 因子级别设置  
results$Variable1 <- factor(results$Variable1, levels = colnames(cor_matrix_ordered) , ordered = TRUE)  
results$Variable2 <- factor(results$Variable2, levels = rownames(cor_matrix_ordered), ordered = TRUE)




black_transparent <- rgb(0, 0, 0, alpha = 0.7)  
p_cor_sig <- ggplot(results, aes(x = Variable2, y = Variable1)) +  
  theme_classic() +  
  scale_fill_gradientn(colors = colorRampPalette(c("navy", "white", "firebrick3"))(100)) +  
  geom_tile(aes(fill = Correlation), color = black_transparent) +  
  labs(fill = NULL) +  
  geom_text(aes(label = text), size = 3, color = "black") +  
  theme(  
    axis.text.x = element_text(size = 9, face = "bold", angle = 30, hjust = 1.02, vjust = 1.1),
    axis.ticks.x = element_blank(),  
    axis.line.x = element_blank(),  
    axis.title = element_blank(),  
    axis.text.y = element_text(face = "bold.italic", size = 9),  
    axis.ticks.y = element_blank(),  
    legend.text = element_text(size = 9, color = "black", face = "bold"),  
    legend.key.width = unit(0.15, "in"),  
    legend.key.height = unit(0.2, "in"),  
    legend.box.margin = margin(l = 0.1, unit = "in"),  
    axis.line = element_blank()  
  )  
p_cor_sig
ggsave("results/figures/cor_sig_2.png", p_cor_sig, width = 15, height = 4.5, units = "in",limitsize = FALSE)




#gene表达热图和微生物热图
# 提取基因数据  
gene_data <- merge_data_temp[, c(1, 25:ncol(merge_data_temp))]  
rownames(gene_data) <- gene_data$sample  
gene_data <- gene_data[, -1]  

# 提取微生物数据  
micro_data_temp <- merge_data_temp[, 1:24]  
rownames(micro_data_temp) <- micro_data_temp$sample  
micro_data_temp <- micro_data_temp[, -1]  

# 匹配注释列  
annotation_columns <- meta_data_WGS_tumor[match(merge_data_temp$sample, meta_data_WGS_tumor$sample),   
                                          c("age_group_45_60", "PAM50", "age")]  
rownames(annotation_columns) <- merge_data_temp$sample  

# 按年龄排序并去除缺失值  
annotation_columns <- annotation_columns[order(annotation_columns$age), ]  
annotation_columns <- na.omit(annotation_columns)  

# 根据注释列筛选基因和微生物数据  
gene_data <- gene_data[rownames(annotation_columns), ]  
micro_data_temp <- micro_data_temp[rownames(annotation_columns), ]  

# 转置微生物数据并进行log10转换  
micro_data_temp <- t(micro_data_temp) %>% as.data.frame()  
micro_data_temp <- log10(min(micro_data_temp[micro_data_temp != 0]) + micro_data_temp)  
micro_data_temp$Microbe <- rownames(micro_data_temp)
# 转置基因数据并进行log10转换  
gene_data <- gene_data %>% t() %>% as.data.frame()  
gene_data <- log10(min(gene_data[gene_data != 0]) + gene_data)  
gene_data$Gene <- rownames(gene_data)
#heatmap - ggplot2
# 创建一个数据框用于存储数据
long_micro_data <- melt(micro_data_temp, varnames = c("Microbe"), value.name = "log10 Abundance")
long_gene_data <- melt(gene_data, varnames = c("Gene"), value.name = "log10 TPM")

#设置顺序
long_micro_data$variable <- factor(long_micro_data$variable, levels = rownames(annotation_columns), ordered = TRUE)
long_gene_data$variable <- factor(long_gene_data$variable, levels = rownames(annotation_columns), ordered = TRUE)
#根据聚类设置微生物和基因的顺序
# 计算距离矩阵和层次聚类
dist_matrix_micro <- vegdist(micro_data_temp[,-ncol(micro_data_temp)], method = "euclidean")
hclust_result_micro <- hclust(dist_matrix_micro, method = "average")
dist_matrix_gene <- vegdist(gene_data[,-ncol(gene_data)], method = "euclidean")
hclust_result_gene <- hclust(dist_matrix_gene, method = "average")

# 重新排序相关性矩阵
micro_data_temp<- micro_data_temp[hclust_result_micro$order,]
gene_data_ordered <- gene_data[hclust_result_gene$order,]

long_micro_data$Microbe <- factor(long_micro_data$Microbe, levels = rownames(micro_data_temp), ordered = TRUE)
long_gene_data$Gene <- factor(long_gene_data$Gene, levels = rownames(gene_data_ordered), ordered = TRUE)

# 绘制热图
annotation_columns$sample <- rownames(annotation_columns)
annotation_columns$sample <- factor(annotation_columns$sample, levels = rownames(annotation_columns), ordered = TRUE)
age_group_45_60 <- ggplot(annotation_columns, aes(x = sample, y = "Age_group_45_60")) +  
  geom_tile(aes(fill = age_group_45_60)) +  
  scale_fill_manual(values = c("0~41" = "#1f77b4", "42~60" = "#2ca02c", ">60" = "#d62728")) +  
  theme_minimal() +  
  theme(  
    axis.text.x = element_blank(),  
    axis.title = element_blank(),  
    axis.ticks = element_blank(),  
    axis.title.y.right = element_blank(),  
    axis.text.y = element_text(face = "bold", size = 9),
    legend.title = element_text(size = 9, hjust = 0.5),  
    panel.spacing = unit(0, "lines"),  # 设置面板间距为0  
    panel.border = element_blank()  # 去掉面板边框
    ) +  
  scale_y_discrete(position = "right")  
age_group_45_60
age_p <- ggplot(annotation_columns, aes(x = sample, y = "Age")) +  
  geom_tile(aes(fill = age)) +  
  theme_minimal() +  
  theme(  
    axis.text.x = element_blank(),  
    axis.ticks = element_blank(),  
    axis.title = element_blank(),  
    axis.ticks.y = element_blank(),  
    axis.text.y = element_text(face = "bold", size = 9),
    axis.title.y.right = element_blank(),  
    panel.spacing = unit(0, "lines")  # 设置面板间距为0  
  ) +  
  scale_y_discrete(position = "right")  

PAM50_p <- ggplot(annotation_columns, aes(x = sample, y = "PAM50")) +  
  geom_tile(aes(fill = PAM50)) +  
  scale_fill_manual(values = my_pal) +  
  theme_minimal() +  
  theme(  
    axis.text.x = element_blank(),  
    axis.ticks = element_blank(),  
    axis.title = element_blank(),  
    axis.line = element_blank(),  
    axis.text.y = element_text(face = "bold", size = 9),
    panel.spacing = unit(0, "lines")  # 设置面板间距为0  
  ) +  
  scale_y_discrete(position = "right")  

# 绘制热图  
micro_temp_p <- ggplot(long_micro_data, aes(x = variable, y = Microbe)) +  
  theme_classic() +  
  scale_fill_gradientn(colors = my_color_palette(100)) +  
  geom_tile(aes(fill = `log10 Abundance`)) +  
  theme(  
    axis.text.x = element_blank(),  
    axis.ticks = element_blank(),  
    axis.title = element_blank(),  
    axis.text.y = element_text(face = "bold.italic", size = 9),  
    legend.title = element_text(size = 9, hjust = 0.5),  
    axis.line = element_blank(),  
    panel.spacing = unit(0, "lines")  # 设置面板间距为0  
  ) +   
  scale_y_discrete(position = "right")  

# 合并图形  
diff_gene_micro_p <- micro_temp_p %>%  
  aplot::insert_top(age_group_45_60, height = 0.05) %>%  
  aplot::insert_top(age_p, height = 0.05) %>%  
  aplot::insert_top(PAM50_p, height = 0.05)  
diff_gene_micro_p

#绘制基因热图
gene_temp_p <- ggplot(long_gene_data, aes(x = variable, y = Gene)) +  
  theme_classic() +  
  scale_fill_gradientn(colors = my_color_palette(100)) +  
  geom_tile(aes(fill = `log10 TPM`)) +  
  theme(  
    axis.text.x = element_blank(),  
    axis.ticks = element_blank(),  
    axis.title = element_blank(),  
    axis.text.y = element_text(face = "bold.italic", size = 9),  
    legend.title = element_text(size = 9, hjust = 0.5),  
    axis.line = element_blank(),  
    panel.spacing = unit(0, "lines")  # 设置面板间距为0  
  ) +   
  scale_y_discrete(position = "right")

# 合并图形
diff_gene_micro_p <- diff_gene_micro_p %>%  
  aplot::insert_bottom(gene_temp_p, height = 1)
ggsave("results/figures/diff_gene_micro_p.png", diff_gene_micro_p, width = 10, height = 8, units = "in")



#文献记载人类与细胞增殖相关的gene
gene_list_cell_Proliferation <- data.frame(  
  Gene_Symbol = c("CCND1", "CDK4", "TP53", "RB1", "MYC", "PTEN", "AKT1", "ERK1/2", "CDK2", "CCNE1",   
                  "MDM2", "BCL2", "FOS", "JUN", "RAC1", "STAT3", "WNT", "BRCA1", "BRCA2", "ERBB2",   
                  "ESR1", "PGR"),  
  Entrez_Gene_ID = c(595, 983, 7157, 5925, 4609, 5728, 207, 5594, 1017, 891,   
                     4193, 596, 2353, 3725, 5870, 6772, 7460, 672, 675, 2064,   
                     2099, 5241)  
)  


TPM_data_tumor_selected <- TPM_data_tumor[TPM_data_tumor$GeneID %in% gene_list_cell_Proliferation$Entrez_Gene_ID,] %>% t() %>% as.data.frame()
TPM_data_tumor_selected$sample <- gsub("R", "", rownames(TPM_data_tumor_selected))

#合并RNA-seq数据和meta数据、差异分析结果
temp_data <- filtered_micro_data_age %>% t() %>% as.data.frame()
temp_data$sample <- rownames(temp_data)
merge_data_temp <- merge(temp_data, TPM_data_tumor_selected, by = "sample",all = F) 

# 假设 merge_data_temp 已经存在  
# 检查并转换非数值列  
merge_data_temp[, 2:45] <- lapply(merge_data_temp[, 2:45], function(x) {  
  if (!is.numeric(x)) as.numeric(as.character(x)) else x  
})  

# 计算相关性矩阵  
cor_matrix_all <- cor(merge_data_temp[, 25:45], merge_data_temp[, 2:24], method = "spearman", use = "pairwise.complete.obs")  

cor_all <- as.data.frame(cor_matrix_all)

# 创建一个测试数据框用于存储p值
cor_all_test <- cor_all
# 创建一个测试数据框用于存储R^2值
cor_all_test_R2 <- cor_all

# 循环遍历所有列
for (i in 1:ncol(cor_all)) {
  # 循环遍历所有行
  for (j in 1:nrow(cor_all)) {
    # 获取merged_data中与cor_all_test的列名相匹配的列索引
    a <- which(colnames(merge_data_temp) == colnames(cor_all_test)[i])
    # 获取merged_data中与cor_all_test的行名相匹配的列索引
    b <- which(colnames(merge_data_temp) == rownames(cor_all_test)[j])
    
    # 打印列索引a
    print(a)
    # 打印行索引b
    print(b)
    
    # 使用spearman方法计算merged_data中第a列和第b列的相关性,并提取p值
    p_value <- cor.test(merge_data_temp[, a], merge_data_temp[, b], method = "spearman")$p.value
    # 使用spearman方法计算merged_data中第a列和第b列的相关性,并提取R^2值
    R2 <- cor.test(merge_data_temp[, a], merge_data_temp[, b], method = "spearman")$estimate^2
    
    # 将p值存储到cor_all_test的对应位置
    cor_all_test[j, i] <- p_value
    # 将R^2值存储到cor_all_test_R2的对应位置
    cor_all_test_R2[j, i] <- R2
  }
}

# 转换为数据框  
cor_all <- as.data.frame(cor_all)  
cor_all_test <- as.data.frame(cor_all_test)  
cor_all_test_R2 <- as.data.frame(cor_all_test_R2)  

# 添加因子列  
cor_all$factor <- rownames(cor_all)  
cor_all_test$factor <- rownames(cor_all_test)  
cor_all_test_R2$factor <- rownames(cor_all_test_R2)  

# 合并数据  
long_matrix_all <- melt(cor_all, id.vars = "factor")  
# 使用 match 来获取 p_value 和 R2  
p_values <- melt(cor_all_test, id.vars = "factor")
R2_values <- melt(cor_all_test_R2, id.vars = "factor")  

# 将 p_value 和 R2 添加到 long_matrix_all  
long_matrix_all <- merge(long_matrix_all, p_values, by = c("factor", "variable"), suffixes = c("", "_p"))  
long_matrix_all <- merge(long_matrix_all, R2_values, by = c("factor", "variable"), suffixes = c("", "_R2"))  

# 重命名列  
colnames(long_matrix_all)[3] <- "correlation"  # 将相关性列重命名  
colnames(long_matrix_all)[4] <- "p_value"      # 将p值列重命名  
colnames(long_matrix_all)[5] <- "R2"           # 将R^2列重命名  

# 调整p值和标记  
long_matrix_all$padj <- p.adjust(long_matrix_all$p_value, method = "BH")  
long_matrix_all <- long_matrix_all %>%  
  mutate(text = case_when(  
    padj < 0.001 ~ "***",  
    padj < 0.01 ~ "**",  
    padj < 0.05 ~ "*",  
    TRUE ~ ""  
  ))  

# 计算距离矩阵和层次聚类  
dist_matrix <- vegdist(t(cor_matrix_all), method = "euclidean")  
hclust_result <- hclust(dist_matrix, method = "average")  

# 绘制树状图  
plot(hclust_result, main = "Complete-linkage Hierarchical Clustering Dendrogram")  

# 重新排序相关性矩阵  
cor_all <- cor_all[, hclust_result$order]  
temp_order <- colnames(cor_all)  

# 因子级别设置  
long_matrix_all$factor <- gsub("\\.", "-", long_matrix_all$factor)  
long_matrix_all$factor <- factor(long_matrix_all$factor, levels = unique(long_matrix_all$factor), ordered = TRUE)  
long_matrix_all$variable <- factor(long_matrix_all$variable, levels = temp_order, ordered = TRUE)  

# 绘制热图  
black_transparent <- rgb(0, 0, 0, alpha = 0.7)  
p_cor_all <- ggplot(long_matrix_all, aes(x = factor, y = variable)) +  
  theme_classic() +  
  scale_fill_gradientn(colors = colorRampPalette(c("navy", "white", "firebrick3"))(100)) +  
  geom_tile(aes(fill = correlation), color = black_transparent) +  
  labs(fill = NULL) +  
  geom_text(aes(label = text), size = 3, color = "black") +  
  theme(  
    axis.text.x = element_text(face = "bold", size = 9, angle = 30, hjust = 0.8, vjust = 1),  
    axis.ticks.x = element_blank(),  
    axis.line.x = element_blank(),  
    axis.title.y = element_blank(),  
    axis.text.y = element_text(face = "bold.italic", size = 9),  
    axis.ticks.y = element_blank(),  
    legend.text = element_text(size = 9, color = "black", face = "bold"),  
    legend.key.width = unit(0.15, "in"),  
    legend.key.height = unit(0.2, "in"),  
    legend.box.margin = margin(l = 0.1, unit = "in"),  
    axis.line = element_blank()  
  ) +
  scale_y_discrete(position = "right")

# 保存图形  
ggsave("results/figures/cor_cell_Proliferation.png", p_cor_all, width = 8, height = 3, units = "in")

#散点图
cor_test_result <- cor.test(merge_data_temp$`s__Sphingopyxis_sp._USTB-05`, merge_data_temp$TP53)  
cor_coefficient <- round(cor_test_result$estimate, 2)  # 相关系数  
p_value <- cor_test_result$p.value           # p值  

# 创建散点图  
p_cor_scatter <- ggplot(merge_data_temp, aes(x = `s__Sphingopyxis_sp._USTB-05`, y = TP53)) +  
  geom_point(size = 3) +  
  geom_smooth(method = "lm", se = TRUE) +  
  theme_classic() +  
  theme(legend.position = "none") +  
  labs(x = "s__Sphingopyxis_sp._USTB-05", y = "TPM_TP53") +  
  annotate("text", x = Inf, y = Inf,   
           label = paste("Correlation:", cor_coefficient, "\nP-value:", p_value),   
           hjust = 1.1, vjust = 1.5, size = 4, color = "black")  

# 显示图形  
p_cor_scatter
ggsave("results/figures/cor_scatter.png", p_cor_scatter, width = 4, height = 4, units = "in")

