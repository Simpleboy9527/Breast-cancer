#export data to HAllA
# 获取属于 "Tumor_tissue" 组的样本名称  
tumor_samples <- meta_data_WGS_inter$sample[meta_data_WGS_inter$Group == "Tumor_tissue"]  

# 从 micro_data_tumor_tss 中选择与肿瘤样本名称匹配的列  
micro_data_tumor_inter <- micro_data_tumor_tss[rownames(micro_data_tumor_tss) %in% tumor_samples, ]  %>% t() %>% as.data.frame()

colnames(rna_data_TPM) <- gsub("R", "", colnames(rna_data_TPM))
rna_data_TPM_tumor <- rna_data_TPM[, colnames(rna_data_TPM) %in% tumor_samples]

# 保存数据
write.table(micro_data_tumor_inter, "data/micro_data_tumor_inter.txt", sep = "\t", quote = F)
write.table(rna_data_TPM_tumor, "data/rna_data_TPM_tumor.txt", sep = "\t", quote = F)

#读取hallA数据