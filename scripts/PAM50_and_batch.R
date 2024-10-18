# 读取数据
feature_counts <- read.table('data/all.sample.rawcount.all.txt',head=T,row.names=1)
count_data <- feature_counts[,-1]
# 处理列名
sample_id <- sub("X", "", colnames(count_data))
colnames(count_data) <- sample_id
meta_data <- read.csv("data/metadata.csv")
# 转置基因表达数据并生成注释
pam_50_data <- t(count_data)
dannot=data.frame(probe=colnames(pam_50_data),
                  "Gene.Symbol" =colnames(pam_50_data),
                  "EntrezGene.ID"=feature_counts$GeneID)
# 加载 PAM50 数据并执行分型
data(pam50.robust)
pam50_result<-molecular.subtyping(sbt.model = "pam50",data=pam_50_data,
                                  annot=dannot,do.mapping=TRUE)
PAM50_res <- as.data.frame(pam50_result$subtype)

# 将标准化后的结果转换为数据框  
PAM50_res_z <- as.data.frame(pam50_result_z$subtype)  

comparison <- data.frame(Original = PAM50_res$`pam50_result$subtype`, Z_Score = PAM50_res_z$`pam50_result_z$subtype`)  
table(comparison)  

# 添加 PAM50 结果到 meta_data
rownames(PAM50_res) <- gsub("R", "", rownames(PAM50_res))
meta_data$PAM50 <- PAM50_res$`pam50_result$subtype`[match(meta_data$T,rownames(PAM50_res))]

write.csv(meta_data,"data/metadata.csv",row.names = F)

meta_data_metagenome <- meta_data[!is.na(meta_data$WGS_T_tissue) | 
                                    !is.na(meta_data$WGS_N_tissue) | 
                                    !is.na(meta_data$WGS_N_blood), ]
meta_data_metagenome_1 <- meta_data_metagenome %>%
  mutate(sample = ifelse(!is.na(WGS_N_blood), WGS_N_blood, WGS_N_tissue))
meta_data_metagenome_2 <- meta_data_metagenome %>%
  mutate(sample = WGS_T_tissue)
meta_data_metagenome <- rbind(meta_data_metagenome_1, meta_data_metagenome_2)

batch_data <- readxl::read_xlsx("data/高深度WGS的任务单批次.xlsx")

meta_data_metagenome$batch <- batch_data$任务单号[match(meta_data_metagenome$sample, batch_data$样品名称)]

# 保存带批次信息的 meta_data_metagenome
write.csv(meta_data_metagenome,"data/metadata_metagenome.csv",row.names = F)
