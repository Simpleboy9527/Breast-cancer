#load metadata
meta_data <- read.csv("data/metadata_2.csv")
meta_data <- meta_data[!is.na(meta_data$age),]
meta_data <- meta_data %>% 
  mutate(
    age_group_45_60 = cut(age,breaks = c(0, 41, 60, Inf),labels = c("0~41","42~60",">60"))
  )
#样本筛选（去除metadata信息不足的样本）
meta_data <- meta_data[!is.na(meta_data$PAM50),]
meta_data <- meta_data[!is.na(meta_data$stage),]

meta_data_metagenome <- read.csv("data/metadata_metagenome.csv")
#add group column
meta_data_metagenome <- meta_data_metagenome %>%
  mutate(Group = case_when(
    grepl("T$", meta_data_metagenome$sample) ~ "Tumor_tissue",
    grepl("N$", meta_data_metagenome$sample) ~ "Adjacent_tissue",
    TRUE ~ "Blood"
  ),
  age_group_55 = cut(age,breaks = c(0, 55, Inf), labels = c("<55", ">=55")),
  age_group_45_60 = cut(age,breaks = c(0, 41, 60, Inf),labels = c("0~41","42~60",">60"))
  )
meta_data_metagenome <- meta_data_metagenome[!is.na(meta_data_metagenome$PAM50),]
meta_data_metagenome <- meta_data_metagenome[!is.na(meta_data_metagenome$age),]


meta_data_WGS_tumor <- meta_data_metagenome %>%
  filter(Group == "Tumor_tissue")
meta_data_WGS_tumor <- meta_data_WGS_tumor[!(is.na(meta_data_WGS_tumor$age) & is.na(meta_data_WGS_tumor$PAM50)),]

meta_data_WGS_inter <- meta_data_metagenome[(!is.na(meta_data_metagenome$WGBS_T_tissue) | 
                                           !is.na(meta_data_metagenome$WGBS_N_tissue) | 
                                           !is.na(meta_data_metagenome$WGBS_N_blood))&(
                                             !is.na(meta_data_metagenome$RNA_NR_Tissue)|
                                               !is.na(meta_data_metagenome$RNA_TR_Tissue)),] #选取有WGBS和RNA-seq数据的样本





#load data
##load mapped reads data
map_reads <- read_excel("data/combine_reads.xlsx")
##load total reads data
total_reads <- read_xlsx("data/WGS_RNAseq_QC.xlsx")


#load microbiome data
micro_data <- read.csv("data/microdata_Not_rmCT.txt",sep = "\t",row.names = 1,check.names = F)
colnames(micro_data) <- gsub("_report_bracken_species.txt", "",colnames(micro_data))
micro_data <- micro_data[rownames(micro_data) > 20 ,] #filter Kraken reads > 20

rownames(meta_data_metagenome) <- meta_data_metagenome$sample

tse_all <- phyloseq_creat(micro_data, meta_data_metagenome)

meta_data_metagenome_temp <- tse_all@sam_data 
meta_data_metagenome <- as.matrix(meta_data_metagenome_temp) %>% as.data.frame()
meta_data_metagenome$age <- as.numeric(meta_data_metagenome$age)
meta_data_metagenome$OS <- as.numeric(meta_data_metagenome$OS)

micro_data_tumor <- micro_data[,colnames(micro_data) %in% meta_data_WGS_tumor$sample]

tse_tumor <- phyloseq_creat(micro_data_tumor,meta_data_WGS_tumor,taxa_are_rows = T)
meta_data_WGS_tumor_temp <- tse_tumor@sam_data
meta_data_WGS_tumor <- as.matrix(meta_data_WGS_tumor_temp) %>% as.data.frame()
meta_data_WGS_tumor$age <- as.numeric(meta_data_WGS_tumor$age)
meta_data_WGS_tumor$OS <- as.numeric(meta_data_WGS_tumor$OS)

#load rna-seq data tpm
rna_data_TPM <- read.table("data/all.sample.TPM.all.txt", header = T, row.names = 1)
colnames(rna_data_TPM) <- gsub("X", "", colnames(rna_data_TPM))
rna_data_TPM <- rna_data_TPM[,-1]


#color_data
my_color_palette <- colorRampPalette(brewer.pal(9, "OrRd"))
col = pal_d3("category20")(20)
col2 = pal_d3("category20", alpha = 0.5)(20)
my_pal = c(col, col2[-8])

