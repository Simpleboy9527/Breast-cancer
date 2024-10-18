alpha_plot_fun <- function(data, meta_data, group_var = "group") { 
  data <- data %>%  
    mutate(group = meta_data[[group_var]][match(rownames(data), meta_data$sample)]) %>%  
    filter(!is.na(group))  
  
  # 定义组别信息  
  group_levels <- unique(data[["group"]])  
  
  # 提前定义显著性比较的空列表  
  significant_comparisons_observed <- NULL  
  significant_comparisons_shannon <- NULL  
  significant_comparisons_evenness <- NULL  
  
  # 如果组别数大于等于2，进行显著性检验  
  if (length(group_levels) >= 2) {  
    comparisons <- combn(group_levels, 2, simplify = FALSE)  
    
    # 定义显著性阈值  
    significance_level <- 0.05  
    
    # 计算显著性检验并过滤显著的比较  
    get_significant_comparisons <- function(variable) {  
      significant_comparisons <- lapply(comparisons, function(comp) {  
        group1_data <- data[data[["group"]] == comp[1], variable]  
        group2_data <- data[data[["group"]] == comp[2], variable]  
        
        # 检查是否有足够的数据进行检验  
        if (length(group1_data) > 1 && length(group2_data) > 1) {  
          test_result <- t.test(group1_data, group2_data)  
          if (test_result$p.value < significance_level) {  
            return(comp)  
          }  
        }  
        return(NULL)  
      })  
      Filter(Negate(is.null), significant_comparisons)  
    }  
    
    significant_comparisons_observed <- get_significant_comparisons("observed")  
    significant_comparisons_shannon <- get_significant_comparisons("diversity_shannon")  
    significant_comparisons_evenness <- get_significant_comparisons("evenness_pielou")  
  }  
  
  # 创建绘图函数  
  create_plot <- function(y_var, significant_comparisons = NULL) {  
    p <- ggplot(data, aes(x = group, y = .data[[y_var]])) +   
      geom_boxplot() +   
      theme_minimal() +
      geom_beeswarm(alpha = 0.5, cex = 0.8) +  
      xlab(NULL) +  
      ylab(y_var) +  
      ggtitle(y_var) +  
      theme(plot.title = element_text(hjust = 0.5),  
            axis.text.x = element_text(size = 12, family = "bold", angle = 45, hjust = 1))  
    
    # 如果有显著性比较，添加显著性标记  
    if (!is.null(significant_comparisons) && length(significant_comparisons) > 0) {  
      p <- p + geom_signif(comparisons = significant_comparisons,   
                           test = "t.test",   
                           map_signif_level = TRUE,   
                           step_increase = 0.1)  
    }  
    
    return(p)  
  }  
  
  # 生成所有多样性指数的图  
  observed_species_plot <- create_plot("observed", significant_comparisons_observed)  
  Shannon_plot <- create_plot("diversity_shannon", significant_comparisons_shannon)  
  Evenness_plot <- create_plot("evenness_pielou", significant_comparisons_evenness)  
  
  # 组合图  
  combined_plots <- observed_species_plot + Shannon_plot + Evenness_plot  
  return(combined_plots) 
}



#beta diversity
beta_diversity_fun <- function(micro_data, metadata){
  micro_data_prop <- micro_data / colSums(micro_data)
  micro_data_prop <- as.matrix(micro_data_prop)
  bray_curtis <- rbiom::beta.div(micro_data_prop, method = "bray")
  jaccard <- rbiom::beta.div(micro_data_prop, method = "jaccard", weighted = F)
  return(list(bray_curtis, jaccard))
}

beta_plot_fun <- function(distance_matrix, meta_data, variables_to_plot = c("Group", "age", "batch", "PAM50")) {
  # Perform PCoA
  data <- cmdscale(distance_matrix, k = 3, eig = TRUE)
  # Prepare PCoA points
  dune_pcoa_points <- as.data.frame(data$points)
  sum_eig <- sum(data$eig)
  eig_percent <- round(data$eig / sum_eig * 100, 1)
  colnames(dune_pcoa_points) <- paste0("PCoA", 1:3)
  
  
  # Initialize plot list
  plot_list <- list()
  
  # Function to calculate and plot each variable
  plot_variable <- function(var_name) {
    # Check if the variable has more than one level (otherwise PERMANOVA cannot be performed)
    if (length(unique(meta_data[[var_name]])) < 2) {
      message(paste("Variable", var_name, "has less than two levels, skipping..."))
      return(NULL)
    }
    
    # Run PERMANOVA
    adonis2_result <- vegan::adonis2(as.formula(paste0("distance_matrix ~ ", var_name)), 
                                     data = meta_data, na.action = na.omit)
    
    # Add grouping information to PCoA points
    dune_pcoa_points[[var_name]] <- meta_data[[var_name]][match(rownames(dune_pcoa_points), meta_data$sample)]
    dune_pcoa_points$Group <- meta_data$Group[match(rownames(dune_pcoa_points), meta_data$sample)]
    dune_pcoa_points$Group <- factor(dune_pcoa_points$Group, levels = c("Tumor_tissue", "Adjacent_tissue", "Blood"))
    
    # Plot
    plot <- ggplot(dune_pcoa_points, aes(x = PCoA1, y = PCoA2, color = !!sym(var_name), shape = Group)) +  
      geom_point(size = 3) +  
      xlab(paste0("PCoA1 (", eig_percent[1], "%)\n",   
                  "PERMANOVA R2: ", round(adonis2_result$R2[1], 3),   
                  "; P-value: ", adonis2_result$`Pr(>F)`[1])) +  
      ylab(paste0("PCoA2 (", eig_percent[2], "%)")) +  
      theme(  
        axis.title = element_text(size = 8, face = "bold"),  
        axis.text = element_text(size = 8, face = "bold"),  
        legend.text = element_text(size = 8, face = "bold"),  
        legend.title = element_text(size = 8, face = "bold"),  
        panel.border = element_rect(fill = NA, color = "black", linewidth = .5, linetype = "solid")  
      )  
    
    return(plot)
  }
  
  # Loop over specified variables to plot
  for (var in variables_to_plot) {
    plot <- plot_variable(var)
    if (!is.null(plot)) {
      plot_list[[var]] <- plot
    }
  }
  
  # Combine selected plots in one row
  p <- patchwork::wrap_plots(plot_list, ncol = length(plot_list))
  return(p)
}


#adnonis
adonis_fun <- function(distance_matrix, meta_data, selected_variables) {  
  adonis_summary <- data.frame(  
    Variable = character(),  
    R2 = numeric(),  
    PrF = numeric(),  
    stringsAsFactors = FALSE  
  )  
  
  for (i in selected_variables) {  
    if (i %in% colnames(meta_data)) {  
      # 检查列 i 中的因子水平数量  
      if (length(unique(meta_data[[i]])) > 1) {  
        distance_matrix <- as.matrix(distance_matrix)  
        matched_samples <- intersect(rownames(distance_matrix), meta_data$sample)  
        distance_matrix <- distance_matrix[matched_samples, matched_samples]  
        
        meta_data <- meta_data[meta_data$sample %in% matched_samples, ]  
        
        adonis_result <- vegan::adonis2(  
          as.formula(paste0("distance_matrix ~ ", i)),  
          data = meta_data, na.action = na.omit  
        )  
        
        # 检查结果是否包含有效的 R² 和 Pr(>F)  
        if (!is.null(adonis_result$R2[1]) && !is.null(adonis_result$`Pr(>F)`[1])) {  
          r2 <- adonis_result$R2[1]  # 提取第一个变量的R2  
          prf <- adonis_result$`Pr(>F)`[1]  # 提取第一个变量的Pr(>F)  
          
          # 将结果添加到数据框  
          adonis_summary <- rbind(adonis_summary, data.frame(  
            Variable = i,  
            R2 = r2,  
            PrF = prf  
          ))  
        } else {  
          message(paste("未能提取有效的R2或Pr(>F)值，跳过列:", i))  
        }  
      } else {  
        message(paste("跳过因子水平少于两个的列:", i))  
      }  
    } else {  
      message(paste("变量", i, "不在meta_data中，跳过..."))  
    }  
  }  
  
  return(adonis_summary)  
}

#diff_abundance
#phyloseq object build
phyloseq_creat <- function(micro_data, meta_data, taxa_are_rows = TRUE) {
  # 检查物种是行还是列
  if (taxa_are_rows) {
    taxa_data <- data.frame(taxonomy = rownames(micro_data), row.names = rownames(micro_data))
    taxa_data <- taxa_data %>%
      mutate(
        Domain = str_extract(taxonomy, "d__([^|]+)"),
        Phylum = str_extract(taxonomy, "p__([^|]+)"),
        Class = str_extract(taxonomy, "c__([^|]+)"),
        Order = str_extract(taxonomy, "o__([^|]+)"),
        Family = str_extract(taxonomy, "f__([^|]+)"),
        Genus = str_extract(taxonomy, "g__([^|]+)"),
        # 调整species的正则表达式，以匹配最后一部分直到字符串结束
        Species = str_extract(taxonomy, "s__([^|]*(?!.*s__))$")
      )
  } else {
    taxa_data <- data.frame(taxonomy = colnames(micro_data), row.names = colnames(micro_data))
    taxa_data <- taxa_data %>%
      mutate(
        Domain = str_extract(taxonomy, "d__([^|]+)"),
        Phylum = str_extract(taxonomy, "p__([^|]+)"),
        Class = str_extract(taxonomy, "c__([^|]+)"),
        Order = str_extract(taxonomy, "o__([^|]+)"),
        Family = str_extract(taxonomy, "f__([^|]+)"),
        Genus = str_extract(taxonomy, "g__([^|]+)"),
        # 调整species的正则表达式，以匹配最后一部分直到字符串结束
        Species = str_extract(taxonomy, "s__([^|]*(?!.*s__))$")
      )
  }
  taxa_data <- taxa_data[, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
  # 确保 OTU 名称一致
  if (taxa_are_rows) {
    otu_names <- rownames(micro_data)
  } else {
    otu_names <- colnames(micro_data)
  }
  
  taxa_names <- rownames(taxa_data)
  
  if (!all(otu_names == taxa_names)) {
    stop("OTU names in micro_data and taxa_data do not match.")
  }
  
  # 创建 phyloseq 对象
  rownames(meta_data) <- meta_data$sample
  otu_table <- otu_table(micro_data, taxa_are_rows = taxa_are_rows)
  tax_table <- tax_table(as.matrix(taxa_data))  # 强制转换为矩阵
  sample_data <- sample_data(meta_data)
  physeq_data <- phyloseq(otu_table, tax_table, sample_data)
  
  return(physeq_data)
}




#creat-heatmap-lm
create_heatmaps <- function(lm_summary, micro_data_tumor_tss, meta_data_tumor, baseline_var, q_threshold = 0.05, color_palette) {
  # 选择基准变量时的摘要
  lm_summary_select <- lm_summary[lm_summary$var_control == baseline_var,] 
  
  # 获取需要调整的列名
  adjust_cols <- setdiff(colnames(lm_summary_select), c("taxa", "var_control", paste0("PAM50", baseline_var), paste0("age:PAM50", baseline_var)))
  
  # 对这些列进行 p 值调整
  lm_summary_select <- lm_summary_select %>%
    mutate(across(all_of(adjust_cols), ~ p.adjust(., method = "BH"), .names = "q_{col}"))
  
  # 根据 q 值进行筛选，获取 q < 0.05 的行
  q_cols <- grep("^q_", colnames(lm_summary_select), value = TRUE)
  lm_summary_filtered <- lm_summary_select %>%
    filter(if_any(all_of(q_cols), ~ . < q_threshold))
  
  write.csv(lm_summary_filtered, file = paste0("results/", baseline_var, "_lm_summary_filtered.csv"), row.names = FALSE)
  
  # 创建 age_diff 热图
  micro_data_tumor_tss_filtered <- micro_data_tumor_tss[,colnames(micro_data_tumor_tss) %in% lm_summary_filtered$taxa[lm_summary_filtered$q_age < q_threshold]] %>% t() %>% as.data.frame()
  annotation_col <- meta_data_tumor[order(meta_data_tumor$age),c("age","PAM50","batch","age_group_55","age_group_45_60")]
  annotation_col <- annotation_col[!is.na(annotation_col$age),]
  micro_data_tumor_tss_filtered <- micro_data_tumor_tss_filtered[,rownames(annotation_col)]
  rownames(micro_data_tumor_tss_filtered) <- gsub(".*\\|g__","g__",rownames(micro_data_tumor_tss_filtered))
  age_heatmap <- pheatmap(micro_data_tumor_tss_filtered, cluster_rows = TRUE, cluster_cols = FALSE, 
                          annotation_col = annotation_col, show_colnames = FALSE, color = color_palette)
  ggsave(paste0("results/figures/", baseline_var, "_age_diff.png"), plot = age_heatmap, width = 8, height = 7, units = "in", dpi = 300)
  
  # 创建 PAM_diff 热图
  q_cols <- grep("^q_PAM50", colnames(lm_summary_filtered), value = TRUE)
  micro_data_tumor_tss_filtered <- data.frame()
  
  for (i in q_cols) {
    temp <- micro_data_tumor_tss[,colnames(micro_data_tumor_tss) %in% lm_summary_filtered$taxa[lm_summary_filtered[[i]] < q_threshold]] %>% t() %>% as.data.frame() 
    temp_char <- gsub("q_PAM50","",i)
    temp <- temp %>% mutate(Contrast = paste0(baseline_var, "_", temp_char))
    micro_data_tumor_tss_filtered <- rbind(micro_data_tumor_tss_filtered, temp)
  }
  
  annotation_row <- micro_data_tumor_tss_filtered$Contrast %>% as.data.frame()
  rownames(micro_data_tumor_tss_filtered) <- gsub(".*\\|g__","g__",rownames(micro_data_tumor_tss_filtered))
  rownames(annotation_row) <- rownames(micro_data_tumor_tss_filtered)
  colnames(annotation_row) <- "Contrast"
  
  annotation_col <- annotation_col[order(annotation_col$PAM50),]
  annotation_col <- annotation_col[!is.na(annotation_col$PAM50),]
  
  micro_data_tumor_tss_filtered <- micro_data_tumor_tss_filtered[,rownames(annotation_col)]
  pam_heatmap <- pheatmap(micro_data_tumor_tss_filtered, cluster_rows = TRUE, cluster_cols = FALSE, 
                          annotation_col = annotation_col, annotation_row = annotation_row,
                          show_colnames = FALSE, color = color_palette)
  ggsave(paste0("results/figures/", baseline_var, "_PAM_diff.png"), plot = pam_heatmap, width = 8, height = 7, units = "in", dpi = 300)
  
  # 创建 PAM_diff 的 log10 热图
  micro_data_tumor_tss_filtered <- log10(micro_data_tumor_tss_filtered + min(micro_data_tumor_tss_filtered[micro_data_tumor_tss_filtered != 0]))
  pam_heatmap_log10 <- pheatmap(micro_data_tumor_tss_filtered, cluster_rows = TRUE, cluster_cols = FALSE, 
                                annotation_col = annotation_col, annotation_row = annotation_row,
                                show_colnames = FALSE, color = color_palette)
  ggsave(paste0("results/figures/", baseline_var, "_PAM_diff_log10.png"), plot = pam_heatmap_log10, width = 8, height = 7, units = "in", dpi = 300)
  
  # 创建 age_PAM_diff 热图
  q_cols <- grep("^q_age_PAM50", colnames(lm_summary_filtered), value = TRUE)
  micro_data_tumor_tss_filtered <- data.frame()
  
  for (i in q_cols) {
    temp <- micro_data_tumor_tss[,colnames(micro_data_tumor_tss) %in% lm_summary_filtered$taxa[lm_summary_filtered[[i]] < q_threshold]] %>% t() %>% as.data.frame() 
    temp_char <- gsub("q_age_PAM50","",i)
    temp <- temp %>% mutate(Contrast = paste0(baseline_var, "_", temp_char))
    micro_data_tumor_tss_filtered <- rbind(micro_data_tumor_tss_filtered, temp)
  }
  
  annotation_row <- micro_data_tumor_tss_filtered$Contrast %>% as.data.frame()
  rownames(micro_data_tumor_tss_filtered) <- gsub(".*\\|s__","s__",rownames(micro_data_tumor_tss_filtered))
  rownames(annotation_row) <- rownames(micro_data_tumor_tss_filtered)
  colnames(annotation_row) <- "Contrast"
  
  micro_data_tumor_tss_filtered <- micro_data_tumor_tss_filtered[,rownames(annotation_col)]
  pam_heatmap <- pheatmap(micro_data_tumor_tss_filtered, cluster_rows = TRUE, cluster_cols = FALSE, 
                          annotation_col = annotation_col, annotation_row = annotation_row,
                          show_colnames = FALSE, color = color_palette)
  ggsave(paste0("results/figures/", baseline_var, "_age_PAM_diff.png"), plot = pam_heatmap, width = 8, height = 7, units = "in", dpi = 300)
  
  # 创建 age_PAM_diff 的 log10 热图
  micro_data_tumor_tss_filtered <- log10(micro_data_tumor_tss_filtered + min(micro_data_tumor_tss_filtered[micro_data_tumor_tss_filtered != 0]))
  pam_heatmap_log10 <- pheatmap(micro_data_tumor_tss_filtered, cluster_rows = TRUE, cluster_cols = FALSE, 
                                annotation_col = annotation_col, annotation_row = annotation_row,
                                show_colnames = FALSE, color = color_palette)
  ggsave(paste0("results/figures/", baseline_var, "_age_PAM_diff_log10.png"), plot = pam_heatmap_log10, width = 8, height = 7, units = "in", dpi = 300)

  }

#物种堆积图
taxa_stack_plot <- function(physeq_data, group_var = "Group", taxa_level = "Genus", top_n = 15,faced_group_ordered) {
  # 提取物种丰度数据
  physeq_data <- norm_tss(physeq_data)
  taxa_abundance <- physeq_data@otu_table %>%
    as.data.frame() 
 
  
  # 提取物种注释数据
  taxa_annotation <- physeq_data@tax_table %>%
    as.data.frame() 
  
  # 提取样本数据
  sample_data <- phyloseq::sample_data(physeq_data) %>%
    as.data.frame()
  
  # 提取分组变量
  group_levels <- unique(sample_data[[group_var]])
  
  # 提取 top_n 物种
  temp_taxa_abundance <- taxa_abundance 
  temp_taxa_abundance$taxa_special_level <- taxa_annotation[[taxa_level]][match(rownames(temp_taxa_abundance), rownames(taxa_annotation))]
  temp_taxa_abundance <- temp_taxa_abundance %>%
    group_by(taxa_special_level) %>%
    summarise_all(sum) 
  
  temp_taxa_abundance <- temp_taxa_abundance[order(rowSums(temp_taxa_abundance[,-1]),decreasing = T),] 
  temp_taxa_abundance$taxa_special_level[is.na(temp_taxa_abundance$taxa_special_level)] <- "Others"
  if (nrow(temp_taxa_abundance) > top_n) {
    temp_taxa_abundance$taxa_special_level[top_n : nrow(temp_taxa_abundance)] <- "Others"
  }
  temp_taxa_abundance <- temp_taxa_abundance %>%
    group_by(taxa_special_level) %>%
    summarise_all(sum) %>%
    as.data.frame()
  
  temp_taxa_abundance <- temp_taxa_abundance[order(rowSums(temp_taxa_abundance[,-1])),]
  factor_levels <- temp_taxa_abundance$taxa_special_level
  #others 放到最后
  if ("Others" %in% factor_levels) {
  others_index <- which(factor_levels == "Others")
  factor_levels <- c(rev(factor_levels[-others_index]), factor_levels[others_index])
  factor_levels <- gsub(".*__", "", factor_levels)
  }else{
    factor_levels <- gsub(".*__", "", rev(factor_levels))
  }
  
  #wide to long
  temp_taxa_abundance_long <- temp_taxa_abundance %>%
    pivot_longer(cols = -taxa_special_level, names_to = "sample", values_to = "abundance")
  temp_taxa_abundance_long$taxa_special_level <- gsub(".*__", "", temp_taxa_abundance_long$taxa_special_level)
  temp_taxa_abundance_long$taxa_special_level <- factor(temp_taxa_abundance_long$taxa_special_level, levels = factor_levels)
  temp_taxa_abundance_long$faced_group <-  sample_data[[group_var]][match(temp_taxa_abundance_long$sample, sample_data$sample)]
  temp_taxa_abundance_long <- temp_taxa_abundance_long[!is.na(temp_taxa_abundance_long$faced_group),]
  
  #sample_order
  # 提取 'Others' 以外的物种丰度数据
  tempdata <- temp_taxa_abundance[temp_taxa_abundance$taxa_special_level != "Others", ]
  
  # 找到丰度最高的物种行索引
  select_row_index <- which.max(rowSums(tempdata[,-1]))
  
  # 根据每个分组（faced_group）排序样本
  sample_order <- unlist(lapply(unique(temp_taxa_abundance_long$faced_group), function(f) {
    temp_data <- tempdata[, colnames(tempdata) %in% temp_taxa_abundance_long$sample[temp_taxa_abundance_long$faced_group == f], drop = FALSE]
    temp_order <- order(t(temp_data[select_row_index, ]), decreasing = TRUE)
    colnames(temp_data)[temp_order]
  }))
  
  # 对样本按排序后的顺序重新排列
  temp_taxa_abundance_long$sample <- factor(temp_taxa_abundance_long$sample, levels = sample_order)
  
  temp_taxa_abundance_long$faced_group <- factor(temp_taxa_abundance_long$faced_group, levels = faced_group_ordered,ordered = T)
  #plot
  bar_plot <- ggplot(temp_taxa_abundance_long, aes(x = sample, weight = abundance , fill = taxa_special_level)) +
    geom_bar(position = "stack", width = 1) +
    theme_bw() +
    scale_fill_manual(values = my_pal[1:(top_n+1)]) +
    facet_grid(~ faced_group, scales = "free_x",space='free') +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      strip.text.x = element_text(face = "bold", size = 8, color = "black", family = "Arial"),
      panel.grid = element_blank(),
      legend.position = "right",
      axis.title = element_text(face = "bold", size = 8, color = "black", family = "Arial"),
      legend.title = element_text(face = "bold", size = 8, color = "black", family = "Arial"),
      legend.text = element_text(face = "bold.italic", size = 8, color = "black", family = "Arial"),
    ) +
    xlab("Groups") +
    scale_y_continuous(
      name = "Relative abundance (%)",
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = paste(seq(0, 100, 25)),
      expand = c(0, 0)
    ) +
    labs(x = NULL, fill = taxa_level[1])

  return(bar_plot)
}
  
  
# 筛选微生物的函数
filter_microbes <- function(micro_data, abundance_threshold = 0, prevalence_threshold = 0) {
  # 计算每个微生物的丰度和流行度
  abundance <- rowMeans(micro_data, na.rm = TRUE) # 计算丰度
  prevalence <- rowSums(micro_data > 0, na.rm = TRUE) / ncol(micro_data) # 计算流行度
  
  # 筛选微生物
  filtered_microbes <- micro_data[abundance >= abundance_threshold & prevalence >= prevalence_threshold, ]
  
  return(filtered_microbes)
}

remove_non_functions <- function() {  
  # 获取当前环境中的所有对象名称  
  all_objects <- ls(envir = .GlobalEnv)  
  
  # 筛选出所有函数对象  
  functions <- sapply(all_objects, function(x) is.function(get(x, envir = .GlobalEnv)))  
  
  # 保留函数对象的名称  
  function_names <- all_objects[functions]  
  
  # 保留要保留的对象  
  keep_objects <- c("meta_data", "meta_data_metagenome", "meta_data_WGS_tumor",   
                    "meta_data_WGS_inter", "micro_data", "tse_all",   
                    "micro_data_tumor", "tse_tumor", "rna_data_TPM",   
                    "my_color_palette", "my_pal")  
  
  # 合并函数对象名称和要保留的对象  
  objects_to_keep <- union(function_names, keep_objects)  
  
  # 删除所有非保留对象  
  rm(list = setdiff(all_objects, objects_to_keep), envir = .GlobalEnv)  
  
  # 返回当前环境中的对象名称，供检查  
  return(ls(envir = .GlobalEnv))  
}  

