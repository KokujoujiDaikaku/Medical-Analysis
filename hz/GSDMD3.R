# 清理环境
rm(list = ls())
gc()
getwd()
setwd(dir="C:/Users/hz/Desktop/Treg/33") #设置工作路径

# 加载包
required_packages <- c("GEOquery", "limma", "ggplot2", "dplyr", "tidyr", "ggpubr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("GEOquery", "limma")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# 检查下载的文件
cat("检查下载的文件...\n")
manual_files <- c(
  "GSE99254_NSCLC.TCell.S12346.TPM.txt.gz"
)

existing_files <- manual_files[file.exists(manual_files)]
cat("找到的文件:", existing_files, "\n")

# 重要修改1: 直接读取TPM文件而不是count文件
cat("\n=== 步骤1: 直接读取TPM矩阵 ===\n")
con <- gzfile("GSE99254_NSCLC.TCell.S12346.TPM.txt.gz", "rt")
# 先读取一行判断列名
first_line <- readLines(con, n = 1)
close(con)

# 根据第一行判断是否包含"symbol"列
has_symbol_col <- grepl("symbol", first_line, ignore.case = TRUE)

# 重新读取数据，正确处理基因名列
con <- gzfile("GSE99254_NSCLC.TCell.S12346.TPM.txt.gz", "rt")
if (has_symbol_col) {
  # 如果有symbol列，先读取数据再设置行名
  tpm_data <- read.delim(con, check.names = FALSE)
  
  # 检查并处理symbol列中的缺失值
  if ("symbol" %in% colnames(tpm_data)) {
    # 将空白字符串转换为NA
    tpm_data$symbol[tpm_data$symbol == ""] <- NA
    
    # 统计缺失值数量
    na_count <- sum(is.na(tpm_data$symbol))
    if (na_count > 0) {
      cat(sprintf("警告: symbol列中存在%d个缺失值，将这些行移除\n", na_count))
      tpm_data <- tpm_data[!is.na(tpm_data$symbol), ]
    }
    
    # 检查是否有重复的symbol
    if (any(duplicated(tpm_data$symbol))) {
      dup_count <- sum(duplicated(tpm_data$symbol))
      cat(sprintf("警告: 存在%d个重复的symbol，将添加后缀处理\n", dup_count))
      tpm_data$symbol <- make.unique(as.character(tpm_data$symbol))
    }
    
    # 将symbol列设为行名并删除该列
    rownames(tpm_data) <- tpm_data$symbol
    tpm_data <- tpm_data[, !colnames(tpm_data) %in% "symbol"]
    cat("已将'symbol'列设为基因名并移除该列\n")
  } else {
    stop("错误: 检测到有symbol列，但实际数据中未找到")
  }
} else {
  # 如果没有symbol列，假设第一列是基因名
  # 先读取数据查看第一列是否有缺失值
  temp_data <- read.delim(con, check.names = FALSE)
  first_col <- temp_data[, 1]
  
  # 处理第一列中的缺失值
  na_count <- sum(is.na(first_col) | first_col == "")
  if (na_count > 0) {
    cat(sprintf("警告: 第一列（基因名）中存在%d个缺失值，将这些行移除\n", na_count))
    temp_data <- temp_data[!is.na(first_col) & first_col != "", ]
  }
  
  # 检查重复基因名
  if (any(duplicated(first_col))) {
    dup_count <- sum(duplicated(first_col))
    cat(sprintf("警告: 存在%d个重复的基因名，将添加后缀处理\n", dup_count))
    temp_data[, 1] <- make.unique(as.character(temp_data[, 1]))
  }
  
  # 设置行名并删除第一列
  rownames(temp_data) <- temp_data[, 1]
  tpm_data <- temp_data[, -1, drop = FALSE]
  cat("已将第一列设为基因名\n")
}
close(con)
# 数据清理：将非数值转换为NA
cat("数据清理前 - TPM矩阵维度:", dim(tpm_data), "\n")
cat("检查数据类型...\n")

# 将所有列转换为数值型
for (col in colnames(tpm_data)) {
  tpm_data[[col]] <- as.numeric(as.character(tpm_data[[col]]))
}

# 移除包含NA的列（样本）
cols_before <- ncol(tpm_data)
tpm_data <- tpm_data[, colSums(is.na(tpm_data)) == 0]
cols_after <- ncol(tpm_data)
cat(sprintf("数据清理后 - 移除了 %d 个包含非数值的样本\n", cols_before - cols_after))

cat("清理后的TPM矩阵维度:", dim(tpm_data), "\n")
cat("TPM值范围:", range(as.matrix(tpm_data), na.rm = TRUE), "\n")
cat("前5行前5列数据:\n")
print(head(tpm_data[, 1:5]))
cat("前20个基因名:\n")
print(head(rownames(tpm_data), 20))
# 重要修改2: 改进的GSDMD基因查找方法
cat("\n=== 检查GSDMD基因表达 ===\n")

# 方法1: 尝试不同的基因名
gsdmd_aliases <- c("GSDMD", "GSMD", "DFNA5L", "C17orf59", "27106")
gsdmd_found <- FALSE
gsdmd_index <- NA

# 方法1: 精确匹配
for (alias in gsdmd_aliases) {
  if (alias %in% rownames(tpm_data)) {
    gsdmd_index <- which(rownames(tpm_data) == alias)
    gsdmd_found <- TRUE
    cat(sprintf("方法1: 找到精确匹配的基因 '%s'\n", alias))
    break
  }
}

# 方法2: 如果方法1失败，尝试模糊匹配
if (!gsdmd_found) {
  for (alias in gsdmd_aliases) {
    gsdmd_rows <- grep(alias, rownames(tpm_data), ignore.case = TRUE)
    if (length(gsdmd_rows) > 0) {
      gsdmd_index <- gsdmd_rows[1]
      gsdmd_found <- TRUE
      cat(sprintf("方法2: 找到模糊匹配的基因 '%s' -> '%s'\n", 
                  alias, rownames(tpm_data)[gsdmd_index]))
      break
    }
  }
}

# 方法3: 如果方法2失败，显示更多信息并让用户选择
if (!gsdmd_found) {
  cat("方法3: 显示包含'GSD'的基因供选择...\n")
  gsd_genes <- grep("GSD", rownames(tpm_data), ignore.case = TRUE)
  
  if (length(gsd_genes) > 0) {
    cat("找到以下包含'GSD'的基因:\n")
    for (i in seq_along(gsd_genes)) {
      cat(sprintf("%d: %s\n", i, rownames(tpm_data)[gsd_genes[i]]))
    }
    
    # 自动选择第一个匹配的基因
    selection <- 1
    if (selection >= 1 && selection <= length(gsd_genes)) {
      gsdmd_index <- gsd_genes[selection]
      gsdmd_found <- TRUE
      cat(sprintf("自动选择基因: %s\n", rownames(tpm_data)[gsdmd_index]))
    }
  }
}

# 如果找到GSDMD基因，提取其表达数据
if (gsdmd_found && !is.na(gsdmd_index)) {
  gsdmd_expression <- tpm_data[gsdmd_index, ]
  cat(sprintf("\n成功提取GSDMD基因表达数据\n"))
  cat(sprintf("基因名: %s\n", rownames(tpm_data)[gsdmd_index]))
  
  # 确保数据是数值型
  gsdmd_expression_numeric <- as.numeric(gsdmd_expression)
  
  cat(sprintf("表达值范围: %.3f - %.3f\n", min(gsdmd_expression_numeric), max(gsdmd_expression_numeric)))
  cat(sprintf("表达值中位数: %.3f\n", median(gsdmd_expression_numeric)))
  cat("前10个样本的表达值:\n")
  print(head(gsdmd_expression_numeric, 10))
} else {
  stop("错误: 未找到GSDMD基因，请检查基因名或数据文件")
}

# 重要修改3: 对GSDMD进行log2(TPM + 1)转换
cat("\n=== 步骤2: 对GSDMD进行log2(TPM + 1)转换 ===\n")
# 使用之前转换的数值型数据进行log2转换
gsdmd_expression_log <- log2(gsdmd_expression_numeric + 1)
cat(sprintf("转换后的表达值范围: %.3f - %.3f\n", min(gsdmd_expression_log), max(gsdmd_expression_log)))
cat(sprintf("转换后的表达值中位数: %.3f\n", median(gsdmd_expression_log)))

# 重要修改4: 严格基于PTH、TTH、PTR、TTR标记进行样本分组
cat("\n=== 步骤3: 基于PTH、TTH、PTR、TTR标记进行样本分组 ===\n")
sample_names <- colnames(tpm_data)
target_groups <- c("PTH", "TTH", "PTR", "TTR")

# 检查每个样本名中是否包含目标分组标记
sample_groups <- rep("Unclassified", length(sample_names))
for (group in target_groups) {
  sample_groups[grepl(group, sample_names, ignore.case = TRUE)] <- group
}

# 统计每个分组的样本数量
group_counts <- table(sample_groups)
cat("样本分组结果:\n")
print(group_counts)

# 检查是否有未分类的样本
unclassified_count <- sum(sample_groups == "Unclassified")
if (unclassified_count > 0) {
  cat(sprintf("\n注意: 有 %d 个样本未被分类\n", unclassified_count))
  cat("未分类的样本名:\n")
  print(sample_names[sample_groups == "Unclassified"])
}

# 创建绘图数据
plot_data <- data.frame(
  Sample = sample_names,
  GSDMD_Expression = gsdmd_expression_log,
  Group = sample_groups,
  stringsAsFactors = FALSE
)

# 过滤数据，只保留目标分组
plot_data_filtered <- plot_data[plot_data$Group %in% target_groups, ]

cat("\n用于分析的目标分组统计:\n")
final_counts <- table(plot_data_filtered$Group)
print(final_counts)

# 检查是否有足够的数据进行分析
if (nrow(plot_data_filtered) == 0) {
  stop("错误: 过滤后没有目标分组数据，请检查样本ID格式")
}

# 检查每种细胞类型是否都有数据
# 根据用户提供的映射关系：PTH(CD4+Tconv外周血)、TTH(CD4+Tconv肿瘤组织)、PTR(Treg外周血)、TTR(Treg肿瘤组织)
cat("\n根据用户提供的映射关系确定细胞类型和组织来源:\n")
cat("- PTH: CD4+Tconv细胞，PBMC组织\n")
cat("- TTH: CD4+Tconv细胞，肿瘤组织\n")
cat("- PTR: Treg细胞，PBMC组织\n")
cat("- TTR: Treg细胞，肿瘤组织\n")

# 重要修改4: 创建灵活的可视化图形
cat("\n=== 创建可视化图形 ===\n")

# 根据分组名称确定细胞类型和组织来源
plot_data_separated <- plot_data_filtered %>%
  mutate(
    CellType = case_when(
      Group == "PTH" ~ "CD4+Tconv",
      Group == "TTH" ~ "CD4+Tconv",
      Group == "PTR" ~ "Treg",
      Group == "TTR" ~ "Treg",
      TRUE ~ "Unknown"  # 添加默认情况，避免NA
    ),
    Tissue = case_when(
      Group == "PTH" ~ "PBMC",
      Group == "TTH" ~ "Tumor",
      Group == "PTR" ~ "PBMC",
      Group == "TTR" ~ "Tumor",
      TRUE ~ "Unknown"  # 添加默认情况，避免NA
    )
  ) %>%
  # 转换为因子并设置水平
  mutate(
    CellType = factor(CellType, levels = c("CD4+Tconv", "Treg", "Unknown")),
    Tissue = factor(Tissue, levels = c("PBMC", "Tumor", "Unknown"))
  )

# 显示分组转换结果
cat("\n分组转换结果预览:\n")
print(head(plot_data_separated[, c("Sample", "Group", "CellType", "Tissue", "GSDMD_Expression")]))

# 检查转换后的数据
cat("\n转换后的细胞类型分布:\n")
print(table(plot_data_separated$CellType))
cat("\n转换后的组织来源分布:\n")
print(table(plot_data_separated$Tissue))
cat("\n转换后的细胞类型和组织来源组合分布:\n")
print(table(plot_data_separated$CellType, plot_data_separated$Tissue))

# 过滤掉Unknown类型
plot_data_separated <- plot_data_separated %>%
  filter(CellType != "Unknown" & Tissue != "Unknown")

cat("\n过滤后的细胞类型和组织来源组合分布:\n")
print(table(plot_data_separated$CellType, plot_data_separated$Tissue))

# 设置颜色 - 与参考图一致的灰色调
plot_colors <- c("#CCCCCC", "#999999", "#666666", "#333333")

# 图形1：主要分析图（专注于同一种细胞类型在不同组织间的比较）
p1 <- ggplot(plot_data_separated, 
             aes(x = CellType, y = GSDMD_Expression, fill = Tissue)) +
  # 添加数据分布的点图
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), 
              alpha = 0.3, size = 1, color = "black") +
  # 添加条形图显示均值
  stat_summary(fun = mean, geom = "bar", position = position_dodge(0.9), alpha = 0.8) +
  # 添加误差线（标准误）
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(0.9), width = 0.2, color = "black") +
  # 设置颜色 - 使用对比明显的颜色
  scale_fill_manual(values = c("PBMC" = "#999999", "Tumor" = "#333333")) +
  # 设置坐标轴标签
  labs(
    title = "GSDMD Expression in Treg and CD4+Tconv Cells",
    x = "Cell Type",
    y = "GSDMD log2(TPM + 1)",
    fill = "Tissue Source"
  ) +
  # 设置主题
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title.y = element_text(angle = 90, vjust = 0.5, size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

# 添加样本数量标签
sample_counts <- plot_data_separated %>%
  group_by(CellType, Tissue) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    label = paste0("n = ", n)
  )

p1 <- p1 + geom_text(data = sample_counts, 
                     aes(x = CellType, y = -0.5, label = label, group = Tissue),
                     position = position_dodge(0.9),
                     vjust = 1, size = 3.5)

# 计算实际的统计显著性
cat("\n计算统计显著性...\n")
stat_results <- list()

# 分析所有可能的比较
possible_comparisons <- list(
  list(name = "Treg_PBMC_vs_Tumor", 
       group1 = filter(plot_data_separated, CellType == "Treg" & Tissue == "PBMC")$GSDMD_Expression,
       group2 = filter(plot_data_separated, CellType == "Treg" & Tissue == "Tumor")$GSDMD_Expression),
  list(name = "CD4Tconv_PBMC_vs_Tumor", 
       group1 = filter(plot_data_separated, CellType == "CD4+Tconv" & Tissue == "PBMC")$GSDMD_Expression,
       group2 = filter(plot_data_separated, CellType == "CD4+Tconv" & Tissue == "Tumor")$GSDMD_Expression)
)

# 执行可能的统计检验
for (comp in possible_comparisons) {
  if (length(comp$group1) >= 2 && length(comp$group2) >= 2) {
    test_result <- t.test(comp$group1, comp$group2)
    stat_results[[comp$name]] <- test_result$p.value
    cat(sprintf("%s: p = %.4e\n", comp$name, test_result$p.value))
  } else {
    stat_results[[comp$name]] <- NA
    cat(sprintf("%s: 样本量不足，无法进行t检验\n", comp$name))
    cat(sprintf("  - 组1样本数: %d\n", length(comp$group1)))
    cat(sprintf("  - 组2样本数: %d\n", length(comp$group2)))
  }
}

# 添加可能的统计显著性标注
max_y <- max(plot_data_separated$GSDMD_Expression, na.rm = TRUE)
y_position <- max_y * 1.05
y_increment <- max_y * 0.1

# 为每个有效的比较添加标注
annotation_count <- 0
for (comp in possible_comparisons) {
  if (!is.na(stat_results[[comp$name]])) {
    if (comp$name == "Treg_PBMC_vs_Tumor") {
      # 在分组条形图中，Treg细胞的两个组织样本本位于x轴位置1的两侧
      p1 <- p1 + 
        annotate("segment", x = 1 - 0.225, xend = 1 + 0.225, 
                 y = y_position + annotation_count * y_increment, 
                 yend = y_position + annotation_count * y_increment, color = "black") +
        annotate("text", x = 1, 
                 y = y_position + (annotation_count + 0.1) * y_increment,
                 label = sprintf("P = %.2e", stat_results[[comp$name]]), size = 3.5)
      annotation_count <- annotation_count + 1
    } else if (comp$name == "CD4Tconv_PBMC_vs_Tumor") {
      # 在分组条形图中，CD4+Tconv细胞的两个组织样本位于x轴位置2的两侧
      p1 <- p1 + 
        annotate("segment", x = 2 - 0.225, xend = 2 + 0.225, 
                 y = y_position + annotation_count * y_increment, 
                 yend = y_position + annotation_count * y_increment, color = "black") +
        annotate("text", x = 2, 
                 y = y_position + (annotation_count + 0.1) * y_increment,
                 label = sprintf("P = %.2e", stat_results[[comp$name]]), size = 3.5)
      annotation_count <- annotation_count + 1
    }
  }
}

print(p1)

# 保存结果
cat("\n=== 保存结果 ===\n")
backup_dir <- "GSDMD_analysis_pth_tth_ptr_ttr"
if (!dir.exists(backup_dir)) {
  dir.create(backup_dir)
}

tryCatch({
  ggsave(file.path(backup_dir, "GSDMD_Treg_CD4Tconv_expression_final.png"), 
         p1, width = 8, height = 6, dpi = 300)
  write.csv(plot_data_filtered, file.path(backup_dir, "GSDMD_expression_data.csv"), row.names = FALSE)
  
  # 保存未分类的样本信息
  unclassified_data <- plot_data[plot_data$Group == "Unclassified", ]
  if (nrow(unclassified_data) > 0) {
    write.csv(unclassified_data, file.path(backup_dir, "Unclassified_samples.csv"), row.names = FALSE)
  }
  
  cat("结果已保存到:", backup_dir, "目录\n")
}, error = function(e) {
  cat("保存失败:", e$message, "\n")
})

# 统计分析
cat("\n=== 统计分析 ===\n")

# 手动计算统计量
calculate_stats <- function(data) {
  if (length(data) == 0) return(list(mean = NA, sd = NA, sem = NA, n = 0))
  
  valid_data <- data[is.finite(data)]
  n <- length(valid_data)
  
  if (n == 0) return(list(mean = NA, sd = NA, sem = NA, n = 0))
  if (n == 1) return(list(mean = mean(valid_data), sd = NA, sem = NA, n = n))
  
  return(list(
    mean = mean(valid_data),
    sd = sd(valid_data),
    sem = sd(valid_data) / sqrt(n),
    n = n
  ))
}

# 分组统计
cell_types <- unique(plot_data_separated$CellType)
tissues <- unique(plot_data_separated$Tissue)

cat("描述性统计 (log2(TPM+1)):\n")
for (ct in cell_types) {
  for (tissue in tissues) {
    subset_data <- plot_data_separated %>% 
      filter(CellType == ct & Tissue == tissue) %>% 
      pull(GSDMD_Expression)
    
    stats <- calculate_stats(subset_data)
    
    cat(sprintf("- %s在%s中: ", ct, tissue))
    if (stats$n > 0) {
      cat(sprintf("均值=%.3f", stats$mean))
      if (!is.na(stats$sd)) cat(sprintf(", SD=%.3f", stats$sd))
      cat(sprintf(", N=%d\n", stats$n))
    } else {
      cat("无数据\n")
    }
  }
}

# T检验比较
cat("\nT检验结果:\n")
for (ct in cell_types) {
  tumor_data <- plot_data_separated %>% 
    filter(CellType == ct & Tissue == "Tumor") %>% 
    pull(GSDMD_Expression)
  
  pbmc_data <- plot_data_separated %>% 
    filter(CellType == ct & Tissue == "PBMC") %>% 
    pull(GSDMD_Expression)
  
  if (length(tumor_data) >= 2 && length(pbmc_data) >= 2) {
    t_test <- t.test(tumor_data, pbmc_data)
    cat(sprintf("%s细胞 - 肿瘤组织 vs PBMC: p = %.4e\n", ct, t_test$p.value))
    
    # 计算均值差异百分比
    mean_tumor <- mean(tumor_data)
    mean_pbmc <- mean(pbmc_data)
    diff_percent <- ((mean_tumor - mean_pbmc) / mean_pbmc) * 100
    
    if (diff_percent > 0) {
      cat(sprintf("  - 肿瘤组织中的表达比PBMC高 %.1f%%\n", diff_percent))
    } else if (diff_percent < 0) {
      cat(sprintf("  - 肿瘤组织中的表达比PBMC低 %.1f%%\n", abs(diff_percent)))
    } else {
      cat("  - 肿瘤组织和PBMC中的表达水平相同\n")
    }
  } else {
    cat(sprintf("%s细胞 - 样本量不足，跳过t检验\n", ct))
    cat(sprintf("  - 肿瘤组织样本数: %d\n", length(tumor_data)))
    cat(sprintf("  - PBMC样本数: %d\n", length(pbmc_data)))
  }
}

# 在控制台输出完整结果
cat("\n=== 完整分析结果 ===\n")
cat("数据处理流程:\n")
cat("1. 读取TPM数据并进行数据清理\n")
cat("2. 使用多种方法查找GSDMD基因\n")
cat("3. 对GSDMD进行log2(TPM + 1)转换\n")
cat("4. 仅基于PTH、TTH、PTR、TTR标记进行样本分组\n\n")

cat("关键发现:\n")
cat("- 成功清理了TPM数据，移除了非数值和NA值\n")
cat("- 使用多种方法成功找到了GSDMD基因\n")
cat("- 严格基于用户提供的标记（PTH、TTH、PTR、TTR）进行分组\n")
cat("- 可视化结果根据实际可用数据动态调整\n\n")

cat("样本信息:\n")
for (group in target_groups) {
  count <- ifelse(group %in% names(final_counts), final_counts[group], 0)
  cat(sprintf("- %s: %d 个样本\n", group, count))
}
cat(sprintf("- 未分类样本: %d 个\n", unclassified_count))

if (unclassified_count > 0) {
  cat("\n注意: 有未分类的样本，可能是因为样本ID中不包含PTH、TTH、PTR、TTR标记\n")
  cat("建议检查样本ID的命名规则，确保包含正确的细胞类型标记\n")
}

cat("\n=== 分析完成 ===\n")
