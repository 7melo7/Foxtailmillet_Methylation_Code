library("tidyverse")
library("ggplot2")

setwd("D:\\Projects\\methy_revison_round_2\\Sanger_PCR")

judge_methylation <- function(df, 
                              unmethyl_col = "Number_of_unmethylated_C",
                              total_col = "Sample_number",
                              threshold = 0.007,
                              alpha = 0.05,
                              alternative = "greater") {
  
  # ---- 1. 参数检查 ----
  if (!is.data.frame(df)) stop("df 必须是数据框")
  if (!unmethyl_col %in% colnames(df)) stop(paste("列", unmethyl_col, "不存在"))
  if (!total_col %in% colnames(df)) stop(paste("列", total_col, "不存在"))
  if (threshold <= 0 || threshold >= 1) stop("threshold 必须在 (0,1) 之间")
  if (alpha <= 0 || alpha >= 1) stop("alpha 必须在 (0,1) 之间")
  if (!alternative %in% c("greater", "less", "two.sided")) {
    stop("alternative 必须是 'greater', 'less' 或 'two.sided'")
  }
  
  # ---- 2. 提取数据 ----
  unmethyl <- df[[unmethyl_col]]
  total <- df[[total_col]]
  
  # 确保数值型
  if (!is.numeric(unmethyl)) stop(paste(unmethyl_col, "必须为数值型"))
  if (!is.numeric(total)) stop(paste(total_col, "必须为数值型"))
  
  # ---- 3. 计算甲基化Clone数 ----
  methylated <- total - unmethyl
  
  # 检查是否有负数（异常值）
  if (any(methylated < 0, na.rm = TRUE)) {
    warning("存在甲基化Clone数为负值，请检查数据（未甲基化数 > 总数）")
  }
  
  # ---- 4. 将计算结果加入数据框 ----
  df$Methylated_C <- methylated
  
  # ---- 5. 初始化状态列 ----
  df$mC_status <- NA_character_
  
  # ---- 6. 循环进行二项检验 ----
  for (i in seq_len(nrow(df))) {
    n <- total[i]
    x <- methylated[i]
    
    # 处理总数为0或缺失值的情况
    if (is.na(n) || n == 0 || is.na(x)) {
      df$mC_status[i] <- NA_character_
      next
    }
    
    # 处理 x 为 NA 的情况
    if (is.na(x)) {
      df$mC_status[i] <- NA_character_
      next
    }
    
    # 二项检验（使用 pbinom 快速计算）
    if (alternative == "greater") {
      # P(X >= x) = 1 - P(X <= x-1)
      p_value <- pbinom(x - 1, size = n, prob = threshold, lower.tail = FALSE)
    } else if (alternative == "less") {
      # P(X <= x)
      p_value <- pbinom(x, size = n, prob = threshold, lower.tail = TRUE)
    } else { # two.sided
      # 双侧检验（近似）
      p_lower <- pbinom(x, size = n, prob = threshold)
      p_upper <- pbinom(x - 1, size = n, prob = threshold, lower.tail = FALSE)
      p_value <- 2 * min(p_lower, p_upper)
      p_value <- min(p_value, 1)  # 截断到1
    }
    
    # 根据 p 值判断
    if (!is.na(p_value) && p_value < alpha) {
      df$mC_status[i] <- "mC"
    } else {
      df$mC_status[i] <- "C"
    }
  }
  
  # ---- 7. 返回结果 ----
  return(df)
}

get_methylation_summary_with_ci <- function(df,
                                            acc_col = "Acc",
                                            context_col = "site.type",
                                            status_col = "mC_status",
                                            conf_level = 0.95) {
  
  library(dplyr)
  
  df_valid <- df %>% filter(!is.na(.data[[status_col]]))
  
  if (nrow(df_valid) == 0) {
    warning("没有有效的 mC_status 数据")
    return(data.frame())
  }
  
  # 按 Acc + site.type 统计
  result_by_group <- df_valid %>%
    group_by(.data[[acc_col]], .data[[context_col]]) %>%
    summarise(
      methylated_count = sum(.data[[status_col]] == "mC", na.rm = TRUE),
      total_count = n(),
      methylation_rate = methylated_count / total_count,
      # 计算二项分布置信区间（Wilson方法）
      ci_lower = ifelse(total_count > 0,
                        qbeta((1 - conf_level) / 2, methylated_count + 1, total_count - methylated_count + 1),
                        NA),
      ci_upper = ifelse(total_count > 0,
                        qbeta((1 + conf_level) / 2, methylated_count + 1, total_count - methylated_count + 1),
                        NA),
      .groups = "drop"
    )
  
  # 按 Acc 汇总
  result_overall <- df_valid %>%
    group_by(.data[[acc_col]]) %>%
    summarise(
      methylated_count = sum(.data[[status_col]] == "mC", na.rm = TRUE),
      total_count = n(),
      methylation_rate = methylated_count / total_count,
      ci_lower = qbeta((1 - conf_level) / 2, methylated_count + 1, total_count - methylated_count + 1),
      ci_upper = qbeta((1 + conf_level) / 2, methylated_count + 1, total_count - methylated_count + 1),
      .groups = "drop"
    ) %>%
    mutate(!!context_col := "Overall")
  
  result_final <- bind_rows(result_by_group, result_overall) %>%
    arrange(.data[[acc_col]], ifelse(.data[[context_col]] == "Overall", 1, 0))
  
  return(result_final)
}

Loci1_Acc.Lst <- c("C13", "L29", "L32", "Q14", "Q18", "Q24")

Loci1.Dat <- lapply(Loci1_Acc.Lst, function(acc) {
  file_path <- paste0("Share_NonTE_Cultivar_Loci2/Hyper_Loci2_", acc, ".csv")
  
  if (file.exists(file_path)) {
    temp_data <- read.csv(file_path, stringsAsFactors = FALSE)
    temp_data$Acc <- acc  # 添加品种名列
    return(temp_data)
  } else {
    warning(paste("文件不存在:", file_path))
    return(NULL)  # 如果文件不存在，返回NULL
  }
})

Loci1.Dat <- Loci1.Dat[!sapply(Loci1.Dat, is.null)]
Loci1.Dat <- do.call(rbind, Loci1.Dat)

Loci1.Dat$Methylated_C <- Loci1.Dat$Sample_number - Loci1.Dat$Number_of_unmethylated_C
Loci1.Dat <- judge_methylation(Loci1.Dat)
Loci1.Summary <- get_methylation_summary_with_ci(Loci1.Dat)

Loci2_Acc.Lst <- c("L6", "L32", "L36", "Q24", "Q28", "Q37")

Loci2.Dat <- lapply(Loci2_Acc.Lst, function(acc) {
  file_path <- paste0("Share_NonTE_Wild_Loci1/Hypo_Loci1_", acc, ".csv")
  
  if (file.exists(file_path)) {
    temp_data <- read.csv(file_path, stringsAsFactors = FALSE)
    temp_data$Acc <- acc  # 添加品种名列
    return(temp_data)
  } else {
    warning(paste("文件不存在:", file_path))
    return(NULL)  # 如果文件不存在，返回NULL
  }
})

Loci2.Dat <- Loci2.Dat[!sapply(Loci2.Dat, is.null)]
Loci2.Dat <- do.call(rbind, Loci2.Dat)

Loci2.Dat$Methylated_C <- Loci2.Dat$Sample_number - Loci2.Dat$Number_of_unmethylated_C
Loci2.Dat <- judge_methylation(Loci2.Dat)
Loci2.Summary <- get_methylation_summary_with_ci(Loci2.Dat)


Loci3_Acc.Lst <- c("L6", "L29", "L32", "Q18", "Q24", "Q28")

Loci3.Dat <- lapply(Loci3_Acc.Lst, function(acc) {
  file_path <- paste0("Share_NonTE_Wild_Loci2/Hypo_Loci2_", acc, ".csv")
  
  if (file.exists(file_path)) {
    temp_data <- read.csv(file_path, stringsAsFactors = FALSE)
    temp_data$Acc <- acc  # 添加品种名列
    return(temp_data)
  } else {
    warning(paste("文件不存在:", file_path))
    return(NULL)  # 如果文件不存在，返回NULL
  }
})

Loci3.Dat <- Loci3.Dat[!sapply(Loci3.Dat, is.null)]
Loci3.Dat <- do.call(rbind, Loci3.Dat)

Loci3.Dat$Methylated_C <- Loci3.Dat$Sample_number - Loci3.Dat$Number_of_unmethylated_C
Loci3.Dat <- judge_methylation(Loci3.Dat)
Loci3.Summary <- get_methylation_summary_with_ci(Loci3.Dat)


Loci4_Acc.Lst <- c("L6", "L29", "L32", "Q18", "Q24", "Q28")

Loci4.Dat <- lapply(Loci4_Acc.Lst, function(acc) {
  file_path <- paste0("Share_NonTE_Wild_Loci3/Hypo_Loci3_", acc, ".csv")
  
  if (file.exists(file_path)) {
    temp_data <- read.csv(file_path, stringsAsFactors = FALSE)
    temp_data$Acc <- acc  # 添加品种名列
    return(temp_data)
  } else {
    warning(paste("文件不存在:", file_path))
    return(NULL)  # 如果文件不存在，返回NULL
  }
})

Loci4.Dat <- Loci4.Dat[!sapply(Loci4.Dat, is.null)]
Loci4.Dat <- do.call(rbind, Loci4.Dat)

Loci4.Dat$Methylated_C <- Loci4.Dat$Sample_number - Loci4.Dat$Number_of_unmethylated_C
Loci4.Dat <- judge_methylation(Loci4.Dat)
Loci4.Summary <- get_methylation_summary_with_ci(Loci4.Dat)



Loci5_Acc.Lst <- c("C13", "L29", "L32", "Q14", "Q24", "Q28")

Loci5.Dat <- lapply(Loci5_Acc.Lst, function(acc) {
  file_path <- paste0("Share_TE_Cultivar_Loci1/Hyper_Loci1_", acc, ".csv")
  
  if (file.exists(file_path)) {
    temp_data <- read.csv(file_path, stringsAsFactors = FALSE)
    temp_data$Acc <- acc  # 添加品种名列
    return(temp_data)
  } else {
    warning(paste("文件不存在:", file_path))
    return(NULL)  # 如果文件不存在，返回NULL
  }
})

Loci5.Dat <- Loci5.Dat[!sapply(Loci5.Dat, is.null)]
Loci5.Dat <- do.call(rbind, Loci5.Dat)

Loci5.Dat$Methylated_C <- Loci5.Dat$Sample_number - Loci5.Dat$Number_of_unmethylated_C
Loci5.Dat <- judge_methylation(Loci5.Dat)
Loci5.Summary <- get_methylation_summary_with_ci(Loci5.Dat)



Loci6_Acc.Lst <- c("C13", "L29", "L32", "Q14", "Q24", "Q28")

Loci6.Dat <- lapply(Loci6_Acc.Lst, function(acc) {
  file_path <- paste0("Share_TE_Cultivar_Loci2/Hyper_Loci2_", acc, ".csv")
  
  if (file.exists(file_path)) {
    temp_data <- read.csv(file_path, stringsAsFactors = FALSE)
    temp_data$Acc <- acc  # 添加品种名列
    return(temp_data)
  } else {
    warning(paste("文件不存在:", file_path))
    return(NULL)  # 如果文件不存在，返回NULL
  }
})

Loci6.Dat <- Loci6.Dat[!sapply(Loci6.Dat, is.null)]
Loci6.Dat <- do.call(rbind, Loci6.Dat)
Loci6.Dat$Methylated_C <- Loci6.Dat$Sample_number - Loci6.Dat$Number_of_unmethylated_C
Loci6.Dat <- judge_methylation(Loci6.Dat)
Loci6.Summary <- get_methylation_summary_with_ci(Loci6.Dat)


Loci7_Acc.Lst <- c("L6", "L29", "L36", "Q14", "Q18", "Q24")

Loci7.Dat <- lapply(Loci7_Acc.Lst, function(acc) {
  file_path <- paste0("Share_TE_Cultivar_Loci3/Hyper_Loci3_", acc, ".csv")
  
  if (file.exists(file_path)) {
    temp_data <- read.csv(file_path, stringsAsFactors = FALSE)
    temp_data$Acc <- acc  # 添加品种名列
    return(temp_data)
  } else {
    warning(paste("文件不存在:", file_path))
    return(NULL)  # 如果文件不存在，返回NULL
  }
})

Loci7.Dat <- Loci7.Dat[!sapply(Loci7.Dat, is.null)]
Loci7.Dat <- do.call(rbind, Loci7.Dat)

Loci7.Dat$Methylated_C <- Loci7.Dat$Sample_number - Loci7.Dat$Number_of_unmethylated_C
Loci7.Dat <- judge_methylation(Loci7.Dat)
Loci7.Summary <- get_methylation_summary_with_ci(Loci7.Dat)

Loci8_Acc.Lst <- c("L29", "L32", "L36", "Q18", "Q28", "Q37")

Loci8.Dat <- lapply(Loci8_Acc.Lst, function(acc) {
  file_path <- paste0("Share_TE_Wild_Loci1/Hypo_Loci1_", acc, ".csv")
  
  if (file.exists(file_path)) {
    temp_data <- read.csv(file_path, stringsAsFactors = FALSE)
    temp_data$Acc <- acc  # 添加品种名列
    return(temp_data)
  } else {
    warning(paste("文件不存在:", file_path))
    return(NULL)  # 如果文件不存在，返回NULL
  }
})

Loci8.Dat <- Loci8.Dat[!sapply(Loci8.Dat, is.null)]
Loci8.Dat <- do.call(rbind, Loci8.Dat)

Loci8.Dat$Methylated_C <- Loci8.Dat$Sample_number - Loci8.Dat$Number_of_unmethylated_C
Loci8.Dat <- judge_methylation(Loci8.Dat)
Loci8.Summary <- get_methylation_summary_with_ci(Loci8.Dat)

Loci9_Acc.Lst <- c("L29", "L32", "L36", "Q14", "Q18", "Q24")

Loci9.Dat <- lapply(Loci9_Acc.Lst, function(acc) {
  file_path <- paste0("Share_TE_Wild_Loci2/Hypo_Loci2_", acc, ".csv")
  
  if (file.exists(file_path)) {
    temp_data <- read.csv(file_path, stringsAsFactors = FALSE)
    temp_data$Acc <- acc  # 添加品种名列
    return(temp_data)
  } else {
    warning(paste("文件不存在:", file_path))
    return(NULL)  # 如果文件不存在，返回NULL
  }
})

Loci9.Dat <- Loci9.Dat[!sapply(Loci9.Dat, is.null)]
Loci9.Dat <- do.call(rbind, Loci9.Dat)

Loci9.Dat$Methylated_C <- Loci9.Dat$Sample_number - Loci9.Dat$Number_of_unmethylated_C
Loci9.Dat <- judge_methylation(Loci9.Dat)
Loci9.Summary <- get_methylation_summary_with_ci(Loci9.Dat)



Loci10_Acc.Lst <- c("L6", "L29", "L32", "Q24", "Q28", "Q37")

Loci10.Dat <- lapply(Loci10_Acc.Lst, function(acc) {
  file_path <- paste0("Share_TE_Wild_Loci3/Hypo_Loci3_", acc, ".csv")
  
  if (file.exists(file_path)) {
    temp_data <- read.csv(file_path, stringsAsFactors = FALSE)
    temp_data$Acc <- acc  # 添加品种名列
    return(temp_data)
  } else {
    warning(paste("文件不存在:", file_path))
    return(NULL)  # 如果文件不存在，返回NULL
  }
})

Loci10.Dat <- Loci10.Dat[!sapply(Loci10.Dat, is.null)]
Loci10.Dat <- do.call(rbind, Loci10.Dat)

Loci10.Dat$Methylated_C <- Loci10.Dat$Sample_number - Loci10.Dat$Number_of_unmethylated_C
Loci10.Dat <- judge_methylation(Loci10.Dat)
Loci10.Summary <- get_methylation_summary_with_ci(Loci10.Dat)


Dat <- rbind(Loci1.Summary, Loci2.Summary, Loci3.Summary,
             Loci4.Summary, Loci5.Summary, Loci6.Summary,
             Loci7.Summary, Loci8.Summary, Loci9.Summary, Loci10.Summary)

Dat$Pop <- sapply(Dat$Acc, function(x){
  if(str_starts(x, "Q")){return("Wild")}
  return("Cultivar")
})
Dat$Pop <- factor(Dat$Pop, levels = c("Wild", "Cultivar"))

ggplot(filter(Dat, site.type=="Overall"), 
       aes(x = Pop, y = 100*methylation_rate, fill = Pop)) +
  stat_boxplot(geom="errorbar", width=0.3, size=0.5) + 
  geom_boxplot(width = 0.6, size = 0.5, outlier.shape = NA) +   # 箱线
  geom_point(position = position_jitter(width = 0.3, height = 0.1),  # 散点
             alpha = 1, size = 0.3) +
  scale_fill_manual(values = c("#62baa1", "#2e6b9a")) +
  coord_cartesian(ylim = c(0, 100)) +
  #scale_y_continuous(limits = c(0,100)) +
  labs(x = "Group", y = "Value") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 8)) +
  facet_grid()

ggsave("Fig/Share.mC_C.pdf", width = 2.7, height = 2.2)


ggplot(filter(Dat, site.type!="Overall"), 
       aes(x = Pop, y = 100*methylation_rate, fill = Pop)) +
  stat_boxplot(geom="errorbar", width=0.3, size=0.5) + 
  geom_boxplot(width = 0.6, size = 0.5, outlier.shape = NA) +   # 箱线
  geom_point(position = position_jitter(width = 0.3, height = 0.1),  # 散点
             alpha = 1, size = 0.3) +
  scale_fill_manual(values = c("#62baa1", "#2e6b9a")) +
  coord_cartesian(ylim = c(0, 100)) +
  #scale_y_continuous(limits = c(0,100)) +
  labs(x = "Group", y = "Value") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 8)) +
  facet_grid(.~site.type)

ggsave("Fig/Share.mC_C.Context.pdf", width = 4, height = 2)

write_tsv(filter(Dat, site.type=="Overall"), "Fig1h.Shared.txt")

wilcox.test(filter(Dat, site.type=="Overall", Pop=="Cultivar", site.type=="Overall")$methylation_rate,
            filter(Dat, site.type=="Overall", Pop=="Wild", site.type=="Overall")$methylation_rate)



