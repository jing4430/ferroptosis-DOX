# 安装必要包（若未安装）
# install.packages("limma")
# install.packages("dplyr")
# install.packages("ggplot2")

# 加载包
library(limma)
library(dplyr)
library(ggplot2)

# 读取数据（假设文件以制表符分隔，第一列为基因名）
fpkm_data <- read.table("limma_CUT_ZERO.txt", header = TRUE, sep = "\t", row.names = 1)

# 将行名（基因名）转换为列，便于合并
fpkm_processed <- fpkm_data %>%
  tibble::rownames_to_column("Gene_Symbol") %>%
  group_by(Gene_Symbol) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%  # 合并重复基因，取均值
  tibble::column_to_rownames("Gene_Symbol")

# 查看处理后的维度（确保重复基因已合并）
dim(fpkm_processed)

# 定义分组（假设前3列为对照组，后3列为处理组）
group <- factor(c(rep("control", 13), rep("DOX", 13)))

# 过滤：保留在至少2个样本中 FPKM ≥ 1 的基因
keep_genes <- rowSums(fpkm_processed > 0) >= 7
fpkm_filtered <- fpkm_processed[keep_genes, ]

# 加1伪计数后取 log2
#fpkm_log2 <- log2(fpkm_filtered + 1)

control_mean <- rowMeans(fpkm_filtered[, 1:13])
DOX_mean <- rowMeans(fpkm_filtered[, 14:26])
log2FC <- log2((DOX_mean + 1e-3) / (control_mean + 1e-3))  # 加伪计数避免除零

# 创建设计矩阵
design <- model.matrix(~ group)

# 拟合线性模型
fit <- lmFit(fpkm_filtered, design)
fit <- eBayes(fit)

# 提取结果
results <- topTable(fit, coef = "groupDOX", number = Inf, sort.by = "none")
p_values <- results$P.Value
q_values <- p.adjust(p_values, method = "fdr")

result_table <- data.frame(
  Gene_Symbol = rownames(fpkm_filtered),
  fpkm_filtered,
  log2FC = log2FC,
  p_value = p_values,
  q_value = q_values,
  stringsAsFactors = FALSE
) %>%
  arrange(Gene_Symbol)  # 按GENE_SYMBOL排序

# 保存结果
write.csv(result_table, "merged_and_recalculated_results.csv", row.names = FALSE)

