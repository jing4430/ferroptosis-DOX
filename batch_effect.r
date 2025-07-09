#使用sva
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("sva")

library(sva)

# 读取FPKM数据
fpkm_data <- read.table("fpkm_data.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# 查看数据
head(fpkm_data)

# 读取样本信息
sample_info <- read.csv("sample_info.csv", header = TRUE)

# 查看样本信息
head(sample_info)

# 确保样本顺序一致
if (!all(colnames(fpkm_data) == sample_info$sample)) {
  stop("样本顺序不一致！请检查数据。")
}

# 提取批次和分组信息
batch <- sample_info$batch
group <- sample_info$group

# 将FPKM数据转换为矩阵
fpkm_matrix <- as.matrix(fpkm_data)

# 使用ComBat去除批次效应
adjusted_data <- ComBat(dat = fpkm_matrix, batch = batch, mod = model.matrix(~group))

# 将调整后的数据转换回数据框
adjusted_data <- as.data.frame(adjusted_data)

#数据平移
adjusted_data <- adjusted_data - min(adjusted_data)

# 保存调整后的数据
write.table(adjusted_data, file = "adjusted_fpkm.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# 安装并加载ggplot2和ggfortify
if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2")
if (!requireNamespace("ggfortify", quietly = TRUE))
    install.packages("ggfortify")

library(ggplot2)
library(ggfortify)

# 检查是否存在常数列或全为零列
zero_variance_genes <- apply(adjusted_data, 1, function(x) var(x) == 0)

# 打印常数列或全为零列的基因
print(names(zero_variance_genes[zero_variance_genes]))

# 去除常数列或全为零列
adjusted_data <- adjusted_data[!zero_variance_genes, ]

# 进行PCA分析
pca_result <- prcomp(t(adjusted_data), scale. = TRUE)

# 可视化PCA结果
autoplot(pca_result, data = sample_info, colour = 'group', shape = 'batch')





#使用 Harmony

# 安装并加载必要的 R 包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("harmony")

library(harmony)
library(ggplot2)
library(ggfortify)

# 读取 FPKM 数据
fpkm_data <- read.table("fpkm_data.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# 读取样本信息
sample_info <- read.csv("batch_info.csv", header = TRUE)

# 对 FPKM 数据进行对数变换
fpkm_matrix <- as.matrix(fpkm_data)
fpkm_matrix <- log2(fpkm_matrix + 1)

# 确保样本顺序一致
if (!all(colnames(fpkm_matrix) == sample_info$sample)) {
  stop("样本顺序不一致！请检查数据。")
}

# 提取分组和批次信息
group <- as.factor(sample_info$group)
batch <- as.factor(sample_info$batch)

# 检查批次与分组的分布
table(batch, group)

# 使用 Harmony 去除批次效应
harmony_data <- HarmonyMatrix(
  data_mat = fpkm_matrix,
  meta_data = sample_info,
  vars_use = 'batch',
  theta = 2,  # 增加 theta 值
  lambda = 1, # 调整 lambda 值
  do_pca = TRUE
)

# 检查 Harmony 输出的样本数
if (ncol(harmony_data) != ncol(fpkm_matrix)) {
  stop("Harmony 输出的样本数与输入数据不一致！")
}

#数据平移
harmony_data <- harmony_data - min(harmony_data)

# 将 Harmony 输出转换为数据框
harmony_data <- as.data.frame(t(harmony_data))  # 转置矩阵，使行是样本，列是主成分
harmony_data$sample <- rownames(harmony_data)  # 添加样本名

# 合并样本信息
harmony_data <- merge(harmony_data, sample_info, by = "sample")

# 检查并去除常数列或全为零列
zero_variance_cols <- apply(harmony_data[, 1:(ncol(harmony_data) - 3)], 2, function(x) var(x) == 0)
harmony_data <- harmony_data[, !zero_variance_cols]

# 检查数据是否为空
if (ncol(harmony_data) == 0) {
  stop("去除常数列或全为零列后，数据为空！")
}

# 进行 PCA 分析
pca_result <- prcomp(harmony_data[, 1:(ncol(harmony_data) - 3)], scale. = TRUE)  # 排除样本名、分组和批次列

# 可视化 PCA 结果
autoplot(pca_result, data = harmony_data, colour = 'group', shape = 'batch')

# 保存调整后的数据
write.table(harmony_data, file = "harmony_adjusted_fpkm.txt", sep = "\t", quote = FALSE, row.names = FALSE)





#使用limma

# 检查FPKM矩阵（示例）
fpkm_matrix <- read.table("fpkm_data.txt", header = TRUE, row.names = 1)
summary(fpkm_matrix)  # 确认无非负值或异常值

fpkm_log2 <- log2(fpkm_matrix + 1)  # 避免log(0)

library(limma)

# 读取样本信息
sample_info <- read.csv("batch_info.csv", header = TRUE)

# 定义批次变量和主变量（实验分组）
batch <- factor(sample_info$batch)  # 批次信息（需为因子）
group <- factor(sample_info$group)  # 实验组别（如control/DOX）

# 创建设计矩阵（主变量为实验分组）
design <- model.matrix(~ group)

# 去除批次效应
corrected_data <- removeBatchEffect(
  x = fpkm_log2,
  batch = batch,
  design = design
)

#数据平移
corrected_data <- corrected_data - min(corrected_data)