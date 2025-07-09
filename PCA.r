# 加载必要包
library(data.table)
library(ggplot2)
library(ggrepel)

# 读取已校正的FPKM数据（假设列为样本，行为基因）
fpkm_data <- fread("batch_corrected_data.txt", sep = "\t") %>% 
  as.data.frame()
rownames(fpkm_data) <- fpkm_data[,1]
fpkm_data <- fpkm_data[,-1]

# 创建批次元数据
batch_info <- data.frame(
  Sample = colnames(fpkm_data),
  Batch = factor(rep(paste0("Batch",1:3), times = c(3,4,3))),
  Group = factor(rep(c("Control","Treatment"), times = c(5,5)))  # 假设实验组信息
)

# 对数转换（适用于校正后的FPKM值）
log_fpkm <- log2(fpkm_data + 0.1)  # 添加伪计数处理负值

# 转置数据（样本×基因）
pca_input <- t(log_fpkm)

# 执行PCA分析（无需再缩放）
pca_result <- prcomp(pca_input, center = TRUE, scale. = FALSE)

# 提取主成分得分并合并分组信息
pca_scores <- data.frame(pca_result$x[, 1:2])  # 取前两个主成分
pca_scores$sample <- rownames(pca_scores)
pca_scores <- left_join(pca_scores, sample_info, by = "sample")

# ----------------------------- 3. 绘制PCA图 -----------------------------
ggplot(pca_scores, aes(x = PC1, y = PC2, color = group, shape = batch)) +
  geom_point(size = 4, alpha = 0.8) +  # 调整点大小和透明度
  labs(title = "PCA of FPKM Expression",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")) +
  scale_color_manual(values = c("control" = "#1f77b4", "DOX" = "#ff7f0e")) +  # 自定义颜色
  scale_shape_manual(values = c(16, 17, 15, 3, 8)) +  # 定义5种批次形状
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right")

# 保存高清图（可选）
ggsave("PCA_plot.pdf", width = 8, height = 6, dpi = 300)



# 加载必要包
library(data.table)
library(ggplot2)
library(ggrepel)

# 读取已校正的FPKM数据（假设列为样本，行为基因）
fpkm_data <- fread("batch_corrected_data.txt", sep = "\t") %>% 
  as.data.frame()
rownames(fpkm_data) <- fpkm_data[,1]
fpkm_data <- fpkm_data[,-1]

# 创建批次元数据
batch_info <- data.frame(
  Sample = colnames(fpkm_data),
  Batch = factor(rep(paste0("Batch",1:3), times = c(3,4,3))),
  Group = factor(rep(c("Control","Treatment"), times = c(5,5)))  # 假设实验组信息
)

# 对数转换（适用于校正后的FPKM值）
log_fpkm <- log2(fpkm_data + 0.1)  # 添加伪计数处理负值

# 转置数据（样本×基因）
pca_input <- t(log_fpkm)

# 执行PCA分析（无需再缩放）
pca_result <- prcomp(pca_input, center = TRUE, scale. = FALSE)

# 计算方差贡献
var_exp <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 1)

ggplot_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Batch = sample_info$Batch,
  Group = sample_info$Group
)

ggplot(ggplot_df, aes(PC1, PC2)) +
  geom_point(aes(color=Group, shape=Batch), size=5, alpha=0.8) +
  stat_ellipse(aes(group=Group), type="t", level=0.9, linewidth=0.8) +
  geom_text_repel(aes(label=rownames(pca_input)), size=3) +
  scale_color_manual(values=c("#1b9e77", "#d95f02")) +
  labs(x=paste("PC1 (", var_exp[1], "%)", sep=""),
       y=paste("PC2 (", var_exp[2], "%)", sep=""),
       title="Batch-Corrected PCA Analysis") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")