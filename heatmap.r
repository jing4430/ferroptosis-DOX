# 加载必要包
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

# ------------------------- 数据读取与预处理 -------------------------
# 读取Excel文件（假设基因名为第一列，样本为后续列）
fpkm_data <- read_excel("sig_data.xls") %>%
  tibble::column_to_rownames("Gene_Symbol")  # 基因名设为行名

# 定义分组信息（假设前3列Control，后3列DOX）
group_col <- data.frame(Group = factor(rep(c("Control", "DOX"), each = 13)))
rownames(group_col) <- colnames(fpkm_data)

# 过滤低表达基因（保留在至少50%样本中FPKM≥1的基因）
#keep_genes <- rowSums(fpkm_data >= 1) >= 0.5 * ncol(fpkm_data)
#fpkm_filtered <- fpkm_data[keep_genes, ]

# 数据标准化（Z-score按行标准化）
fpkm_zscore <- t(scale(t(fpkm_data)))

# ------------------------- 热图绘制 -------------------------
# 定义颜色方案
heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)  # 红蓝渐变
group_colors <- list(Group = c(Control = "#4DAF4A", DOX = "#E41A1C"))  # 分组颜色

# 生成热图
pheatmap(
  mat = fpkm_zscore,
  color = heatmap_colors,
  cluster_rows = TRUE,        # 基因聚类
  cluster_cols = FALSE,        # 关闭样本聚类（保持分组顺序）
  show_rownames = FALSE,       # 不显示所有基因名
  show_colnames = TRUE,        # 显示样本名
  annotation_col = group_col,  # 列注释条
  annotation_colors = group_colors,
  border_color = NA,           # 移除单元格边框
  fontsize_col = 10,           # 列名字体大小
  main = "DEGs Expression Heatmap (Control vs DOX)",
  filename = "heatmap.pdf",    # 输出PDF矢量图
  width = 8,                   # 图像宽度（英寸）
  height = 10                  # 图像高度（英寸）
)