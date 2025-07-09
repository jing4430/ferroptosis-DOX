library(readxl)
library(ggplot2)
library(dplyr)

# ------------------------- 数据读取与预处理 -------------------------
# 读取 Excel 文件（确保列名匹配）
df <- read_excel("data_all.xls") %>%
  select(Gene_Symbol, log2FC, q_value) %>%
  na.omit() %>%
  mutate(log10q = -log10(q_value))

# ------------------------- 基因分类逻辑 -------------------------
df <- df %>%
  mutate(
    category = case_when(
      # 显著上调基因（双阈值满足）
      log2FC > 1 & q_value < 0.05 ~ "Sig_Up",
      # 显著下调基因（双阈值满足）
      log2FC < -1 & q_value < 0.05 ~ "Sig_Down",
      # 上调相关基因（满足单阈值且方向向上）
      (log2FC > 1 & q_value >= 0.05) | (q_value < 0.05 & log2FC > 0) ~ "Up_Related",
      # 下调相关基因（满足单阈值且方向向下）
      (log2FC < -1 & q_value >= 0.05) | (q_value < 0.05 & log2FC < 0) ~ "Down_Related",
      # 不显著基因
      TRUE ~ "Non_Sig"
    )
  )

# ------------------------- 颜色映射系统 -------------------------
color_palette <- c(
  "Sig_Up" = "#00BFFF",      # 天蓝色
  "Sig_Down" = "#FF69B4",    # 粉红色
  "Up_Related" = "#87CEFA",  # 浅天蓝
  "Down_Related" = "#FFB6C1",# 浅粉红
  "Non_Sig" = "#D3D3D3"      # 浅灰色
)

# ------------------------- 可视化实现 -------------------------
ggplot(df, aes(x = log2FC, y = log10q)) +
  # 散点图层（实心圆点）
  geom_point(
    aes(color = category),
    shape = 19,
    size = 2.5,
    alpha = 0.8
  ) +
  # 参考线系统
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  # 颜色映射
  scale_color_manual(
    name = "Expression Category",
    values = color_palette,
    labels = c(
      "Sig_Up" = "Significant Up (log2FC >1 & q <0.05)",
      "Sig_Down" = "Significant Down (log2FC < -1 & q <0.05)",
      "Up_Related" = "Up-regulated Features",
      "Down_Related" = "Down-regulated Features",
      "Non_Sig" = "Not Significant"
    )
  ) +
  # 坐标轴与标题
  labs(
    x = expression(bold(log[2]("Fold Change"))),
    y = expression(bold(-log[10]("q"))),
    title = "control vs DOX"
  ) +
  # 主题定制
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key = element_rect(fill = "white")
  )

# ------------------------- 图形保存 -------------------------
ggsave(
  "volcano_plot.png",
  width = 28,  # 厘米单位
  height = 18, 
  units = "cm",
  dpi = 600
)