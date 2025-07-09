library(ggplot2)
library(ggpubr)
library(magrittr)
library(ggsignif)
library(cowplot)
library(tidyverse)
library(ggsci)
library(gridExtra)



df = read.table(file='E:\\test\\results\\find_gene.txt',sep ='\t',header=TRUE)

##################################################################################################################
group_box <- function(gene=gene,df=df,length=length){
p <- ggboxplot(df, #数据对象
          x = 'group', # 选择x轴用那一列数据
          y = gene, #选择y轴用什么数据
          fill = 'group', #颜色根据哪一列决定
          bxp.errorbar = T, #是否添加error bar
          bxp.errorbar.width = 0.2, #error bar的长度
          palette = c('#EC6D60','#2BA7AA'), #颜色风格
          add = 'point', # 是否添加boxplot上面的点点
          
)+labs(
      x = gene, # x轴的名字
      y = 'FPKM' # y轴的名字
      )+
  geom_signif(comparisons = list(c('U', 'T')), # 设置要对比的组
             # y_position = c(34,36,38), #设置3个显著性标记的高度
              tip_length = length, #设置显著性那条横线两头向下的长度
            #  map_signif_level = T, #设置是否标记显著性的*号，还是直接标记数值
              test = t.test #设置显著性计算方式
  ) +
  theme(
    axis.text.x = element_text(color = 'black', size = 16, angle = 0),
    axis.text.y = element_text(color = 'black', size = 16, angle = 0),
    axis.title.x = element_text(color = 'black', size = 16, angle = 0),
    axis.title.y = element_text(color = 'black', size = 16, angle = 90),
    legend.title = element_text(color = 'black', size = 16),
    legend.text = element_text(color = 'black', size = 16),
    axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
    axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
    panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
  )}




p1<-group_box(gene='GCH1',df=df,length=c(0.5,0.03))
p2<-group_box(gene='AKR1C3',df=df,length=c(0.02,0.6))
p3<-group_box(gene='NR1H4',df=df,length=c(0.02,0.6))
p4<-group_box(gene='UCP1',df=df,length=c(0.02,0.6))
p5<-group_box(gene='POLQ',df=df,length=c(0.02,0.6))

grid.arrange(p1,p2,p3,p4,p5,ncol=3)
