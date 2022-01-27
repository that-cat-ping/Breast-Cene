rm(list=ls())
library(DESeq2)

#设置工作目录

setwd("C:\\Users\\thecat\\Desktop\\thecat")

#表达量数据
mycount <- read.csv("差异分析.csv",header = T, row.names = 1)
mycount<-t(mycount)

#分组数据
group<-read.csv("label.csv",header = T, row.names = 1)
attach(group)
coldata <- data.frame(row.names = colnames(mycount), condition)

# 正式构建dds矩阵
dds <- DESeqDataSetFromMatrix(mycount,coldata, design = ~ condition)
dds <- DESeq(dds)  #对原始dds进行normalize
#显示dds信息
dds

# 查看结果的名称，本次实验中是 "Intercept"，"condition_akap95_vs_control"
resultsNames(dds)
# 将结果用results()函数来获取，赋值给res变量
#res = results(dds, contrast=c("stage1", "early", "later"))
res = results(dds, contrast=c("condition", "turmor", "paracancerous"))
head(res)
# summary一下，看一下结果的概要信息
summary(res)
## 输出图片
plotMA(res)	#画火山图，横轴是标准化后的平均readscount，纵轴是差异倍数，大于0是上调，小于0是下调,红色点表示显著差异的基因

plotMA(res, alpha = 0.05, colSig = 'red', colLine = 'skyblue')
# alpha:p-value
# colSig:显著性基因的颜色
# colLine:y=0的水平线

#write.csv(res, file="1vs3_DESeq_results.csv")
write.csv(res, file="2vs2_DESeq_results.csv")
diff_gene_Group <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

dim(diff_gene_Group)
#显示结果的首信息
head(diff_gene_Group)
#将结果进行输出保存
write.csv(diff_gene_Group, file = "./differential_gene1_1vs3.csv")




filter_up <- subset(res, padj < 0.05 & log2FoldChange > 1) #过滤上调基因

filter_down <- subset(res, padj < 0.05 & log2FoldChange < -1) #过滤下调基因
print(paste('差异上调基因数量: ', nrow(filter_up)))  #打印上调基因数量
print(paste('差异下调基因数量: ', nrow(filter_down)))  #打印下调基因数量


write.csv(filter_up, file="./filter_up_gene.csv", quote = F)  
write.csv(filter_down, file="./filter_down_gene.csv", quote = F)

