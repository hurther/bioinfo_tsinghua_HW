library(VennDiagram)

#载入数据
deseq2Data <- read.table("deseq2.mutant.light.vs.dark.txt",header = T,sep = '\t')
deseq2Data <- subset(deseq2Data,padj < 0.05 & abs(log2FoldChange) > 1)
edgerData <- read.table("edger.mutant.light.vs.dark.txt",header = T,sep = '\t')
edgerData <- edgerData[abs(edgerData$logFC) > 1 & edgerData$FDR < 0.05,]
group1 <- rownames(deseq2Data)
group2 <- rownames(edgerData)

#绘制Venn图
venn.plot <- venn.diagram(
  x = list(
    DEseq2=group1,
    edgeR=group2         
  ),
  filename = "DEseq2.vs.edgeR.Venn.png", imagetype = "png",
  lwd = 3, #边框线宽度
  col = "transparent",
  fontface = "bold", #标签字体
  fill = c("green","yellow"), #填充色
  alpha = 0.6, #透明度
  cex = 2, #标签字体大小
  cat.cex = 1.5, #类名字体大小
  cat.fontface = "bold", #类名字体
  margin = 0.05, #边际距离
  cat.dist = c(0.04, 0.04),
  cat.pos = c(-20, 20)
)
