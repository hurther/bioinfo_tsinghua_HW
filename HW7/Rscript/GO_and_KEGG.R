data <- read.table("wt.light.vs.dark.all.txt",header = T,sep = '\t')
data.filtered <- subset(data, log2FoldChange > 1 & data$padj < 0.05)
genes <- rownames(data.filtered)
write.table(genes,"GO_and_KEGG.txt",sep = "\n",quote = F,col.names = F,row.names = F)
