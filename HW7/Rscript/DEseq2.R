library(DESeq2)

#准备输入数据
raw.counts <- read.table("count_exon.txt",sep='\t',header = T,row.names = 1)
mu.raw.counts <- raw.counts[,c("UD1_1","UD1_2","UD1_3","UD0_1","UD0_2","UD0_3")]
#过滤掉表达量低的基因
mu.filtered.counts <- mu.raw.counts[rowMeans(mu.raw.counts)>5,]

#提供样本条件信息
conditions <- factor(c(rep("Control",3),rep("Treatment",3)),levels = c("Control","Treatment"))
colData <- data.frame(row.names = colnames(mu.filtered.counts),conditions=conditions)

#差异分析
dds <- DESeqDataSetFromMatrix(mu.filtered.counts, colData, design = ~conditions)
dds <- DESeq(dds)
res <- results(dds)

#保存结果
#res <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.table(res,"deseq2.mutant.light.vs.dark.txt",sep = '\t',row.names = T,quote = F)