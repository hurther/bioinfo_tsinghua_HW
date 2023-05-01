library(pheatmap)
library(edgeR)

#准备输入数据
raw.counts <- read.table("count_exon.txt",sep='\t',header = T,row.names = 1)
mu.raw.counts <- raw.counts[,c("UD1_1","UD1_2","UD1_3","UD0_1","UD0_2","UD0_3")]
edgerData <- read.table("edger.mutant.light.vs.dark.txt",sep = '\t',header = T,row.names = 1)
edgerData <- edgerData[edgerData$FDR < 0.05,]
#这里比较log2Foldchange比较实际值
edgerData <- edgerData[order(edgerData$logFC),]
target <- rbind(head(edgerData,n=10),tail(edgerData,n=10))
counts.matrix <- mu.raw.counts[rownames(target),] 

#计算log10CPM
y <- DGEList(counts = counts.matrix)
CPM.matrix <- edgeR::cpm(y,log=F)
log10.CPM.matrix <- log10(CPM.matrix+1)

#计算z-scores
z.scores <- (log10.CPM.matrix - rowMeans(log10.CPM.matrix))/apply(log10.CPM.matrix,1,sd)

#作heatmap
annotation_col = data.frame(Group = factor(c(rep("Control",3),rep("Treatment",3))))
rownames(annotation_col) = colnames(z.scores)
annotation_row = data.frame(Source = factor(c(rep("Tail",10),rep("Top",10))))
rownames(annotation_row) = rownames(z.scores)
ann_colors = list(Group = c(Control = "#7FBC41", Treatment = "#DE77AE"),
                  Source = c(Tail = "#807DBA",Top = "#BCBBDC"))
pheatmap(z.scores,cutree_cols = 2,
         annotation_col = annotation_col,annotation_colors = ann_colors,annotation_row = annotation_row,
         cluster_rows=F, show_rownames=TRUE, cluster_cols=TRUE)
