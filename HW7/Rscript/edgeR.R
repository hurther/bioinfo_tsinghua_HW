library(edgeR)

#准备输入数据
raw.counts <- read.table("count_exon.txt",sep='\t',header = T,row.names = 1)
mu.raw.counts <- raw.counts[,c("UD1_1","UD1_2","UD1_3","UD0_1","UD0_2","UD0_3")]
#过滤掉表达量低的基因
mu.filtered.counts <- mu.raw.counts[rowMeans(mu.raw.counts)>5,]

#提供样本条件信息
conditions <- factor(c(rep("Control",3),rep("Treatment",3)),levels = c("Control","Treatment"))
design <- model.matrix(~conditions)

#差异分析
y <- DGEList(counts = mu.filtered.counts)
y <- calcNormFactors(y, method="TMM")
y <- estimateDisp(y,design = design)

fit <- glmFit(y,design = design)
lrt <- glmLRT(fit,coef=2)

#保存结果
diff.table <- topTags(lrt, n = nrow(y))$table
#diff.table <- diff.table[abs(diff.table$logFC) > 1 & diff.table$FDR < 0.05,]
write.table(diff.table,"edger.mutant.light.vs.dark.txt",sep = '\t',quote = F,row.names = T,col.names = T)
