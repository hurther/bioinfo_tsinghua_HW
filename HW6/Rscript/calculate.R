counts.matrix <- read.table("D:/大二下/生物信息学/bioinfo_tsinghua_HW/HW6/output/counts_matrix.txt")
counts.matrix <- counts.matrix[,-c(1)]

CPM.matrix <- t(1000000*t(counts.matrix)/colSums(counts.matrix))

log10.CPM.matrix <- log10(CPM.matrix+1) 
# 1为pseudocount, 避免count为0时对数未定义的情况

z.scores <- (log10.CPM.matrix - rowMeans(log10.CPM.matrix))/apply(log10.CPM.matrix,1,sd)
# apply(log10.CPM.matrix,1,sd)表示计算每行(1表示行,2表示列)的标准差(sd函数)
# rowMeans(log10.CPM.matrix)和apply(log10.CPM.matrix,1,mean)效果是一样的

# write.table(z.scores,"D:/大二下/生物信息学/bioinfo_tsinghua_HW/HW6/output/z_scores.txt",quote = F)