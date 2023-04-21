library(pheatmap)
df = read.table("D:/大二下/生物信息学/bioinfo_tsinghua_HW/HW6/output/z_scores.txt")
df[is.na(df[,1:ncol(df)])]=0  # NaN值置为0
df[df > 2]=2
df[df < -2]=-2
annotation_col = data.frame(CellType = factor(rep(c("COAD", "ESCA" ,"READ"), c(50,50,50))))
rownames(annotation_col) = colnames(df)
ann_colors = list(CellType = c(COAD = "#7FBC41", ESCA = "#DE77AE", READ = "#807DBA"))
pheatmap(df,cutree_cols = 3,annotation_col = annotation_col,annotation_colors = ann_colors,
       cluster_rows=F, cluster_cols=F, show_rownames=F,show_colnames = F)