# 处理COAD中的数据
setwd("D:/CODE/R_projects/tumor-transcriptome-demo/COAD/")
file_list <- dir()
COAD_data <- data.frame()
order = 1

for (i in file_list){
  df <- read.table(i, header = T,sep = "")
  colnames(df)[7] <- paste("COAD",as.character(order),sep="")
  df <- df[,-(2:6)]
  if (order == 1){
    COAD_data <- df
  }
  else{
    COAD_data <- merge(COAD_data,df,by="Geneid")
  }
  order=order+1
}

# 处理ESCA的数据
setwd("D:/CODE/R_projects/tumor-transcriptome-demo/ESCA/")
file_list <- dir()
ESCA_data <- data.frame()
order = 1

for (i in file_list){
  df <- read.table(i, header = T,sep = "")
  colnames(df)[7] <- paste("ESCA",as.character(order),sep="")
  df <- df[,-(2:6)]
  if (order == 1){
    ESCA_data <- df
  }
  else{
    ESCA_data <- merge(ESCA_data,df,by="Geneid")
  }
  order=order+1
}

# 处理READ的数据
setwd("D:/CODE/R_projects/tumor-transcriptome-demo/READ/")
file_list <- dir()
READ_data <- data.frame()
order = 1

for (i in file_list){
  df <- read.table(i, header = T,sep = "")
  colnames(df)[7] <- paste("READ",as.character(order),sep="")
  df <- df[,-(2:6)]
  if (order == 1){
    READ_data <- df
  }
  else{
    READ_data <- merge(READ_data,df,by="Geneid")
  }
  order=order+1
}
rm(df,order)

counts.matrix <- merge(COAD_data,ESCA_data,by="Geneid")
counts.matrix <- merge(counts.matrix,READ_data,by="Geneid")
rm(COAD_data,ESCA_data,READ_data)
write.table(counts.matrix,
            "D:/大二下/生物信息学/bioinfo_tsinghua_HW/HW6/output/counts_matrix.txt",
            quote = FALSE)
