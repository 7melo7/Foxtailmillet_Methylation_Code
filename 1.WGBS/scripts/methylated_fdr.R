library(graphics)
args=commandArgs(T)
file1 <- args[1]
file2 <- args[2]
data <- read.table(file1,header=F)
p <- data$V8
fdr=p.adjust(p,"BH")
data$V9 <- fdr
write.table(data,file=file2,quote=F,sep="\t",row.names=F,col.names=F)
