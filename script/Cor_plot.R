library(pheatmap)

argv <- commandArgs(T)
count <- argv[1]
output <- argv[2]

count_file <- read.table(count,header = T,row.names =1)
matrix <- cor(count_file)

png(paste(output,"/","All_sample_Cor.png",sep=""),width = 800,height = 800)
pheatmap(matrix,cluster_rows=F,cluster_cols = F,display_numbers = T,fontsize = 12)
dev.off()

png(paste(output,"/","All_sample_Cor.pdf",sep=""))
pheatmap(matrix,cluster_rows=F,cluster_cols = F,display_numbers = T,fontsize = 12)
dev.off()
