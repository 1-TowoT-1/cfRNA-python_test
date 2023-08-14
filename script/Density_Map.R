argv <-  commandArgs(T)
bam_name <- argv[3]

coverage_count <- read.table(argv[1],header = F)

coverage_count <- data.frame(coverage_count,"region"=coverage_count$V3-coverage_count$V2)
coverage_count <- coverage_count[,c(1,4,8)]

colnames(coverage_count) <- c("chr","counts","region")
coverage_count$normalize_count <- coverage_count$counts/coverage_count$region
coverage_count <- data.frame(coverage_count,"group"=rep(c(1:100),nrow(coverage_count)/100))


sum_counts <- aggregate(normalize_count ~ group, coverage_count, sum)
sum_counts$normalize_count <- sum_counts$normalize_count/nrow(coverage_count)

png(paste(argv[2],"/","Density_Map.png",sep=""),width=1200,height = 800)
plot(x=sum_counts$group,y=sum_counts$normalize_count,type='l',xlab="gene region",ylab="normalized_count",main=bam_name)
dev.off()
