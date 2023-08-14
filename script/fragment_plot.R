argv <- commandArgs(T)
par(mai = c(1.2, 1.5, 1.2, 1))


data <- read.table(argv[1], header = F)
sample <- argv[2]
outdir <- argv[3]

# 设置插入片段长度的阈值，过滤掉太长的片段
length_cutoff <- 1200
fragment <- data$V1[data$V1 <= length_cutoff]
# 利用直方图统计频数分布，设置柱子个数
breaks_num <- 500
res <- hist(fragment, breaks = breaks_num, plot = FALSE)
# 添加坐标原点
png(paste(outdir,"/",sample,"_fragment.png",sep=""),width=1200,height=800)
par(mai = c(1.2, 1.5, 1.2, 1))
plot(x = c(0, res$breaks),
     y = c(0, 0, res$counts) / 10^2,
     type = "l", col = "red",
     xlab = "Fragment length(bp)",
     ylab = expression(Normalized ~ read ~ density ~ 10^2),
     main = paste(sample," Fragment sizes",sep = ""))
dev.off()

pdf(paste(outdir,"/",sample,"_fragment.pdf",sep=""))
par(mai = c(1.2, 1.5, 1.2, 1))
plot(x = c(0, res$breaks),
     y = c(0, 0, res$counts) / 10^2,
     type = "l", col = "red",
     xlab = "Fragment length(bp)",
     ylab = expression(Normalized ~ read ~ density ~ 10^2),
     main = paste(sample," Fragment sizes",sep = ""))
dev.off()
