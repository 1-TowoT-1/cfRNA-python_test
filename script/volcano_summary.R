library(ggplot2)
argv <- commandArgs(T)
input=argv[1]
main_tital=argv[2]
judge=argv[3]
outdir=argv[4]

# input="C:/Users/91361/Desktop/F_liver_cancer_vs_F_health_raw_edgeR.xls"
# outdir="C:/Users/91361/Desktop/"
# judge="edgeR"
# main_tital="F_liver_cancer_vs_F_health"


resdata <- read.table(input,header = T)
resdata <- na.omit(resdata)

if(grepl("edgeR",judge)){
threshold <- as.factor(ifelse(resdata$FDR < 0.05 &abs(resdata$logFC) >= 1 ,ifelse(resdata$logFC >= 1 ,'Up','Down'),'Not'))
summary_deg <- data.frame(table(threshold))
summary_deg <- rbind(summary_deg,data.frame(threshold="All",Freq=sum(summary_deg$Freq)))
write.table(summary_deg,paste(outdir,"/",main_tital,"_DEG_summary.xls",sep=""),row.names = F,quote=F,sep="\t")

p <- ggplot(resdata,aes(x=logFC,y=-log10(FDR),colour=threshold)) +
    xlab("log2(Fold Change)")+ylab("-log10(FDR)") +
    geom_point(size = 0.5,alpha=1) +
    ylim(0,max(-log10(resdata$FDR))) + xlim(-max(abs(resdata$logFC)),max(abs(resdata$logFC)))+
    scale_color_manual(values=c("green","grey", "red"))+
    geom_vline(xintercept = c(-log2(2), log2(2)), lty = 3, color = 'black') +  #添加阈值线
    geom_hline(yintercept = -log10(0.05), lty = 3, color = 'black')+
    ggtitle(main_tital)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
ggsave(paste(outdir,"/",main_tital,"_volcano.png",sep=""),p)
}


if(grepl("DEseq2",judge)){
  threshold <- as.factor(ifelse(resdata$padj < 0.05 &abs(resdata$log10FoldChange) >= 1 ,ifelse(resdata$log10FoldChange >= 1 ,'Up','Down'),'Not'))
  summary_deg <- data.frame(table(threshold))
  summary_deg <- rbind(summary_deg,data.frame(threshold="All",Freq=sum(summary_deg$Freq)))
  write.table(summary_deg,paste(outdir,"/",main_tital,"_DEG_summary.xls",sep=""),row.names = F,quote=F,sep="\t")
  

  p <- ggplot(resdata,aes(x=log2FoldChange,y=-log10(padj),colour=threshold)) +
    xlab("log2(Fold Change)")+ylab("-log10(padj)") +
    geom_point(size = 0.5,alpha=1) +
    ylim(0,max(-log10(resdata$padj))) + xlim(-max(abs(resdata$log2FoldChange)),max(abs(resdata$log2FoldChange)))+
    scale_color_manual(values=c("green","grey", "red"))+
    geom_vline(xintercept = c(-log2(0.5), log2(0.5)), lty = 3, color = 'black') +  #添加阈值线
    geom_hline(yintercept = -log10(0.05), lty = 3, color = 'black')+
    theme_bw()+
    ggtitle(main_tital)+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste(outdir,"/",main_tital,"_volcano.png",sep=""),p)
}