library(ggplot2)
library(dplyr)
library(scales)

args<-commandArgs(T)
pi <- read.table(args[1],head=T)
pi <- na.omit(pi)
str(pi)
pi$SNP1 <- seq(1, nrow(pi), 1)
pi$CHROM <- factor(pi$CHROM, levels = unique(pi$CHROM))
chr <- aggregate(pi$SNP1, by = list(pi$CHROM), FUN = median)
table(pi$CHROM)

str(pi)

library(RColorBrewer)
#
##special regions
pi_chr01 <- pi[pi$CHROM=="chr03",]
pi_chr01$CHROM <- factor(pi_chr01$CHROM, levels = unique(pi_chr01$CHROM))
pi_chr01 <- pi_chr01[pi_chr01[,3]>=6700000,]
pi_chr01 <- pi_chr01[pi_chr01[,3]<=9300000,]
str(pi_chr01)
min(pi_chr01$PBE)
max(pi_chr01$PBE)
#pi_chr01[pi_chr01[,3]==4877500,]$SNP1
#chr <- aggregate(pi_chr01$SNP1, by = list(pi_chr01$CHROM), FUN = median)

p <- ggplot(pi_chr01, aes(POS,PBE,colour=POP)) +
#p <- ggplot(pi_chr01, aes(POS,PBE)) +
  #设置点的大小，透明度
  geom_point(show.legend = FALSE,alpha=0.25, size=0.5) +
  geom_smooth(method = "loess",span=0.5, level = 0.95,alpha=0, linewidth = 0.4) +
  #设置颜色
  #scale_color_manual(values = brewer.pal(5,"Set1")) +
  #设定X轴
  scale_x_continuous(breaks = seq(6700000,9300000, 1000000),labels = label_currency(scale_cut = cut_short_scale()),expand = c(0, 0),limits=c(6700000, 9300000)) +
  #去除绘图区和X轴之间的gap
  #scale_y_continuous(breaks = seq(0, 11000, 5000), labels = as.character(seq(0, 11000, 5000)), expand = c(0, 0),
  #scale_y_continuous(breaks = seq(0, 3000, 1500), labels = as.character(seq(0, 3000, 1500)), expand = c(0, 0),
  scale_y_continuous(breaks = seq(-0.4, 0.45, 0.4), labels = as.character(seq(-0.4, 0.45, 0.4)), expand = c(0, 0),
  #limits = c(0, max(pi_chr01$PBE)+50)) +
  limits = c(-0.2, 0.46)) +
  #limits = c(0, max(pi_chr01$PBE)+0.002)) +
  geom_vline(xintercept = 8031484,col = "gray", lwd = 0.2, lty = 2, alpha=1)+
  geom_vline(xintercept = 8032179,col = "gray", lwd = 0.2, lty = 2, alpha=1)+
  #geom_vline(xintercept = 7734219,col = "gray", lwd = 0.2, lty = 2, alpha=1)+
  #geom_vline(xintercept = 7771853,col = "gray", lwd = 0.2, lty = 2, alpha=1)+
  #geom_vline(xintercept = 7793673,col = "gray", lwd = 0.2, lty = 2, alpha=1)+
  #geom_vline(xintercept = 7813981,col = "gray", lwd = 0.2, lty = 2, alpha=1)+
  theme(legend.title = element_text(size=4),legend.key.size = unit(0.3,"cm"),legend.text = element_text(size=6), panel.grid = element_blank(), axis.line = element_line(colour = 'black'), panel.background = element_rect(fill = 'transparent'),legend.position="right") +
  labs(x = 'chr01', y = expression('PBE')) +
  guides(colour = guide_legend(title = "Populations",nrow = 4, byrow = F))

#png(paste(args[1],"_sdhB_ggplot.png",sep=""),width = 2400, height = 900,res=400)pdf(paste(args[1],"_chr13_ggplot.pdf",sep=""),width = 4, height = 3)
pdf(paste(args[1],"_chr01_sdhD_ggplot3.pdf",sep=""),width = 5, height = 3)
p
dev.off()
