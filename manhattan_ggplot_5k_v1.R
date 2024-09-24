library(ggplot2)
library(dplyr)

args<-commandArgs(T)
#data <- read.table(args[1],head=T)
data <- read.table(paste(args[1],"_5k.fst",sep=""),head=T)
#data$POS <- (data[,2]+data[,3])/2
RS1 <- data$RS1
RS2 <- data$RS2
S1S2 <- data$S1S2

TRS1 <- -log(1-RS1, base = 10)
TRS2 <- -log(1-RS2, base = 10)
TS1S2 <- -log(1-S1S2, base = 10)

PBS <- (TRS1 + TRS2 - TS1S2)/2
PBS_m <- median(PBS,na.rm=T)

TS1S2_m <- median(TS1S2)

PBE <- PBS - (TS1S2*PBS_m/TS1S2_m)
PBE_0 <- data.frame(PBE)
PBE1 <- cbind(data[,c(1,2)],PBE_0[,1])
colnames(PBE1)[3] <- "PBE"
str(PBE1)
write.table(PBE1,paste(args[1],"filtered_PBE",sep=""),sep = "\t",row.names = F,quote = F)

PBE1[PBE1 == "Inf"] = max(PBE1[,3] < "Inf",na.rm=T) 
PBE1[PBE1 == "-Inf"] = min(PBE1[,3] > "-Inf",na.rm=T) 

PBE.thresh=quantile(PBE1[,3],probs=0.99,na.rm=T)
PBE.thresh
max(PBE1[,3],na.rm=T)
min(PBE1[,3],na.rm=T)
outliers <- PBE1[PBE1[,3]>=PBE.thresh,1:3]
str(outliers)
#write.table(outliers,paste(args[1],"_PBE_0.99_outliers"),sep = "\t",row.names = F,quote = F)
options(scipen = 200)
outliers[,4] <-outliers[,2]-5
outliers[,5] <-outliers[,2]+5
outliers1 <- outliers[,c(1,4,5,3)]
outliers1 = na.omit(outliers1)
write.table(outliers1,paste(args[1],"_PBE_5k0.99.range.pbe",sep=""),sep="\t",quote=F,col.name=F,row.name=F)
write.table(outliers1[,c(1,2,3)],paste(args[1],"_PBE_5k0.99.range",sep=""),sep="\t",quote=F,col.name=F,row.name=F)

PBE1$SNP1 <- seq(1, nrow(PBE1), 1)
PBE1$CHROM <- factor(PBE1$CHROM, levels = unique(PBE1$CHROM))
chr <- aggregate(PBE1$SNP1, by = list(PBE1$CHROM), FUN = median)
str(PBE1)
write.table(PBE1,paste(args[1],"filtered_PBE_ID",sep=""),sep = "\t",row.names = F,quote = F)

p <- ggplot(PBE1, aes(SNP1, PBE)) +
  #annotate('rect', xmin = 0, xmax = max(PBE1$SNP1), ymin = min(PBE1[,4],na.rm=T)-0.02, ymax = max(PBE1[,4],na.rm=T)+0.02, fill = 'gray') +
  #设置点的大小，透明度
  geom_point(aes(color = CHROM), show.legend = FALSE,alpha=0.8, size=0.7) +
  #设置颜色
  scale_color_manual(values = rep(c("#E78AC3", "skyblue"), 3)) +
  #设定X轴
  scale_x_continuous(breaks = chr$x, labels = chr$Group.1, expand = c(0, 0))+ #,limits=c(21626000,21816000)) +
  #去除绘图区和X轴之间的gap
  scale_y_continuous(breaks = seq(-2, 2, 0.5), labels = as.character(seq(-2, 2, 0.5)), expand = c(0, 0), 
  #scale_y_continuous(breaks = seq(-2, 2, 0.1), labels = as.character(seq(-2, 2, 0.1)), expand = c(0, 0), 
  limits = c(-0.72, 0.67)) +
  #limits = c(min(PBE1$PBE,na.rm=T)-0.02, max(PBE1$PBE,na.rm=T)+0.02)) +
 # limits = c(-0.05, max(PBE1$PBE,na.rm=T)+0.02)) +
  geom_hline(yintercept = PBE.thresh,col = "black", lwd = 0.2, lty = 2, alpha=1)+
  
  #geom_vline(xintercept = 126894,col = "black", lwd = 0.2, lty = 2, alpha=1)+
  #geom_vline(xintercept = 648433,col = "black", lwd = 0.2, lty = 2, alpha=1)+
  #geom_vline(xintercept = 1694684,col = "black", lwd = 0.2, lty = 2, alpha=1)+
  #geom_vline(xintercept = 1932460,col = "black", lwd = 0.2, lty = 2, alpha=1)+
  geom_vline(xintercept = 716,col = "black", lwd = 0.2, lty = 2, alpha=1)+
  geom_vline(xintercept = 3590,col = "black", lwd = 0.2, lty = 2, alpha=1)+
  geom_vline(xintercept = 9368,col = "black", lwd = 0.2, lty = 2, alpha=1)+
  geom_vline(xintercept = 10521,col = "black", lwd = 0.2, lty = 2, alpha=1)+

  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), panel.background = element_rect(fill = 'transparent')) +
  labs(x = 'Chromosome', y = expression('PBE'),axis.text.x = element_text(size=11),axis.text.y = element_text(size=11))

#pdf(paste(args[1],"_PBE_5k_ggplot.pdf",sep=""),width = 16, height = 3)
png(paste(args[1],"_PBE_5k_ggplot.png",sep=""),width = 4000, height = 750,res=400)
p
dev.off()

##special regions
#PBE2 <- PBE1[PBE1[,1]=="chr01",]
#str(PBE2)
#PBE2$SNP1 <- seq(1, nrow(PBE2), 1)
#PBE2$CHROM <- factor(PBE2$CHROM, levels = unique(PBE2$CHROM))
#PBE2 <- PBE2[PBE2[,2]<=5178000,]
#PBE2 <- PBE2[PBE2[,2]>=4577000,]
#str(PBE2)
#chr <- aggregate(PBE2$SNP1, by = list(PBE2$CHROM), FUN = median)
#
#p <- ggplot(PBE2, aes(POS, PBE)) +
#  #annotate('rect', xmin = 0, xmax = max(PBE1$SNP1), ymin = min(PBE1[,4],na.rm=T)-0.02, ymax = max(PBE1[,4],na.rm=T)+0.02, fill = 'gray') +
#  geom_vline(xintercept = 4877168,col = "gray", lwd = 0.2, lty = 2, alpha=1)+
#  geom_vline(xintercept = 4878068,col = "gray", lwd = 0.2, lty = 2, alpha=1)+
#  #设置点的大小，透明度
#  geom_point(aes(color = CHROM), show.legend = FALSE,alpha=0.8, size=1) +
#  #设置颜色
#  scale_color_manual(values = rep(c("#E78AC3", "skyblue"), 3)) +
#  #设定X轴
#  scale_x_continuous(expand = c(0, 0),limits=c(4577000,5178000)) +
#  #去除绘图区和X轴之间的gap
#  scale_y_continuous(breaks = seq(-1, 1, 0.2), labels = as.character(seq(-1, 1, 0.2)), expand = c(0, 0), limits = c(min(PBE2$PBE,na.rm=T)-0.02, max(PBE2$PBE,na.rm=T)+0.02)) +
#  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), panel.background = element_rect(fill = 'transparent')) +
#  labs(x = 'Chromosome', y = expression('PBE'))
#
#pdf(paste(args[1],"_PBE_10k2k_chr01_sdhB.pdf",sep=""),width = 3.8, height = 4)
#p
#dev.off()
#
#PBE2 <- PBE1[PBE1[,1]=="chr01",]
#PBE2$SNP1 <- seq(1, nrow(PBE2), 1)
#PBE2$CHROM <- factor(PBE2$CHROM, levels = unique(PBE2$CHROM))
#PBE2 <- PBE2[PBE2[,2]<=22100000,]
#PBE2 <- PBE2[PBE2[,2]>=21600000,]
#str(PBE2)
#chr <- aggregate(PBE2$SNP1, by = list(PBE2$CHROM), FUN = median)
#
#p <- ggplot(PBE2, aes(POS, PBE)) +
#  #annotate('rect', xmin = 0, xmax = max(PBE1$SNP1), ymin = min(PBE1[,4],na.rm=T)-0.02, ymax = max(PBE1[,4],na.rm=T)+0.02, fill = 'gray') +
#  geom_vline(xintercept = 21756216,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21757429,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21764043,col = "blue", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21764821,col = "blue", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21769050,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21770562,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21811467,col = "blue", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21812973,col = "blue", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21872146,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21873796,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21885433,col = "blue", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21889529,col = "blue", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21891351,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21893959,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21894157,col = "blue", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21899061,col = "blue", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21899462,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21899909,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21903597,col = "blue", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21905103,col = "blue", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21906199,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21907711,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  #设置点的大小，透明度
#  geom_point(aes(color = CHROM), show.legend = FALSE,alpha=0.8, size=1) +
#  #设置颜色
#  scale_color_manual(values = rep(c("#E78AC3", "skyblue"), 3)) +
#  #设定X轴
#  scale_x_continuous(breaks = c(seq(21600000,22100000,by=1000000)),expand = c(0, 0),limits=c(21600000,22100000)) +
#  #去除绘图区和X轴之间的gap
#  scale_y_continuous(breaks = seq(-1, 1, 0.2), labels = as.character(seq(-1, 1, 0.2)), expand = c(0, 0), limits = c(min(PBE2$PBE,na.rm=T)-0.02, max(PBE2$PBE,na.rm=T)+0.02)) +
#  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), panel.background = element_rect(fill = 'transparent')) +
#  labs(x = 'Chromosome', y = expression('PBE'))
#
#pdf(paste(args[1],"_PBE_10k2k_chr01_P450.pdf",sep=""),width = 3.8, height = 4)
#p
#dev.off()
#
#PBE2 <- PBE1[PBE1[,1]=="chr02",]
#PBE2$SNP1 <- seq(1, nrow(PBE2), 1)
#PBE2$CHROM <- factor(PBE2$CHROM, levels = unique(PBE1$CHROM))
#PBE2 <- PBE2[PBE2[,2]<=15400000,]
#PBE2 <- PBE2[PBE2[,2]>=14600000,]
#str(PBE2)
#chr <- aggregate(PBE2$SNP1, by = list(PBE2$CHROM), FUN = median)
#
#p <- ggplot(PBE2, aes(POS, PBE)) +
#  #annotate('rect', xmin = 0, xmax = max(PBE1$SNP1), ymin = min(PBE1[,4],na.rm=T)-0.02, ymax = max(PBE1[,4],na.rm=T)+0.02, fill = 'gray') +
#  geom_vline(xintercept = 14882795,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 14884376,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 14888405,col = "blue", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 14890009,col = "blue", lwd = 0.1, lty = 2, alpha=1)+
#  #设置点的大小，透明度
#  geom_point(aes(color = CHROM), show.legend = FALSE,alpha=0.8, size=1) +
#  #设置颜色
#  scale_color_manual(values = rep(c( "skyblue"), 3)) +
#  #设定X轴
#  scale_x_continuous(breaks = c(seq(14600000,15400000,by=1000000)),expand = c(0, 0),limits=c(14600000,15400000)) +
#  #去除绘图区和X轴之间的gap
#  scale_y_continuous(breaks = seq(-1, 1, 0.2), labels = as.character(seq(-1, 1, 0.2)), expand = c(0, 0), limits = c(min(PBE2$PBE,na.rm=T)-0.02, max(PBE2$PBE,na.rm=T)+0.02)) +
#  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), panel.background = element_rect(fill = 'transparent')) +
#  labs(x = 'Chromosome', y = expression('PBE'))
#
#pdf(paste(args[1],"_PBE_10k2k_chr02_SH.pdf",sep=""),width = 3.8, height = 4)
#p
#dev.off()
#
#PBE2 <- PBE1[PBE1[,1]=="chr02",]
#PBE2$SNP1 <- seq(1, nrow(PBE2), 1)
#PBE2$CHROM <- factor(PBE2$CHROM, levels = unique(PBE1$CHROM))
#PBE2 <- PBE2[PBE2[,2]<=22200000,]
#PBE2 <- PBE2[PBE2[,2]>=21000000,]
#str(PBE2)
#
#p <- ggplot(PBE2, aes(POS, PBE)) +
#  #annotate('rect', xmin = 0, xmax = max(PBE1$SNP1), ymin = min(PBE1[,4],na.rm=T)-0.02, ymax = max(PBE1[,4],na.rm=T)+0.02, fill = 'gray') +
#  geom_vline(xintercept = 21611323,col = "blue", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21617966,col = "blue", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21751027,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21752720,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21787481,col = "blue", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21801311,col = "blue", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21822594,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 21825383,col = "red", lwd = 0.1, lty = 2, alpha=1)+
#  #设置点的大小，透明度
#  geom_point(aes(color = CHROM), show.legend = FALSE,alpha=0.8, size=1) +
#  #设置颜色
#  scale_color_manual(values = rep(c("skyblue"), 3)) +
#  #设定X轴
#  scale_x_continuous(breaks = c(seq(21000000,22200000,by=1000000)),expand = c(0, 0),limits=c(21000000,22200000)) +
#  #去除绘图区和X轴之间的gap
#  scale_y_continuous(breaks = seq(-1, 1, 0.2), labels = as.character(seq(-1, 1, 0.2)), expand = c(0, 0), limits = c(min(PBE2$PBE,na.rm=T)-0.02, max(PBE2$PBE,na.rm=T)+0.02)) +
#  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), panel.background = element_rect(fill = 'transparent')) +
#  labs(x = 'Chromosome', y = expression('PBE'))
#
#pdf(paste(args[1],"_PBE_10k2k_chr02_CE.pdf",sep=""),width = 3.8, height = 4)
#p
#dev.off()
#
#PBE2 <- PBE1[PBE1[,1]=="chr03",]
#PBE2$SNP1 <- seq(1, nrow(PBE2), 1)
#PBE2$CHROM <- factor(PBE2$CHROM, levels = unique(PBE1$CHROM))
#PBE2 <- PBE2[PBE2[,2]<=8280000,]
#PBE2 <- PBE2[PBE2[,2]>=7780000,]
#str(PBE2)
#chr <- aggregate(PBE2$SNP1, by = list(PBE2$CHROM), FUN = median)
#
#p <- ggplot(PBE2, aes(POS, PBE)) +
#  #annotate('rect', xmin = 0, xmax = max(PBE1$SNP1), ymin = min(PBE1[,4],na.rm=T)-0.02, ymax = max(PBE1[,4],na.rm=T)+0.02, fill = 'gray') +
#  geom_vline(xintercept = 8031484,col = "gray", lwd = 0.1, lty = 2, alpha=1)+
#  geom_vline(xintercept = 8032179,col = "gray", lwd = 0.1, lty = 2, alpha=1)+
#  #设置点的大小，透明度
#  geom_point(aes(color = CHROM), show.legend = FALSE,alpha=0.8, size=1) +
#  #设置颜色
#  scale_color_manual(values = rep(c("#E78AC3"), 3)) +
#  #设定X轴
#  scale_x_continuous(breaks = c(seq(7780000,8280000,by=1000000)),expand = c(0, 0),limits=c(7780000,8280000)) +
#  #去除绘图区和X轴之间的gap
#  scale_y_continuous(breaks = seq(-1, 1, 0.2), labels = as.character(seq(-1, 1, 0.2)), expand = c(0, 0), limits = c(min(PBE2$PBE,na.rm=T)-0.02, max(PBE2$PBE,na.rm=T)+0.02)) +
#  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), panel.background = element_rect(fill = 'transparent')) +
#  labs(x = 'Chromosome', y = expression('PBE'))
#
#pdf(paste(args[1],"_PBE_10k2k_chr03_sdhD.pdf",sep=""),width = 3.8, height = 4)
#p
#dev.off()
#
