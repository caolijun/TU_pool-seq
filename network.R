library(ape)
library(pegas)
library("adegenet")
library(RColorBrewer)

nbin<-fasta2DNAbin("sdhd_all_noLab_BJMY.fas")
h<-pegas::haplotype(nbin,strict=FALSE,trailingGapsAsN=TRUE)
h
ng<-nbin
hname<-paste("H", 1:nrow(h), sep = "")
#nmname<-c(rep("BJ", 180), rep("SD", 165), rep("LabR", 24), rep("NMHH2", 16), rep("ZJXS1", 36), rep("YNKM", 43), rep("Others", 204))
#rownames(ng)<-nmname
hg<-pegas::haplotype(ng, strict = FALSE, labels = hname, trailingGapsAsN = TRUE) #extracts haplotypes from DNAbin object
hg
netg<-haploNet(hg, d = NULL, getProb = TRUE) #constructs the haplotype network
netg

ind.hapg<-with(
  utils::stack(setNames(attr(hg, "index"), rownames(hg))),
  table(hap=ind, individuals=rownames(ng)[values])
)

#brewer.pal(7,"Set2")

#gbg<-c(rep("#66C2A5"), 
#       rep("#FC8D62"), 
#       rep("#8DA0CB"), 
#       #rep("#E78AC3"), 
#       rep("#A6D854"), 
#       rep("#FFD92F"), 
#       rep("#E5C494"))
#,rep("darkseagreen1", 15))

hap_list <- as.data.frame(ind.hapg) ###Group colors in the legend command which gbg argument below, were assigned according to the order of names in the ind.hapg object.
#hapcol<-c("BJ", 
#          "LabR", 
#          "NMHH2", 
#          "Others", 
#          "SD", 
#          "YNKM", 
#         "ZJXS1")
#,"Termessos")


#ubg<-c(rep("#66C2A5",1), 
#       rep("#FC8D62",1), 
#       rep("#8DA0CB",1), 
#       #rep("#E78AC3",1),
#       rep("#A6D854",1),
#       rep("#FFD92F",1), 
#       rep("#E5C494",1))

#rep("olivedrab3",1),

pdf("tu_sdhd_all.pdf", width = 5, height = 8) #save as pdf file
par(mar=c(0.001,0.001,0.001, 0.001))
plot(netg, size=attr(netg, "freq"), scale.ratio = 90, cex = 0.8, labels=TRUE, show.mutation=0, font=2, fast=F)
#plot(netg, size=attr(netg, "freq"), bg = gbg, scale.ratio = 90, cex = 0.8, labels=TRUE, pie = ind.hapg, show.mutation=1, font=2, fast=TRUE)
#legend(x=-400,y=400, colnames(ind.hapg), fill=ubg, cex=0.8, ncol=1, bty = "n", x.intersp = 0.2)
dev.off()

write.table(hap_list,"SDHd_hap_list",row.names = F,quote = F)
