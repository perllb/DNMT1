## Get closest genes

setwd(dir = "~/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/2018.08.30_DNMT1L1_figures_scripts/NearbyGenes/")

##### Set up environment #####
# path to directory with featurecount output
#rm(list=ls())
#gc()

# load deseqAbstraction
library(devtools)
install_github("perllb/deseqAbstraction")
library(deseqAbstraction)
library(DESeq2)
library(ggplot2)
#######
##### initialize the rnaseq object #####
path <- "~/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/RNAseq/PairedParam/Quant/hg38.gencode.exon.primary.txt"

# make deseqAbs object, give name of the featurecount file
colData <- data.frame(condition=c(rep("CTR",3),rep("KO",3)),
                      samples=c(107,108,109,110,111,112))

dnmt <- deseqAbs$new(name="DNMT1KO",filename=path,colData=colData)
dnmt$fullAuto()

#########

### Get close genes #####
bedfile <- "/Users/perludvik/Documents/bioinformatics/genomicData/hg38/NCBI_iGenomes/NCBI_genes.bestRefSeq.transcript.bed"
bed <- read.delim(bedfile,header=F,col.names = c("Chr","Start","End","ID","score","Strand"))
head(bed)

l1hs.up <- read.delim(file = "features/L1HS.up.padj05.sense.bed",header=F,col.names = c("Chr","Start","End","ID","score","Strand"))
head(l1hs.up)
l1pa2.up <- read.delim(file = "features/L1PA2.up.padj05.sense.bed",header=F,col.names = c("Chr","Start","End","ID","score","Strand"))
head(l1pa2.up)
l1pa3.up <- read.delim(file = "features/L1PA3.up.padj05.sense.bed",header=F,col.names = c("Chr","Start","End","ID","score","Strand"))
head(l1pa3.up)
herv <- read.delim(file = "features/hg38.HERV.score300.bed",header=F,col.names = c("Chr","Start","End","ID","Strand","score"))
head(herv)
sva <- read.delim(file = "features/SVA.hg38.uniqID.bed",header=F,col.names = c("Chr","Start","End","ID","score","Strand"))
head(sva)


d <- 50000

plotAll <- function(d) {
    
  if(!dir.exists(paste("output/",d,"/",sep = ""))) {
    dir.create(paste("output/",d,"/",sep = ""))
  }
  
  l1hs.close <- closeGenes(a=bed,l1hs.up,d=d)
  l1pa2.close <- closeGenes(a=bed,l1pa2.up,d=d)
  l1pa3.close <- closeGenes(a=bed,l1pa3.up,d=d)
  herv.close <- closeGenes(a=bed,herv,d=d)
  sva.close <- closeGenes(a=bed,sva,d=d)
  
  # make bed for deeptools analysis: chr - TSS leftmost- TSS rightmost - IDleft:IDright - "." - strand Gene
  
  ## L1HS
  l1hs.close.bed <- data.frame(chr=l1hs.close$A_chr,
                               tssLeft=ifelse(test=l1hs.close$B_TSS>l1hs.close$A_TSS,yes = l1hs.close$A_TSS,no = l1hs.close$B_TSS),
                               tssRight=ifelse(test=l1hs.close$B_TSS<l1hs.close$A_TSS,yes = l1hs.close$A_TSS,no = l1hs.close$B_TSS),
                               id=ifelse(test=l1hs.close$B_TSS<l1hs.close$A_TSS,yes = paste(l1hs.close$B_ID,l1hs.close$A_ID,sep = ":"),no = paste(l1hs.close$A_ID,l1hs.close$B_ID,sep = ":")),
                               score=".",
                               strandgene=ifelse(test=l1hs.close$B_TSS<l1hs.close$A_TSS,yes = "+",no = "-"))
  
  
  
  write.table(x = l1hs.close.bed,file = paste("output/",d,"/L1HS.close.",d,".bed",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(x = data.frame(chr=l1hs.close$A_chr,start=l1hs.close$B_TSS-1,end=l1hs.close$B_TSS,id=l1hs.close$B_ID,score=".",strand=l1hs.close$B_Strand),file = paste("output/",d,"/L1HS.close.",d,".TEonly.bed",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(x = data.frame(chr=l1hs.close$A_chr,start=l1hs.close$A_TSS-1,end=l1hs.close$A_TSS,id=l1hs.close$A_ID,score=".",strand=l1hs.close$A_Strand),file = paste("output/",d,"/L1HS.close.",d,".GeneOnly.bed",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  
  ## L1PA2
  l1pa2.close.bed <- data.frame(chr=l1pa2.close$A_chr,
                               tssLeft=ifelse(test=l1pa2.close$B_TSS>l1pa2.close$A_TSS,yes = l1pa2.close$A_TSS,no = l1pa2.close$B_TSS),
                               tssRight=ifelse(test=l1pa2.close$B_TSS<l1pa2.close$A_TSS,yes = l1pa2.close$A_TSS,no = l1pa2.close$B_TSS),
                               id=ifelse(test=l1pa2.close$B_TSS<l1pa2.close$A_TSS,yes = paste(l1pa2.close$B_ID,l1pa2.close$A_ID,sep = ":"),no = paste(l1pa2.close$A_ID,l1pa2.close$B_ID,sep = ":")),
                               score=".",
                               strandgene=ifelse(test=l1pa2.close$B_TSS<l1pa2.close$A_TSS,yes = "+",no = "-"))
  
  
  
  write.table(x = l1pa2.close.bed,file = paste("output/",d,"/l1pa2.close.",d,".bed",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(x = data.frame(chr=l1pa2.close$A_chr,start=l1pa2.close$B_TSS-1,end=l1pa2.close$B_TSS,id=l1pa2.close$B_ID,score=".",strand=l1pa2.close$B_Strand),file = paste("output/",d,"/l1pa2.close.",d,".TEonly.bed",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(x = data.frame(chr=l1pa2.close$A_chr,start=l1pa2.close$A_TSS-1,end=l1pa2.close$A_TSS,id=l1pa2.close$A_ID,score=".",strand=l1pa2.close$A_Strand),file = paste("output/",d,"/l1pa2.close.",d,".GeneOnly.bed",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  
  ## L1PA3
  l1pa3.close.bed <- data.frame(chr=l1pa3.close$A_chr,
                               tssLeft=ifelse(test=l1pa3.close$B_TSS>l1pa3.close$A_TSS,yes = l1pa3.close$A_TSS,no = l1pa3.close$B_TSS),
                               tssRight=ifelse(test=l1pa3.close$B_TSS<l1pa3.close$A_TSS,yes = l1pa3.close$A_TSS,no = l1pa3.close$B_TSS),
                               id=ifelse(test=l1pa3.close$B_TSS<l1pa3.close$A_TSS,yes = paste(l1pa3.close$B_ID,l1pa3.close$A_ID,sep = ":"),no = paste(l1pa3.close$A_ID,l1pa3.close$B_ID,sep = ":")),
                               score=".",
                               strandgene=ifelse(test=l1pa3.close$B_TSS<l1pa3.close$A_TSS,yes = "+",no = "-"))
  
  
  
  write.table(x = l1pa3.close.bed,file = paste("output/",d,"/l1pa3.close.",d,".bed",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(x = data.frame(chr=l1pa3.close$A_chr,start=l1pa3.close$B_TSS-1,end=l1pa3.close$B_TSS,id=l1pa3.close$B_ID,score=".",strand=l1pa3.close$B_Strand),file = paste("output/",d,"/l1pa3.close.",d,".TEonly.bed",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(x = data.frame(chr=l1pa3.close$A_chr,start=l1pa3.close$A_TSS-1,end=l1pa3.close$A_TSS,id=l1pa3.close$A_ID,score=".",strand=l1pa3.close$A_Strand),file = paste("output/",d,"/l1pa3.close.",d,".GeneOnly.bed",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  
  ## HERV
  herv.close.bed <- data.frame(chr=herv.close$A_chr,
                                tssLeft=ifelse(test=herv.close$B_TSS>herv.close$A_TSS,yes = herv.close$A_TSS,no = herv.close$B_TSS),
                                tssRight=ifelse(test=herv.close$B_TSS<herv.close$A_TSS,yes = herv.close$A_TSS,no = herv.close$B_TSS),
                                id=ifelse(test=herv.close$B_TSS<herv.close$A_TSS,yes = paste(herv.close$B_ID,herv.close$A_ID,sep = ":"),no = paste(herv.close$A_ID,herv.close$B_ID,sep = ":")),
                                score=".",
                                strandgene=ifelse(test=herv.close$B_TSS<herv.close$A_TSS,yes = "+",no = "-"))
  
  
  
  write.table(x = herv.close.bed,file = paste("output/",d,"/herv.close.",d,".bed",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(x = data.frame(chr=herv.close$A_chr,start=herv.close$B_TSS-1,end=herv.close$B_TSS,id=herv.close$B_ID,score=".",strand=herv.close$B_Strand),file = paste("output/",d,"/herv.close.",d,".TEonly.bed",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(x = data.frame(chr=herv.close$A_chr,start=herv.close$A_TSS-1,end=herv.close$A_TSS,id=herv.close$A_ID,score=".",strand=herv.close$A_Strand),file = paste("output/",d,"/herv.close.",d,".GeneOnly.bed",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  
  ## sva
  sva.close.bed <- data.frame(chr=sva.close$A_chr,
                                tssLeft=ifelse(test=sva.close$B_TSS>sva.close$A_TSS,yes = sva.close$A_TSS,no = sva.close$B_TSS),
                                tssRight=ifelse(test=sva.close$B_TSS<sva.close$A_TSS,yes = sva.close$A_TSS,no = sva.close$B_TSS),
                                id=ifelse(test=sva.close$B_TSS<sva.close$A_TSS,yes = paste(sva.close$B_ID,sva.close$A_ID,sep = ":"),no = paste(sva.close$A_ID,sva.close$B_ID,sep = ":")),
                                score=".",
                                strandgene=ifelse(test=sva.close$B_TSS<sva.close$A_TSS,yes = "+",no = "-"))
  
  
  
  write.table(x = sva.close.bed,file = paste("output/",d,"/sva.close.",d,".bed",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(x = data.frame(chr=sva.close$A_chr,start=sva.close$B_TSS-1,end=sva.close$B_TSS,id=sva.close$B_ID,score=".",strand=sva.close$B_Strand),file = paste("output/",d,"/sva.close.",d,".TEonly.bed",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(x = data.frame(chr=sva.close$A_chr,start=sva.close$A_TSS-1,end=sva.close$A_TSS,id=sva.close$A_ID,score=".",strand=sva.close$A_Strand),file = paste("output/",d,"/sva.close.",d,".GeneOnly.bed",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  
  
  ## Fuse with test
  l1hs.close.test <- merge(data.frame(l1hs.close),data.frame(dnmt$test$Default),by.x=1,by.y=0)
  l1pa2.close.test <- merge(data.frame(l1pa2.close),data.frame(dnmt$test$Default),by.x=1,by.y=0)
  l1pa3.close.test <- merge(data.frame(l1pa3.close),data.frame(dnmt$test$Default),by.x=1,by.y=0)
  herv.close.test <- merge(data.frame(herv.close),data.frame(dnmt$test$Default),by.x=1,by.y=0)
  sva.close.test <- merge(data.frame(sva.close),data.frame(dnmt$test$Default),by.x=1,by.y=0)
  
  ########
  
  ### Statistics ###
  lst <- list()
  t <- wilcox.test(l1hs.close.test$log2FoldChange,sva.close.test$log2FoldChange)
  lst <- list(lst,paste("L1HS - SVA: ",format(round(t$p.value,7)),sep = ""))
  t <- wilcox.test(l1hs.close.test$log2FoldChange,herv.close.test$log2FoldChange)
  lst <- list(lst,paste("L1HS - HERV: ",format(round(t$p.value,7)),sep = ""))
  
  t <- wilcox.test(l1pa2.close.test$log2FoldChange,sva.close.test$log2FoldChange)
  lst <- list(lst,paste("L1PA2 - SVA: ",format(round(t$p.value,7)),sep = ""))
  t <- wilcox.test(l1pa2.close.test$log2FoldChange,herv.close.test$log2FoldChange)
  lst <- list(lst,paste("L1PA2 - HERV: ",format(round(t$p.value,7)),sep = ""))
  
  t <- wilcox.test(l1pa3.close.test$log2FoldChange,sva.close.test$log2FoldChange)
  lst <- list(lst,paste("L1PA3 - SVA: ",format(round(t$p.value,7)),sep = ""))
  t <- wilcox.test(l1pa3.close.test$log2FoldChange,herv.close.test$log2FoldChange)
  lst <- list(lst,paste("L1PA3 - HERV: ",format(round(t$p.value,7)),sep = ""))
  
  print(lst,row.names=F)
  ### plot Mean #####
  plotCloseMean.ma <- function(close.df=NULL,class=NULL,col="yellow",cexbase=1.1) {
  
    plot(x = log2(dnmt$test$Default$baseMean),y = dnmt$test$Default$log2FoldChange,pch=16,cex=.2,ylab="log2(FC)",xlab="log2(baseMean)")
    abline(h = 0,lty=2,col="grey")
    points(x = log2(close.df$baseMean),y = close.df$log2FoldChange,pch=16,cex=cexbase+.2)
    points(x = log2(close.df$baseMean),y = close.df$log2FoldChange,pch=16,cex=cexbase,col=col)
    mtext(text = class)  
  }
  
  ## get same colors as boxplots/vioplots..
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  col <- gg_color_hue(5)
  
  pdf(file = paste("output/",d,"/maPlot_genes.",d/1000,"kb.L1HS.PA2.PA3_up.pdf",sep = ""),width = 10, height = 5)
  par(mfrow=c(1,3))
  plotCloseMean.ma(l1hs.close.test,"L1HS up",col=col[1],cexbase = 1)
  plotCloseMean.ma(l1pa2.close.test,"L1PA2 up",col=col[2],cexbase = 1)
  plotCloseMean.ma(l1pa3.close.test,"L1PA3 up",col=col[3],cexbase = 1)
  dev.off()
  
  pdf(file = paste("output/",d,"/maPlot_genes.",d/1000,"kb.L1HS.PA2.PA3.HERV.SVA_up.pdf",sep = ""),width = 15, height = 5)
  par(mfrow=c(1,5))
  plotCloseMean.ma(l1hs.close.test,"L1HS up",col=col[1])
  plotCloseMean.ma(l1pa2.close.test,"L1PA2 up",col=col[2])
  plotCloseMean.ma(l1pa3.close.test,"L1PA3 up",col=col[3])
  plotCloseMean.ma(herv.close.test,"HERV",col=col[4],cexbase = .5)
  plotCloseMean.ma(sva.close.test,"SVA",col=col[5],cexbase = .5)
  dev.off()
  
  plotCloseHeat <- function(close.df=NULL) {
    
    close.id <- close.df$A_ID  
    
    heatGenes(data = dnmt$VST,genes = close.id,z = T,cluster_col = F)
  }
  
  pdf(file = paste("output/",d,"/heatMap_genes.",d/1000,"kb_L1HSup.pdf",sep = ""))
  plotCloseHeat(close.df = l1hs.close)
  dev.off()
  pdf(file = paste("output/",d,"/heatMap_genes.",d/1000,"kb_L1PA2up.pdf",sep = ""))
  plotCloseHeat(close.df = l1pa2.close)
  dev.off()
  pdf(file = paste("output/",d,"/heatMap_genes.",d/1000,"kb_L1PA3up.pdf",sep = ""))
  plotCloseHeat(close.df = l1pa3.close)
  dev.off()
  pdf(file = paste("output/",d,"/heatMap_genes.",d/1000,"kb_HERV.pdf",sep = ""))
  plotCloseHeat(close.df = herv.close)
  dev.off()
  pdf(file = paste("output/",d,"/heatMap_genes.",d/1000,"kb_SVA.pdf",sep = ""))
  plotCloseHeat(close.df = sva.close)
  dev.off()
  
  
  
  ### FC ####
  l1hs.close.fc <- l1hs.close.test$log2FoldChange
  l1pa2.close.fc <- l1pa2.close.test$log2FoldChange
  l1pa3.close.fc <- l1pa3.close.test$log2FoldChange
  herv.close.fc <- herv.close.test$log2FoldChange
  sva.close.fc <- sva.close.test$log2FoldChange
  
  fc.df <- data.frame(fc=l1hs.close.fc,class="L1HS up")
  fc.df <- rbind(fc.df,data.frame(fc=l1pa2.close.fc,class="L1PA2 up"))
  fc.df <- rbind(fc.df,data.frame(fc=l1pa3.close.fc,class="L1PA3 up"))
  fc.df <- rbind(fc.df,data.frame(fc=herv.close.fc,class="HERV"))
  fc.df <- rbind(fc.df,data.frame(fc=sva.close.fc,class="SVA"))
  
  dev.off()
  
  pdf(file = paste("output/",d,"/FC_genes_",d/1000,"kb.violinplot.pdf",sep=""))
  ggplot(fc.df, aes(y = fc, x = class,fill=class)) +
    geom_violin(draw_quantiles = T) +
    theme(panel.grid.major = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    stat_summary(fun.y="median", geom="point") +
    ylab(label = "log2(FC DNMT1 KO / CTR)") +
    ggtitle(label = paste("FC genes <",d,"kb",sep = "")) +
    #  geom_jitter(width = .1,alpha=.05,cex=0.2) +
    geom_hline(yintercept = 0,linetype="dashed")
  dev.off()
  
  pdf(file = paste("output/",d,"/FC_genes_",d/1000,"kb.boxplot.pdf",sep = ""))
  ggplot(fc.df, aes(y = fc, fill = class,x=class)) +
    geom_boxplot(outlier.shape = NA) +
    theme(panel.grid.major = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ylab(label = "log2(FC DNMT1 KO / CTR)") +
    ggtitle(label = paste("FC genes <",d/1000,"kb",sep = "")) +
    geom_hline(yintercept = 0,linetype="dashed")
  dev.off()
  #  geom_jitter(width = .1,alpha=.1,cex=0.2)
  
  ########
  
  ####  Base means.. #####
  l1hs.base <- merge(l1hs.close$A_ID,dnmt$baseMean$Mean,by.x=1,by.y=0)
  l1pa2.base <- merge(l1pa2.close$A_ID,dnmt$baseMean$Mean,by.x=1,by.y=0)
  l1pa3.base <- merge(l1pa3.close$A_ID,dnmt$baseMean$Mean,by.x=1,by.y=0)
  herv.base <- merge(herv.close$A_ID,dnmt$baseMean$Mean,by.x=1,by.y=0)
  sva.base <- merge(sva.close$A_ID,dnmt$baseMean$Mean,by.x=1,by.y=0)
  
  base.df <- data.frame(Mean=l1hs.base$CTR,treat="CTR",class="L1HS up")
  base.df <- rbind(base.df,data.frame(Mean=l1hs.base$KO,treat="KO",class="L1HS up"))
  
  base.df <- rbind(base.df,data.frame(Mean=l1pa2.base$CTR,treat="CTR",class="L1PA2 up"))
  base.df <- rbind(base.df,data.frame(Mean=l1pa2.base$KO,treat="KO",class="L1PA2 up"))
  
  base.df <- rbind(base.df,data.frame(Mean=l1pa3.base$CTR,treat="CTR",class="L1PA3 up"))
  base.df <- rbind(base.df,data.frame(Mean=l1pa3.base$KO,treat="KO",class="L1PA3 up"))
  
  base.df <- rbind(base.df,data.frame(Mean=herv.base$CTR,treat="CTR",class="HERV up"))
  base.df <- rbind(base.df,data.frame(Mean=herv.base$KO,treat="KO",class="HERV up"))
  
  base.df <- rbind(base.df,data.frame(Mean=sva.base$CTR,treat="CTR",class="SVA up"))
  base.df <- rbind(base.df,data.frame(Mean=sva.base$KO,treat="KO",class="SVA up"))
  
  head(base.df)
  
  ### Base mean plots
  pdf(file = paste("output/",d,"/BaseMean_genes_",d/1000,"kb.boxplot.pdf",sep = ""))
  ggplot(data = base.df,aes(y = Mean,x=class,fill=treat)) +
    geom_boxplot(outlier.shape = NA) + 
    scale_y_continuous(limits = c(0,200)) +
    theme(panel.grid.major = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ylab(label = "baseMean") +
    ggtitle(label = paste("FC genes <",d/1000,"kb",sep = "")) +
    geom_hline(yintercept = 0,linetype="dashed")
  dev.off()  
  
  base.log.df <- base.df
  base.log.df$Mean <- log2(base.log.df$Mean+1)
  pdf(file = paste("output/",d,"/BaseMean.Log2_genes_",d/1000,"kb.boxplot.pdf",sep = ""))
  ggplot(data = base.log.df,aes(y = Mean,x=class,fill=treat)) +
    geom_boxplot(outlier.shape = NA) + 
    theme(panel.grid.major = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ylab(label = "log2(baseMean)") +
    ggtitle(label = paste("FC genes <",d/1000,"kb",sep = "")) +
    geom_hline(yintercept = 0,linetype="dashed")
  dev.off()
  
  ##### 
  

}

dev.off()
plot(x = 1,y = 1)
plotAll(d = 10000)
plotAll(d = 25000)
plotAll(d = 50000)
plotAll(d = 100000)
