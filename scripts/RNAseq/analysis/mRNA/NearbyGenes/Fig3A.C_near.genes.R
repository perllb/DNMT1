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


# Get close genes (50kb)
d <- 50000
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



## L1PA2
l1pa2.close.bed <- data.frame(chr=l1pa2.close$A_chr,
                             tssLeft=ifelse(test=l1pa2.close$B_TSS>l1pa2.close$A_TSS,yes = l1pa2.close$A_TSS,no = l1pa2.close$B_TSS),
                             tssRight=ifelse(test=l1pa2.close$B_TSS<l1pa2.close$A_TSS,yes = l1pa2.close$A_TSS,no = l1pa2.close$B_TSS),
                             id=ifelse(test=l1pa2.close$B_TSS<l1pa2.close$A_TSS,yes = paste(l1pa2.close$B_ID,l1pa2.close$A_ID,sep = ":"),no = paste(l1pa2.close$A_ID,l1pa2.close$B_ID,sep = ":")),
                             score=".",
                             strandgene=ifelse(test=l1pa2.close$B_TSS<l1pa2.close$A_TSS,yes = "+",no = "-"))



## L1PA3
l1pa3.close.bed <- data.frame(chr=l1pa3.close$A_chr,
                             tssLeft=ifelse(test=l1pa3.close$B_TSS>l1pa3.close$A_TSS,yes = l1pa3.close$A_TSS,no = l1pa3.close$B_TSS),
                             tssRight=ifelse(test=l1pa3.close$B_TSS<l1pa3.close$A_TSS,yes = l1pa3.close$A_TSS,no = l1pa3.close$B_TSS),
                             id=ifelse(test=l1pa3.close$B_TSS<l1pa3.close$A_TSS,yes = paste(l1pa3.close$B_ID,l1pa3.close$A_ID,sep = ":"),no = paste(l1pa3.close$A_ID,l1pa3.close$B_ID,sep = ":")),
                             score=".",
                             strandgene=ifelse(test=l1pa3.close$B_TSS<l1pa3.close$A_TSS,yes = "+",no = "-"))


## HERV
herv.close.bed <- data.frame(chr=herv.close$A_chr,
                              tssLeft=ifelse(test=herv.close$B_TSS>herv.close$A_TSS,yes = herv.close$A_TSS,no = herv.close$B_TSS),
                              tssRight=ifelse(test=herv.close$B_TSS<herv.close$A_TSS,yes = herv.close$A_TSS,no = herv.close$B_TSS),
                              id=ifelse(test=herv.close$B_TSS<herv.close$A_TSS,yes = paste(herv.close$B_ID,herv.close$A_ID,sep = ":"),no = paste(herv.close$A_ID,herv.close$B_ID,sep = ":")),
                              score=".",
                              strandgene=ifelse(test=herv.close$B_TSS<herv.close$A_TSS,yes = "+",no = "-"))



## sva
sva.close.bed <- data.frame(chr=sva.close$A_chr,
                              tssLeft=ifelse(test=sva.close$B_TSS>sva.close$A_TSS,yes = sva.close$A_TSS,no = sva.close$B_TSS),
                              tssRight=ifelse(test=sva.close$B_TSS<sva.close$A_TSS,yes = sva.close$A_TSS,no = sva.close$B_TSS),
                              id=ifelse(test=sva.close$B_TSS<sva.close$A_TSS,yes = paste(sva.close$B_ID,sva.close$A_ID,sep = ":"),no = paste(sva.close$A_ID,sva.close$B_ID,sep = ":")),
                              score=".",
                              strandgene=ifelse(test=sva.close$B_TSS<sva.close$A_TSS,yes = "+",no = "-"))


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

par(mfrow=c(1,3))
plotCloseMean.ma(l1hs.close.test,"L1HS up",col=col[1],cexbase = 1)
plotCloseMean.ma(l1pa2.close.test,"L1PA2 up",col=col[2],cexbase = 1)
plotCloseMean.ma(l1pa3.close.test,"L1PA3 up",col=col[3],cexbase = 1)

par(mfrow=c(2,3))
plotCloseMean.ma(l1hs.close.test,"L1HS up",col=col[1])
plotCloseMean.ma(l1pa2.close.test,"L1PA2 up",col=col[2])
plotCloseMean.ma(l1pa3.close.test,"L1PA3 up",col=col[3])
plotCloseMean.ma(herv.close.test,"HERV",col=col[4],cexbase = .5)
plotCloseMean.ma(sva.close.test,"SVA",col=col[5],cexbase = .5)

plotCloseHeat <- function(close.df=NULL) {
  
  close.id <- close.df$A_ID  
  
  heatGenes(data = dnmt$VST,genes = close.id,z = T,cluster_col = F)
}

plotCloseHeat(close.df = l1hs.close)

plotCloseHeat(close.df = l1pa2.close)

plotCloseHeat(close.df = l1pa3.close)
 
plotCloseHeat(close.df = herv.close)

plotCloseHeat(close.df = sva.close)

### FC violin plot ####
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

ggplot(fc.df, aes(y = fc, x = class,fill=class)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="median", geom="point") +
  ylab(label = "log2(FC DNMT1 KO / CTR)") +
  ggtitle(label = paste("FC genes <",d,"kb",sep = "")) +
  #  geom_jitter(width = .1,alpha=.05,cex=0.2) +
  geom_hline(yintercept = 0,linetype="dashed")


# boxplot
ggplot(fc.df, aes(y = fc, fill = class,x=class)) +
  geom_boxplot(outlier.shape = NA) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab(label = "log2(FC DNMT1 KO / CTR)") +
  ggtitle(label = paste("FC genes <",d/1000,"kb",sep = "")) +
  geom_hline(yintercept = 0,linetype="dashed")

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
ggplot(data = base.df,aes(y = Mean,x=class,fill=treat)) +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(limits = c(0,200)) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab(label = "baseMean") +
  ggtitle(label = paste("baseMean genes <",d/1000,"kb",sep = "")) +
  geom_hline(yintercept = 0,linetype="dashed")

base.log.df <- base.df
base.log.df$Mean <- log2(base.log.df$Mean+1)

ggplot(data = base.log.df,aes(y = Mean,x=class,fill=treat)) +
  geom_boxplot(outlier.shape = NA) + 
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab(label = "log2(baseMean)") +
  ggtitle(label = paste("log2 basemean genes <",d/1000,"kb",sep = "")) +
  geom_hline(yintercept = 0,linetype="dashed")

##### 

