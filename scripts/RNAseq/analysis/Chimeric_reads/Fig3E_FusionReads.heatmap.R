rm(list=ls())
#Import and parse data..
library(gplots)
setwd("~/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/RNAseq/PairedParam/LINE1 optimize/")

mRNA <- read.delim(file = "~/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/RNAseq/PairedParam/Results/NCBI/DESeq.KO.NCBI.txt",col.names = c("Row.names","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"),header=F,skip=1,stringsAsFactors = F)
TSS.TE <- read.delim(file = "~/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/RNAseq/PairedParam/Analyses/TSS analysis/TSS-overlap/data/input/intersected_gencode.v25.TSS.gene_type.proteinCoding.transcript_overlap_Retro.uniqID.txt",header=F,stringsAsFactors = F)
fam <- read.delim("~/Documents/bioinformatics/genomicData/hg38/Repeats/Retro.hg38.features.txt",header=F,stringsAsFactors = F)

head(TSS.TE)

TSS.TE[grep("GDPD4",TSS.TE$V10),]
names.uniq <- TSS.TE$V4
names <- gsub("\\_.*","",names.uniq)

TE.TSS.ID <- cbind(names.uniq,names,TSS.TE$V10)
TE.TSS.ID <- data.frame(TE.TSS.ID[!duplicated(TE.TSS.ID),])
TE.TSS.ID$pairUniq <- paste(TE.TSS.ID[,3],TE.TSS.ID[,1],sep=":")
TE.TSS.ID$pair <- paste(TE.TSS.ID[,3],TE.TSS.ID[,2],sep=":")

TSS.ID.diffex <- merge(TE.TSS.ID,mRNA,by.x=2,by.y=1)
TSS.diffex.fam <- merge(TSS.ID.diffex,fam,by.x=2,by.y=1)

fam.fc <- data.frame(Fam=TSS.diffex.fam[,10],Log2FC=TSS.diffex.fam[,4])
numbers.fam <- table(fam.fc$Fam)

# fusion
fusion.pa <- read.delim("~/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/RNAseq/PairedParam//Analyses/TSS analysis/Chimeric/Quant_Fusion/L1_NCBIgenes.Quant/L1PA.NCBIgenes.sameStrand.quant.ltr.Gene.txt",row.names = NULL)
fusion.hs <- read.delim("~/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/RNAseq/PairedParam//Analyses/TSS analysis/Chimeric/Quant_Fusion/L1_NCBIgenes.Quant/L1HS.NCBIgenes.sameStrand.quant.ltr.Gene.txt",row.names = NULL)
fusion <- rbind(fusion.hs,fusion.pa)
head(fusion)
fusion$pair <- paste(fusion$Gene,fusion$row.names,sep=":")

fusion.pair <- as.character(fusion$pair)
head(fusion.pair)

L1.TSS <- TE.TSS.ID[grep("L1PA2|L1HS|PA3",TE.TSS.ID[,1]),]

head(fusion)
head(L1.TSS)

l1.tss.fusion <- merge(L1.TSS,fusion,by.x=4,by.y=9)


plot.matrix <- data.frame(sapply(l1.tss.fusion[,-c(1:7)],as.numeric))
plot.matrix$name <- paste(l1.tss.fusion$names,l1.tss.fusion$V3,sep=" - ")


library(plyr)
plot.mat <- ddply(plot.matrix,.(name),summarize,
                  wt1=sum(s107),
                  wt2=sum(s108),
                  wt3=sum(s109),
                  ko1=sum(s110),
                  ko2=sum(s111),
                  ko3=sum(s112))

rnames <- plot.mat$name
rnames.uniq <- gsub("\\_.*","",rnames)
plot.mat <- as.matrix(plot.mat[,-1])
rownames(plot.mat) <- rnames
plot.mat <- plot.mat[order(rownames(plot.mat)),]
hmcols<-colorRampPalette(c("black","red","yellow","white"))(256)
hmcols<-colorRampPalette(c("blue","white","red"))(256)
heatmap.2(plot.mat,col=hmcols,scale="row",
          margin=c(10,8),trace='none',Rowv = F,cexRow = 0.7)
