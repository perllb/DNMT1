rm(list=ls())
gc()
#########################
## Clustering TRIM28 KD
########################
source("http://bioconductor.org/biocLite.R")
#biocLite("vsn")
#biocLite("DESeq2")
library(DESeq2)
library("vsn")
library("RColorBrewer")
#install.packages("gplots")
library("gplots")
#install.packages("devtools")
library(devtools)
#install_github("ggbiplot","vqv")
#library(ggbiplot)

#setwd("/home/pbrattaas//Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/RNAseq/PairedParam/Quant/")
setwd("/Users/perludvik/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/RNAseq/PairedParam/Quant/")

##################################################
# Read data and transform to matrix
##################################################

data <- read.csv("hg38.NCBI.genesBestRefSeq.exon.primary.s2.txt",header=T,sep = "\t",skip=1)
head(data)

pos <- data[,1:5]
head(pos)

length <- data[,6]
rnames <- data[,1]
mat.data <- data.matrix(data[,7:ncol(data)])
mat.data.2 <- mat.data
head(mat.data.2)
coln <- c("hNES_CTR_1","hNES_CTR_2","hNES_CTR_3","hNES_KO_1","hNES_KO_2","hNES_KO_3")
colnames(mat.data.2) <- coln
rownames(mat.data.2) <- rnames


countData <- mat.data.2
head(countData)
countData[grep("DNMT1",rownames(countData)),]

condition <- c(rep("CTR",3),rep("KO",3))
colData <- data.frame(condition=condition)
rownames(colData) <- colnames(mat.data.2)
colnames(countData) <- rownames(colData)
rownames(countData) <- rnames
head(countData)

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)

dds <- DESeq(dds)

baseMeanPerLvl <- sapply( levels(dds$condition), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$condition == lvl] ) )
baseSEPerLvl <- sapply( levels(dds$condition), function(lvl) apply( counts(dds,normalized=TRUE)[,dds$condition == lvl],1,sd ) )*1.96/sqrt(3)
head(baseMeanPerLvl)
head(baseSEPerLvl)

baseMeanPerLvl[grep("DNMT1",rownames(baseMeanPerLvl)),]

# log transform, norm
vsd <- varianceStabilizingTransformation(dds,blind=F)
head(assay(vsd))

PCAplotter(dat = vsd,title = "NCBI top 5000",ntop = 5000,color = colData$condition,shape=condition,label = rownames(colData))

##### PCA 
data <- plotPCA(vsd,intgroup = c("condition"),returnData = T)
percentVar <- round(100*attr(data,"percentVar"))
library("ggplot2")
ggplot(data,aes(PC1,PC2,color=condition)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=data$name),size=3)

plot(data$PC1,data$PC2,pch=16,col=data$group,
     xlab=paste0("PC1: ",percentVar[1],"% variance"),
     ylab=paste0("PC2: ",percentVar[2],"% variance"))
#identify(data$PC1,data$PC2,labels = data$name)

head(assay(vsd))
head(counts(dds,normalized=T))

###### Diffex
colData
resSai <- results(dds)
resSaiSort <- resSai[order(resSai$pvalue),]
dim(data.frame(resSaiSort[!is.na(resSaiSort$padj) && resSaiSort$padj<0.05 && resSaiSort$log2FoldChange>0,]))
dim(resSaiSort)

###  tables sorted on log2fc
resSai.FCorder <- resSai[order(-resSai$log2FoldChange),]

up <- resSai[!is.na(resSai$padj) & resSai$padj<0.05 & resSai$log2FoldChange>0,]
down <- resSai[!is.na(resSai$padj) & resSai$padj<0.05 & resSai$log2FoldChange<0,]
notSign <- resSai[!is.na(resSai$padj) & resSai$padj>0.05,]

col <- ifelse(resSai$padj<0.001,yes=ifelse(resSai$log2FoldChange>0,yes="firebrick3",no="steelblue4"),"black")

## FC basemean
cex <- ifelse(resSai$padj<0.001,0.4,0.1)

# Ma plot
plot(log(resSai$baseMean),resSai$log2FoldChange,pch=16,cex=.1)







