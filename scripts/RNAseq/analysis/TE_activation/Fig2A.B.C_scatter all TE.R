######## L1 sense vs sensesense qusensefication 
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

rm(list=ls())
gc()

setwd(dir = "~/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/RNAseq/PairedParam/")

sense <- read.delim(file = "Quant/hg38.unique.ERE.notInNCBI.Exon.txt",header=T,skip=1)
header <- read.delim(file="Quant/hg38.unique.ERE.notInNCBI.Exon.txt",header=F,skip=1,nrows=1,stringsAsFactors = F)
header <- gsub(pattern = "/projects/bmc-scc/user/medpvb/backup/projects/DNMT1KO/Aligned_hg38_STAR_unique/",replacement = "",header)
header <- gsub(pattern = "Aligned.sortedByCoord.out.bam",replacement = "",header)

### Remove first two rows, they are read wrong
colnames(sense) <- header
head(sense)

library(devtools)
install_github("perllb/deseqAbstraction")
library(deseqAbstraction)

# set path and filenames
path <- "/home/pbrattaas/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/RNAseq/PairedParam/Quant/"
te.file <- "hg38.unique.ERE.notInNCBI.Exon.txt"
ncbi.file <- "hg38.NCBI.genesBestRefSeq.exon.primary.s2.txt"

## structure sense data
pos <- sense[,1:5]
head(pos)

length <- sense[,6]
rnames <- sense[,1]

# filter out non-expressed EREs
mat.data <- sense[,7:ncol(sense)]
rownames(mat.data) <- rnames
head(mat.data)
data.filter <- mat.data[rowSums(mat.data)>5,]
head(data.filter)

coln <- c("hNES_CTR_1","hNES_CTR_2","hNES_CTR_3","hNES_KO_1","hNES_KO_2","hNES_KO_3")

colnames(data.filter) <- coln

countData.s <- data.filter
head(countData.s)

condition <- c(rep("CTR",3),rep("KO",3))
colData <- data.frame(condition=condition)
rownames(colData) <- colnames(data.filter)
colnames(countData.s) <- rownames(colData)
rownames(countData.s) <- rownames(data.filter)
head(countData.s)
colData

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = countData.s,
                              colData = colData,
                              design = ~ condition)

## Scale to sequencing depth, using NCBI gene counts
summary <- read.delim("Quant/hg38.NCBI.genesBestRefSeq.exon.primary.s2.txt.summary")

## Get columns with FWD reads only
fwdCol <- grep("FWD",colnames(summary))
revCol <- grep("REV",colnames(summary))

## remove those columns 
remove <- c(fwdCol,revCol)
sum.data <- summary[,-remove]
head(sum.data)

## Scale to read depth, using NCBI gene counts 
read.genome <- colSums(sum.data[c(1,2,4),-1])
sizeFactors(dds) <- read.genome/mean(read.genome)
dds <- DESeq(dds)

# log transform, norm
vsd <- varianceStabilizingTransformation(dds,blind=F)
head(assay(vsd))

# results diffex
res <- results(dds)
rnames.res <- gsub(pattern = "\\_.*","",rownames(res))
head(res)

### plot RetroFamilies

# read repeatmasker family data
#fam <- read.delim(file="~/genomicData/Repeats/hg38/Retro.hg38.features.uniqID.txt",header=F)
fam <- read.delim(file="~/Documents/bioinformatics/genomicData/hg38/Repeats/Retro.hg38.features.uniqID.txt",header=F)

# merge diffex with repeatmasker class/fam
res.fam <- merge(rownames(res),fam,by.x=1,by.y=1,sort=F)
res
head(res.fam)
dim(res.fam)
dim(res)
tail(res.fam)

## color
col <- ifelse(res.fam$V2=="LINE",yes="darkolivegreen",
              no=ifelse(res.fam$V2=="LTR",yes="darkblue",
                        no = ifelse(res.fam$V3=="SVA",yes="red",
                                    no=ifelse(res.fam$V2=="SINE",yes="orange",no="black"))))
col.f <- ifelse(res$padj>0.05,yes = "black",no=col)
#cbind(res.fam,col)
cex <- ifelse(res$padj>0.05,0.1,0.6)
data.frame(res[!is.na(res$padj) & res$padj<0.05,])

# ma plot
plot(x=log2(res$baseMean),y = res$log2FoldChange,col=col.f,pch=16,ce=cex,main="DNMT1KO: Retroelement response",xlim=c(.5,9),ylab="log2(FC)",xlab="log2(mean expression)")
legend("topleft",legend = c("LINE","LTR","SVA","SINE"),col = c("darkolivegreen","darkblue","red","orange"),pch=16,cex = 0.7)

# volcano plot
plot(x=res$log2FoldChange,y = log2(res$baseMean),col=col.f,pch=16,ce=cex,xlim=c(-11,11),main="DNMT1KO: Retroelement response")
legend("topleft",legend = c("LINE","LTR","SVA","SINE"),col = c("darkolivegreen","darkblue","red","orange"),pch=16,cex = 0.7)

######

## stats on changed TEs 
up <- res[!is.na(res$padj) & res$padj<0.05 & res$log2FoldChange>0,]
up.fam <- merge(rownames(up),fam,by.x=1,by.y=1,sort=F)

up.fam.table <- table(up.fam$V3)
up.fam.table[up.fam.table>30]

up.class.table <- table(up.fam$V2)
up.class.table[up.class.table>3]


down <- res[!is.na(res$padj) & res$padj<0.05 & res$log2FoldChange<0,]
down.fam <- merge(rownames(down),fam,by.x=1,by.y=1,sort=F)

down.fam.table <- table(down.fam$V3)
down.fam.table[down.fam.table>30]

down.class.table <- table(down.fam$V2)
down.class.table[down.class.table>0]

unchange <- res[!is.na(res$padj) & res$padj>=0.05,]
unchange.fam <- merge(rownames(unchange),fam,by.x=1,by.y=1,sort=F)

unchange.fam.table <- table(unchange.fam$V3)
unchange.fam.table[unchange.fam.table>30]

unchange.class.table <- table(unchange.fam$V2)
unchange.class.table[unchange.class.table>3]

sum(unchange.class.table[unchange.class.table>3],up.class.table[up.class.table>3],down.class.table[down.class.table>3])

par(mfrow=c(1,1))
pie(up.class.table[up.class.table>3],labels = paste(names(up.class.table[up.class.table>3]),up.class.table[up.class.table>3],sep = ": "),main="Sign up")
pie(down.class.table[down.class.table>3],labels = paste(names(down.class.table[down.class.table>3]),down.class.table[down.class.table>3],sep = ": "),main="Sign down")
pie(unchange.class.table[unchange.class.table>3],labels = paste(names(unchange.class.table[unchange.class.table>3]),unchange.class.table[unchange.class.table>3],sep = ": "),main="unchanged")

#####

# Barplot of number of upregulated TEs in subfams

tab <- table(up.fam.all)
tab.ord <- tab[order(-tab)]

# Read repeatmasker class
fam.nonID <- read.delim(file="~/Documents/bioinformatics/genomicData/hg38/Repeats/Retro.hg38.features.txt",header=F)

# barplot 
tab.fam <- merge(names(tab.ord),fam.nonID,by=1,sort=F)
col <- ifelse(tab.fam$V2=="LINE",yes="darkolivegreen",
              no=ifelse(tab.fam$V2=="LTR",yes="darkblue",
                        no = ifelse(tab.fam$V3=="SVA",yes="red",
                                    no=ifelse(tab.fam$V2=="SINE",yes="orange",no="black"))))
x <- barplot(tab.ord[1:20],las=2,col=col,ylab="number of upregulated elements",main="DNMT1-KO: upregulated elements by family")


## plot only HS-L1PA - sorted by age
l1 <- tab.ord[grep("L1HS|L1PA",names(tab.ord))]
l1 <- l1[order(names(l1))]
cbind(l1,1:length(l1))
l1.ord <- l1[c(1,11:18,2:10)]

numbers <- read.delim("~/Documents/bioinformatics/genomicData/hg38/Repeats/hg38.RM.Numbers.txt",sep=" ")
l1.num <- numbers[grep("L1HS|L1PA",numbers$X.A.n),]

l1.merge <- merge(data.frame(l1.ord),l1.num,by.x=1,by.y=2)
l1.merge <- l1.merge[c(1,11:18,2:10),]
rownames(l1.merge) <- make.names(l1.merge$up.fam.all,unique = T)
x <- barplot(l1.merge$Freq,col="darkolivegreen")
lines(x,l1.merge$X27189/40,lty=2)
axis(side = 4,at = seq(0,400,100),labels = seq(0,400,100)*40)
mtext(text = "# upregulated",side = 2,line=2)
mtext(text = "# total in family",side=4,line=2)
axis(side=1,at = x,labels = l1.merge$up.fam.all,las=2)
