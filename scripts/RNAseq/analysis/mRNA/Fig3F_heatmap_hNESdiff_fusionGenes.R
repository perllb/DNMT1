library(gplots)

rm(list=ls())
setwd("~/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/RNAseq/PairedParam/Results/")
mRNA <- read.delim(file = "~/LocalData/Falk/Quant_Falk/hg38.s2.unique.NCBI.Exon.txt",header=T,skip=1,stringsAsFactors = F)
TSS.TE <- read.delim(file = "~/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/RNAseq/PairedParam/Analyses/TSS analysis/TSS-overlap/data/input/intersected_gencode.v25.TSS.gene_type.proteinCoding.transcript_overlap_Retro.uniqID.txt",header=F,stringsAsFactors = F)
fam <- read.delim("~/genomicData/Repeats/hg38/Retro.hg38.features.txt",header=F,stringsAsFactors = F)

head(TSS.TE)
TSS.TE[grep("GDPD4",TSS.TE$V10),]
names.uniq <- TSS.TE$V4
names <- gsub("\\_.*","",names.uniq)

TE.TSS.ID <- cbind(names,TSS.TE$V10)
TE.TSS.ID <- TE.TSS.ID[!duplicated(TE.TSS.ID),]

#mRNA <- read.delim(file = "/Users/perludvik/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/LINE1 Manuscript Figures - scripts and plots/hNES diff Falk RNAseq/hg38.s2.unique.NCBI.Exon.txt",header=T,skip=1,stringsAsFactors = F)
#TSS.TE <- read.delim(file = "/Users/experludvik//Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/RNAseq/PairedParam/Analyses/TSS analysis/TSS-overlap/data/input/intersected_gencode.v25.TSS.gene_type.proteinCoding.transcript_overlap_Retro.uniqID.txt",header=F,stringsAsFactors = F)
#fam <- read.delim("/Users/perludvik/Documents/bioinformatics/genomicData/hg38//Retro.hg38.features.txt",header=F,stringsAsFactors = F)

## structure mRNA data
head(mRNA)
mRNA.count <- mRNA[,7:ncol(mRNA)]
head(mRNA.count)
cnames <- colnames(mRNA.count)
cnames <- gsub(pattern = "X.projects.fs1.medpvb.backup.projects.Falk.hNES.diff.Aligned_hg38_STAR_unique.hg38.unique.",replacement = "",cnames)
cnames <- gsub(pattern = "Aligned.sortedByCoord.out.bam",replacement = "",cnames)
colnames(mRNA.count) <- cnames
rownames(mRNA.count) <- mRNA$Geneid
head(mRNA.count)
order <- c(5,10,2,7,3,8,4,9,1,6)
mRNA.count <- mRNA.count[,order]
head(mRNA.count)

library(DESeq2)

rep <- rep(c("r1","r2"),5)
d <- c("d0","d0","d25","d25","d50","d50","d75","d75","d100","d100")
colData <- data.frame(condition=d,rep=rep)
rownames(colData) <- colnames(mRNA.count)

dds <- DESeqDataSetFromMatrix(countData = mRNA.count,colData = colData,design =~ condition)
dds <- DESeq(dds)

vsd <- varianceStabilizingTransformation(dds)
assay <- assay(vsd)
head(assay)

fusion.pa <- read.delim("../Analyses/TSS analysis/Chimeric/Quant_Fusion/L1_NCBIgenes.Quant/L1PA.NCBIgenes.sameStrand.quant.ltr.Gene.txt",row.names = NULL)
fusion.hs <- read.delim("../Analyses/TSS analysis/Chimeric/Quant_Fusion/L1_NCBIgenes.Quant/L1HS.NCBIgenes.sameStrand.quant.ltr.Gene.txt",row.names = NULL)
fusion <- rbind(fusion.hs,fusion.pa)
head(fusion)

L1.TSS <- TE.TSS.ID[grep("L1PA2|L1HS|PA3|PA4",TE.TSS.ID[,1]),]
l1.tss.fusion <- merge(L1.TSS,fusion[,-1],by.x=2,by.y=1)
genes.tss.fusion <- as.character(l1.tss.fusion$V2)

dnmt1 <- read.delim("../Results/DNMT1KO_DESeq.NCBI.txt");
head(dnmt1)
dnmt1.up <- dnmt1[!is.na(dnmt1$padj) & dnmt1$log2FoldChange>0 & dnmt1$padj<0.05,]
genes.up <- as.character(dnmt1.up$Row.names)
head(genes.up)

inter <- intersect(genes.tss.fusion,genes.up)

diff.fusion <- setdiff(genes.tss.fusion,genes.up)
diff.up <- setdiff(genes.up,genes.tss.fusion)

head(assay)

plot.assay <- merge(inter,assay,by.x=1,by.y=0)
head(plot.assay)

plot.mat <- as.matrix(plot.assay[,-1])
rownames(plot.mat) <- plot.assay[,1]
hmcols<-colorRampPalette(c("darkblue","white","red","darkred"))(256)
heatmap.2(plot.mat,col=hmcols,scale="row",Colv=F,
          margin=c(10,8),trace='none',Rowv = F,cexRow = 0.7)


## plot means
plot.mean <- data.frame(nes=rowMeans(plot.mat[,1:2]),
                        d25=rowMeans(plot.mat[,3:4]),
                        d50=rowMeans(plot.mat[,5:6]),
                        d75=rowMeans(plot.mat[,7:8]),
                        d100=rowMeans(plot.mat[,9:10]))
plot.mean

heatmap.2(as.matrix(plot.mean),col=hmcols,scale="row",Colv=F,
          margin=c(10,8),trace='none',Rowv = F,cexRow = 0.7)
