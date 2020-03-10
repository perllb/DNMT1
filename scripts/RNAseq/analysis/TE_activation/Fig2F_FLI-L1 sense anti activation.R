###### check FLI L1
rm(list=ls())

setwd("~/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/RNAseq/PairedParam/Results/LINE/")

sense <- read.delim(file = "ERE.DNMT1KO.BaseMean_sense.pos.txt",row.names = 1)[,-1]
head(sense)

anti <- read.delim(file = "ERE.DNMT1KO.BaseMean_anti.pos.txt",row.names = 1)[,-1]
head(anti)


sense.dat <- sense[,5:6]
anti.dat <- anti[,5:6]
sense.l1 <- sense.dat[grep("L1HS|L1PA2",rownames(sense.dat)),]
anti.l1 <- anti.dat[grep("L1HS|L1PA2",rownames(anti.dat)),]

merge <- merge(sense.l1,anti.l1,by=0,all=T)
head(merge)
colnames(merge) <- c("L1","SENSE\nCTR","SENSE\nKO","ANTI\nCTR","ANTI\nKO")

#FLI L1 
fli <- read.delim("~/Documents/bioinformatics/genomicData/hg38/Repeats/l1base2.intersect.L1RM_ID.txt",header=F)
head(fli)

merge.fli <- merge(fli,merge,by=1)
merge.fli[is.na(merge.fli)]<-0
head(merge.fli)

boxplot(merge.fli[,-1],cex=0.8,pch=16,xaxt='n',ylab="Read number")
axis(side = 1,at = 1:4,labels = colnames(merge.fli[-1]),tick = F)
mtext(text = "DNMT1-KO: Response on FLI-L1",side = 3,line = 1,at = 1.3,font=2)

wilcox.test(merge.fli$`SENSE
CTR`,merge.fli$`SENSE
KO`,correct = F)
