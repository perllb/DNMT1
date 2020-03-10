## Violin plot coverage
library(ggplot2)
library(tidyverse)

setwd(dir = "~/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/2018.08.30_DNMT1L1_figures_scripts/WGBS/analysis/compFeatures/")
###### 3x coverage ######

##### Classes ######

dat <- data.frame()

for(class in c("LINE","SINE","LTR","SVA")){
  
  print(class)
  wt <- paste("3x/outputTE/beds/WT_",class,".mean.3x_rmNAN.bed",sep="")
  ko <- paste("3x/outputTE/beds/DNMT1ko_",class,".mean.3x_rmNAN.bed",sep = "")
  wt <- read.delim(file = wt,header=F)
  ko <- read.delim(file = ko,header=F)
  wt.df <- data.frame(value=wt,variable="WT",class=class)
  ko.df <- data.frame(value=ko,variable="KO",class=class)
  dat <- rbind(dat,wt.df,ko.df)

}

head(dat)
dat %>% group_by(class,variable) %>%
  summarize(mean=mean(V1))
  

  

pdf(file = "3x/outputTE/plots/Classes.3x.meanOfTEFragments.boxplot.pdf")
ggplot(dat, aes(y = V1, x = variable,fill=class)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab(label = "% mCpG / CpG") +
  ggtitle(label = "TE class levels",subtitle = "3x coverage") +
  geom_boxplot(width=0.05,outlier.size = .5) +
  facet_wrap(~class, scale="free",nrow = 1) +
  scale_fill_brewer(palette="Pastel1")
dev.off()

pdf(file = "3x/outputTE/plots/Classes.3x.meanOfTEFragments.pdf")
ggplot(dat, aes(y = V1, x = variable,fill=class)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab(label = "% mCpG / CpG") +
  ggtitle(label = "TE class levels",subtitle = "3x coverage") +
  stat_summary(fun.y="mean", geom="point") +
  facet_wrap(~class, scale="free",nrow = 1) +
  scale_fill_brewer(palette="Pastel1")
dev.off()



##### L1HS FL6000 promoter #####
title <- "L1HS.FL6000.promoter"
wt <- "3x/L1/WT_L1HS.FL6000.promoter_mean.3x_rmNAN.bed"
ko <- "3x/L1/DNMT1ko_L1HS.FL6000.promoter_mean.3x_rmNAN.bed"
wt <- read.delim(file = wt,header=F)
ko <- read.delim(file = ko,header=F)
wt.df <- data.frame(value=wt,variable="WT")
ko.df <- data.frame(value=ko,variable="KO")
head(wt.df)
head(ko.df)

dat <- rbind(wt.df,ko.df)
head(dat)

pdf(file = paste("3x/outputTE/plots/",title,"_mean.3x_boxplot.pdf",sep = ""))
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "% mCpG / CpG") +
  ggtitle(label = paste(title," levels",sep = ""),subtitle = "3x coverage") +
  geom_boxplot(width=0.02,outlier.size = .5)
dev.off()

pdf(file = paste("3x/outputTE/plots/",title,"_mean.3x.pdf",sep = ""))
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "% mCpG / CpG") +
  ggtitle(label = paste(title," levels",sep = ""),subtitle = "3x coverage")
dev.off()

pdf(file = paste("3x/outputTE/plots/",title,"_mean.3x.jitter.pdf",sep = ""))
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "% mCpG / CpG") +
  ggtitle(label = paste(title," levels",sep = ""),subtitle = "3x coverage") +
  geom_jitter(width = .1,alpha=.1)
dev.off()

#### Counts ..
title <- "L1HS.FL6000.promoter"
wt <- "3x/L1/WT_L1HS.FL6000.promoter.3x_count.bed"
ko <- "3x/L1/DNMT1ko_L1HS.FL6000.promoter.3x_count.bed"
wt <- read.delim(file = wt,header=F)
ko <- read.delim(file = ko,header=F)
wt.df <- data.frame(value=wt,variable="WT")
ko.df <- data.frame(value=ko,variable="KO")
head(wt.df)
head(ko.df)

dat <- rbind(wt.df,ko.df)
head(dat)

# jitter
pdf(file = paste("3x/outputTE/plots/",title,"_COUNTS.3x_jitter.pdf",sep = ""))
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "Number of CpGs included per element") +
  ggtitle(label = paste(title," CpGs included",sep = ""),subtitle = "Number of CpG included per element (3x coverage)") +
  geom_jitter(width = 0.2,alpha=0.1)
dev.off()

# no jitter
pdf(file = paste("3x/outputTE/plots/",title,"_COUNTS.3x.pdf",sep = ""))
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "Number of CpGs included per element") +
  ggtitle(label = paste(title," CpGs included",sep = ""),subtitle = "Number of CpG included per element (3x coverage)")
dev.off()

pdf(file = paste("3x/outputTE/plots/",title,"_COUNTS.3x.hist.pdf",sep = ""))
ggplot(dat, aes(fill = variable, x = V1)) +
  geom_histogram(binwidth = 1,alpha=.5,position="identity") +
  ggtitle(label = "Number of CpG called",subtitle = paste("each",title,"(with at least 3x coverage)")) +
  xlab(label = "number of CpGs with at least 3x coverage")
dev.off()

###########

##### L1HS FL6000 full fragment #####
title <- "L1HS.FL6000"
ko <- "3x/outputTE/beds/DNMT1ko_L1HS.FL6000mean.3x_rmNAN.bed"
wt <- "3x/outputTE/beds/WT_L1HS.FL6000mean.3x_rmNAN.bed"

wt <- read.delim(file = wt,header=F)
ko <- read.delim(file = ko,header=F)

wt.df <- data.frame(value=wt,variable="WT")
ko.df <- data.frame(value=ko,variable="KO")
head(wt.df)
head(ko.df)

dat <- rbind(wt.df,ko.df)
head(dat)

pdf(file = paste("3x/outputTE/plots/",title,"_mean.3x_boxplot.pdf",sep = ""))
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "% mCpG / CpG") +
  ggtitle(label = paste(title," levels",sep = ""),subtitle = "3x coverage") +
  geom_boxplot(width=0.02,outlier.size = .5)
dev.off()

pdf(file = paste("3x/outputTE/plots/",title,"_mean.3x.pdf",sep = ""))
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "% mCpG / CpG") +
  ggtitle(label = paste(title," levels",sep = ""),subtitle = "3x coverage")
dev.off()

pdf(file = paste("3x/outputTE/plots/",title,"_mean.3x.jitter.pdf",sep = ""))
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "% mCpG / CpG") +
  ggtitle(label = paste(title," levels",sep = ""),subtitle = "3x coverage") +
  geom_jitter(width = .1,alpha=.1)
dev.off()


#### Counts ..
title <- "L1HS.FL6000"
ko <- "3x/outputTE/beds/DNMT1ko_L1HS.FL6000mean.3x_count.bed"
wt <- "3x/outputTE/beds/WT_L1HS.FL6000mean.3x_count.bed"
wt <- read.delim(file = wt,header=F)
ko <- read.delim(file = ko,header=F)
wt.df <- data.frame(value=wt,variable="WT")
ko.df <- data.frame(value=ko,variable="KO")
head(wt.df)
head(ko.df)

dat <- rbind(wt.df,ko.df)
head(dat)

# jitter
pdf(file = paste("3x/outputTE/plots/",title,"_COUNTS.3x_jitter.pdf",sep = ""))
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "Number of CpGs included per element") +
  ggtitle(label = paste(title," CpGs included",sep = ""),subtitle = "Number of CpG included per element (3x coverage)") +
  geom_jitter(width = 0.2,alpha=0.1)
dev.off()

# no jitter
pdf(file = paste("3x/outputTE/plots/",title,"_COUNTS.3x.pdf",sep = ""))
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "Number of CpGs included per element") +
  ggtitle(label = paste(title," CpGs included",sep = ""),subtitle = "Number of CpG included per element (3x coverage)")
dev.off()

pdf(file = paste("3x/outputTE/plots/",title,"_COUNTS.3x.hist.pdf",sep = ""))
ggplot(dat, aes(fill = variable, x = V1)) +
  geom_histogram(binwidth = 1,alpha=.5,position="identity") +
  ggtitle(label = "Number of CpG called",subtitle = paste("each",title,"(with at least 3x coverage)")) +
  xlab(label = "number of CpGs with at least 3x coverage")
dev.off()

###########


##### FLI L1 full fragment #####
## 3x
title <- "FLI-L1"
ko <- "3x/outputTE/beds/DNMT1ko_FLI-L1.mean.3x_rmNAN.bed"
wt <- "3x/outputTE/beds/WT_FLI-L1.mean.3x_rmNAN.bed"

wt <- read.delim(file = wt,header=F)
ko <- read.delim(file = ko,header=F)

wt.df <- data.frame(value=wt,variable="WT")
ko.df <- data.frame(value=ko,variable="KO")
head(wt.df)
head(ko.df)

mean(wt.df$V1)
mean(ko.df$V1)

dat <- rbind(wt.df,ko.df)
head(dat)

pdf(file = paste("3x/outputTE/plots/",title,"_mean.3x_boxplot.pdf",sep = ""))
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "% mCpG / CpG") +
  ggtitle(label = paste(title," levels",sep = ""),subtitle = "3x coverage") +
  geom_boxplot(width=0.02,outlier.size = .5)
dev.off()

pdf(file = paste("3x/outputTE/plots/",title,"_mean.3x.pdf",sep = ""))
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "% mCpG / CpG") +
  ggtitle(label = paste(title," levels",sep = ""),subtitle = "3x coverage")
dev.off()

## with jitter
pdf(file = paste("3x/outputTE/plots/",title,"_mean.3x.jitter.pdf",sep = ""))
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "% mCpG / CpG") +
  ggtitle(label = paste(title," levels",sep = ""),subtitle = "3x coverage") +
  geom_jitter(width = .1,alpha=.1)
dev.off()


#### Counts ..
title <- "FLI-L1"
ko <- "3x/outputTE/beds/DNMT1ko_FLI-L1.3x_count.bed"
wt <- "3x/outputTE/beds/WT_FLI-L1.3x_count.bed"

wt <- read.delim(file = wt,header=F)
ko <- read.delim(file = ko,header=F)
wt.df <- data.frame(value=wt,variable="WT")
ko.df <- data.frame(value=ko,variable="KO")
head(wt.df)
head(ko.df)

dat <- rbind(wt.df,ko.df)
head(dat)

# jitter
pdf(file = paste("3x/outputTE/plots/",title,"_COUNTS.3x_jitter.pdf",sep = ""))
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "Number of CpGs included per element") +
  ggtitle(label = paste(title," CpGs included",sep = ""),subtitle = "Number of CpG included per element (3x coverage)") +
  geom_jitter(width = 0.2,alpha=0.1)
dev.off()

# no jitter
pdf(file = paste("3x/outputTE/plots/",title,"_COUNTS.3x.pdf",sep = ""))
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "Number of CpGs included per element") +
  ggtitle(label = paste(title," CpGs included",sep = ""),subtitle = "Number of CpG included per element (3x coverage)")
dev.off()

pdf(file = paste("3x/outputTE/plots/",title,"_COUNTS.3x.hist.pdf",sep = ""))
ggplot(dat, aes(fill = variable, x = V1)) +
  geom_histogram(binwidth = 1,alpha=.5,position="identity") +
  ggtitle(label = "Number of CpG called",subtitle = paste("each",title,"(with at least 3x coverage)")) +
  xlab(label = "number of CpGs with at least 3x coverage")
dev.off()

######

##### FLI L1 full fragment #####
## 1x
title <- "FLI-L1"
ko <- "outputTE_1x/beds/DNMT1ko_FLI-L1.mean_rmNAN.bed"
wt <- "outputTE_1x/beds/WT_FLI-L1.mean_rmNAN.bed"

wt <- read.delim(file = wt,header=F)
ko <- read.delim(file = ko,header=F)

wt.df <- data.frame(value=wt,variable="WT")
ko.df <- data.frame(value=ko,variable="KO")
head(wt.df)
head(ko.df)

dat <- rbind(wt.df,ko.df)
head(dat)

pdf(file = paste("outputTE_1x/plots/",title,"_mean.1x_boxplot.pdf",sep = ""))
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "% mCpG / CpG") +
  ggtitle(label = paste(title," levels",sep = ""),subtitle = ">=1x coverage") +
  geom_boxplot(width=0.02,outlier.size = .5)
dev.off()

pdf(file = paste("outputTE_1x/plots/",title,"_mean.1x.pdf",sep = ""))
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "% mCpG / CpG") +
  ggtitle(label = paste(title," levels",sep = ""),subtitle = ">=1x coverage")
dev.off()

######
