## Violin plot coverage
library(ggplot2)
setwd(dir = "~/Dropbox (MN)/Per/PhD/Projects/DNAmeth/hNES DNMT1 KO/2018.08.30_DNMT1L1_figures_scripts/WGBS/")

wt <- read.delim(file = "data/WT_bismark_1000bin.mean_only.bed",header=F)
ko <- read.delim(file = "data/DNMT1ko_bismark_1000bin.mean_only.bed",header=F)

wt.df <- data.frame(value=wt,variable="WT")
ko.df <- data.frame(value=ko,variable="KO")
head(wt.df)
head(ko.df)

dat <- rbind(wt.df,ko.df)
head(dat)

pdf(file = "analysis/global/global_vioplot_boxplot_allCpG.pdf")
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T, adjust = 2) +
  theme(panel.grid.major = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "% mCpG / CpG") +
  ggtitle(label = "Global levels (1000bp bins)") +
  geom_boxplot(width=0.02,outlier.size = .5)
dev.off()

pdf(file = "analysis/global/global_vioplot_allCpG.pdf")
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T,adjust=2) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "% mCpG / CpG") +
  ggtitle(label = "Global levels (1000bp bins)")
dev.off()


############

###### 3x coverage ######

wt <- read.delim(file = "data/3x/WT_bismark_1000bin.mean_only.3x.bed",header=F)
ko <- read.delim(file = "data/3x/DNMT1ko_bismark_1000bin.mean_only.3x.bed",header=F)

wt.df <- data.frame(value=wt,variable="WT")
ko.df <- data.frame(value=ko,variable="KO")
head(wt.df)
head(ko.df)

dat <- rbind(wt.df,ko.df)
head(dat)

pdf(file = "analysis/global/3x_global_vioplot_boxplot.allCpG.pdf")
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T,adjust=2.5) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "% mCpG / CpG") +
  ggtitle(label = "Global levels (1000bp bins)",subtitle = "3x coverage") +
  geom_boxplot(width=0.02,outlier.size = .5)
dev.off()

pdf(file = "analysis/global/3x_global_vioplot_allCpG.pdf")
ggplot(dat, aes(y = V1, x = variable)) +
  geom_violin(draw_quantiles = T,adjust=2.5) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun.y="mean", geom="point") +
  ylab(label = "% mCpG / CpG") +
  ggtitle(label = "Global levels (1000bp bins)",subtitle = "3x coverage")
dev.off()
