### plot CHH statistcs
# Take data from reports from bismark (CX)

wt <- data.frame(unm = 3335314391, m = 13731074)
ko <- data.frame(unm = 3513505002, m = 14153937)

tab <- data.frame(rbind(wt,ko))
rownames(tab) <- c("WT","KO")

tab$perc <- 100*tab$m/(tab$unm+tab$m)

x <- barplot(tab$perc,ylim=c(0,1),col = "black",ylab="% mCHH/CHH")
mtext(text = c("CTR","DNMT1-KO"),side = 1,at=x)
mtext(text = "CHH methylation")
