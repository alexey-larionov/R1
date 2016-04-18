
library(gap)

d <- read.table("RegressionResults.txt", header=T, sep=" ")

d.sig <- d[d$P.LRT< 0.00005, ]

r <- gcontrol2(d$P.LRT[d$SE<.5])  # use SE<0.5 as a proxy for common SNPs
r$lambda



