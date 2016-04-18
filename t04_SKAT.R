#############################

Y <- d.covar$cc.x
Q <- cbind(d.covar$Eigen_1.x, d.covar$Eigen_2.x, d.covar$Eigen_3.x, d.covar$Eigen_4.x, d.covar$Eigen_5.x)
W <- cbind(d.covar$sub_dx_age.x, d.covar$treatment.x)

d.a <-read.table("GeneticData/wecare_feb2016_std_filter_VV_biallelic.txt", header=T, sep="\t")

Gene.List <- unique(d.a$SYMBOL)


keep <- !is.na(Y)
CALL.RATE.CUTOFF <- 0.95
MAF.LOWER.LIMIT <- 0.0
MAF.UPPER.LIMIT <- 0.05

write.table(t(c("Gene", "M", "P.SKAT")), file="SKATResults.txt", append=F, quote=F, sep=" ", row.names=F, col.names=F)

system.time(r <- lapply(Gene.List, FUN=function(Gene) { 
  write.table(t(RunSKATAnalysis(Gene)), file="SKATResults.txt", append=T, quote=F, sep=" ", row.names=F, col.names=F) }))

#########################################
