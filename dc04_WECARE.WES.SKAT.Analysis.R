#########################
# This Script performs a SKAT analysis by gene:
# 
# Last updated: 02/27/16
#########################

library(SKAT)

RunSKATAnalysis <- function(Gene) {
	Variant.List <- d.a$RawVarID[d.a$SYMBOL==Gene]
	G <- gt[row.names(gt)%in%Variant.List,]
	if(nrow(G) > 1) { 
		CR <- apply(G, 1, FUN=function(v) { round(sum(!is.na(v[keep]))/length(v[keep]),3) })
		G <- G[CR>CALL.RATE.CUTOFF,]
		if(nrow(G) > 1 ) {
			MAF <- apply(G, 1, FUN=function(v) { round(sum(v[keep],na.rm=T)/(length(v[keep])*2),3) })
			G <- G[unlist(lapply(MAF, FUN=function(v) { (v > MAF.LOWER.LIMIT && v < MAF.UPPER.LIMIT) || (v < (1-MAF.LOWER.LIMIT) && v > (1-MAF.UPPER.LIMIT))})),]
		}
	}
	if(nrow(G) <= 1) { r <- c(as.character(Gene), nrow(G), NA) }
	if(nrow(G) > 1) {
		skat.reg.null <- SKAT_Null_Model(Y ~ W+Q, out_type="D")
		skat.reg <- SKAT(Z=t(as.matrix(G)), obj=skat.reg.null, kernel="linear.weighted", method="optimal.adj")
		M <- skat.reg$param$n.marker
		P.SKAT <- skat.reg$p.value
		r <- c(as.character(Gene), M, format(P.SKAT, scientific=T, digits=3))
	}
	r
} 

############# data input
gt <- read.table("GeneticData/wecare_feb2016_std_filter_GT_biallelic_add.txt", header=T, sep="\t")
dim(gt)

gq <- read.table("GeneticData/wecare_feb2016_std_filter_GQ_biallelic.txt", header=T, sep="\t")
row.names(gq) <- gq$RawVarID
gq <- gq[,2:ncol(gq)]
dim(gq)

dp <- read.table("GeneticData/wecare_feb2016_std_filter_DP_biallelic.txt", header=T, sep="\t")
row.names(dp) <- dp$RawVarID
dp <- dp[,2:ncol(dp)]
dim(dp)

###### filter gt by gq and dp:
NA -> gt[gq == 0]
NA -> gt[dp < 5]


d.IDs <- read.table("GeneticData/samples_ids.txt", header=T, sep="\t")
d.covar <- read.table("DemographicVariables/WECARE.Exome.DemographicVariables.txt", header=T, sep="\t")
dim(d.covar)
d.covar <- d.covar[match(d.IDs$gwas_id, d.covar$labid.x), ]
dim(d.covar)


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

