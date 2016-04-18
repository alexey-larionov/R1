#########################
# This Script performs the following tasks:
# 1. Links the genotype/sequence data to the clinical data
# 2. Creates covariates matrix for analysis
# 3. runs a regression analysis and outputs: 1) P.LRT; 2) P.Score; 3) P.Wald; 4) firth's logistic regression
#
# Last updated: 02/29/16
#########################

library(statmod)
library(logistf)
#library(data.table)

RunAnalysis <- function(v) {
	med.DP <- median(as.numeric(dp[v,]), na.rm=T)
	med.GQ <- median(as.numeric(gq[v,]), na.rm=T)
	v <- gt[v,]
	Variant <- row.names(v)
	#print(Variant)
	CR <- round(sum(!is.na(v[keep]))/length(v[keep]),3)
	v <- as.numeric(v)
	MAF <- round(sum(v[keep],na.rm=T)/(length(v[keep])*2),3)
	if(CR < CALL.RATE.CUTOFF) {
		r <- c(Variant, CR, MAF, med.DP, med.GQ, rep(NA, 8))
	}
	if(CR >= CALL.RATE.CUTOFF) {
		if(!(MAF >= MAF.LOWER.CUTOFF && MAF <= MAF.UPPER.CUTOFF)) {
			r <- c(Variant, CR, MAF, med.DP, med.GQ, rep(NA, 8))
		}
		if(MAF >= MAF.LOWER.CUTOFF && MAF <= MAF.UPPER.CUTOFF) {
			#v <- ifelse(is.na(v), mean(v, na.rm=T), v) # impute expected score for missing values?
			reg.null <- glm(Y ~ W+Q, family=binomial, subset=!is.na(v))
			reg <- glm(Y ~ v + W+Q, family=binomial)
			chi.stat.LRT = 2*(logLik(reg) - logLik(reg.null))
			P.LRT = 1-pchisq(chi.stat.LRT, df=1)
			if(CR < 1) {
				P.Score<-NA
			}
			if(CR == 1) {
				Score.Z=glm.scoretest(reg.null,v[keep],dispersion=NULL)
				P.Score=2*(1-pnorm(abs(Score.Z)))
			}
			# Firth regression
			reg.f <- logistf(Y~v+W+Q, pl=T)
			if(!"v" %in% reg.f$terms) {
				Estimate.firth <- NA
				SE.firth <- NA
				P.firth <- NA
			}
			if("v" %in% reg.f$terms) {
				Estimate.firth <- reg.f$coef["v"]
				SE.firth <- sqrt(diag(reg.f$var))[2]
				P.firth <- reg.f$prob["v"]
			}
			if (!"v" %in% rownames(summary(reg)$coefficients)){
				Estimate = NA
				SE = NA
				P.Wald = NA
			}
			if ("v" %in% rownames(summary(reg)$coefficients)){
				Estimate = coef(summary(reg))["v","Estimate"]
				SE = coef(summary(reg))["v","Std. Error"]
				P.Wald = coef(summary(reg))["v","Pr(>|z|)"]
		    }
		    r <- c(Variant, CR, MAF, med.DP, med.GQ, round(Estimate,3), round(SE,3), format(P.Wald, scientific=T, digits=3), format(P.LRT, scientific=T, digits=3), format(P.Score, scientific=T, digits=3), round(Estimate.firth,3), round(SE.firth,3), format(P.firth, scientific=T, digits=3))
		}
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

keep <- !is.na(Y)
CALL.RATE.CUTOFF <- 0.95
MAF.LOWER.CUTOFF <- 0.01
MAF.UPPER.CUTOFF <- 0.99

write.table(t(c("Variant", "Call.Rate", "MAF", "Median.DP", "Median.GQ", "Estimate", "SE", "P.Wald", "P.LRT", "P.Score", "Estimate.Firth", "SE.Firth", "P.Firth")), file="RegressionResults.txt", append=F, quote=F, sep=" ", row.names=F, col.names=F)
system.time(r <- lapply(1:nrow(gt), FUN=function(G) { 
	write.table(t(RunAnalysis(G)), file="RegressionResults.txt", append=T, quote=F, sep=" ", row.names=F, col.names=F) }))
