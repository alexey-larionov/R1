#########################
# This Script performs the following tasks:
# 1. Links the genotype/sequence data to the clinical data via a "plate_well position" ID
# 2. Creates covariates Eigen 1-5, age, treatment
# 3. runs a regression analysis and outputs: 1) P.LRT; 2) P.Score; 3) P.Wald
#
# Last updated: 11/22/15
#########################

library(statmod)
library(data.table)

RunAnalysis <- function(v) {
	Variant <- v[[1]][1]
	v <- as.numeric(v)[2:513] # first element will be NA in coversion
	v <- ifelse(is.na(v), mean(v, na.rm=T), v) # impute expected score for missing values?
	reg.null <- glm(Y ~ W+Q, family=binomial)
	reg <- glm(Y ~ v + W+Q, family=binomial)
	chi.stat.LRT = 2*(logLik(reg) - logLik(reg.null))
	P.LRT = 1-pchisq(chi.stat.LRT, df=1)
	Score.Z=glm.scoretest(reg.null,v[keep],dispersion=NULL)
	P.Score=2*(1-pnorm(abs(Score.Z)))
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
    r <- c(Variant, round(Estimate,3), round(SE,3), format(P.Wald, scientific=T, digits=3), format(P.LRT, scientific=T, digits=3), format(P.Score, scientific=T, digits=3))
    r
} 


g <- fread("GenotypeData.txt", header=T, sep=" ")
PID <- names(g)
PID <- PID[2:length(PID)]
PID.r <- unlist(lapply(PID, FUN=function(v) { 
	p <- unlist(strsplit(v, "_"))
	w <- as.numeric(substr(p[2], start=2, stop=nchar(p[2])))
	p <- paste(p[1], paste(substr(p[2], 1, 1), w, sep=""), sep="_")
	p
}))


d <- read.table("Clinical.Info.txt", header=T, sep="\t")
d <- d[match(PID.r, d$PID), ]

Y <- d$cc.x
Q <- cbind(d$Eigen_1, d$Eigen_2, d$Eigen_3, d$Eigen_4, d$Eigen_5)
W <- cbind(d$sub_dx_age.x, d$treatment)

keep <- !is.na(Y)


write.table(t(c("Variant", "Estimate", "SE", "P.Wald", "P.LRT", "P.Score")), file="RegressionResults.txt", append=F, quote=F, sep=" ", row.names=F, col.names=F)
system.time(r <- apply(g, 1, FUN=function(G) { 
	write.table(t(RunAnalysis(G)), file="RegressionResults.txt", append=T, quote=F, sep=" ", row.names=F, col.names=F) }))


