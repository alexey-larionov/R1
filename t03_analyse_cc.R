# ------------------------------------------------------------------------ #

# calculate null model (imputing missing values)
v <- ifelse(is.na(v), mean(v, na.rm=T), v)
reg.null <- glm(Y ~ W+Q, family=binomial)

# Test alternative calculations
q <- d[,c("Eigen_1", "Eigen_2", "Eigen_3")]
w <- d[,c("sub_dx_age", "chemo", "hormone")]
str(q)
str(w)

# ------------------------------------------------------------------------ #
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
# ------------------------------------------------------------------------ #

write.table(t(c("Variant", "Estimate", "SE", "P.Wald", "P.LRT", "P.Score")), file="RegressionResults.txt", append=F, quote=F, sep=" ", row.names=F, col.names=F)

system.time(r <- apply(g, 1, FUN=function(G) { 
  write.table(t(RunAnalysis(G)), file="RegressionResults.txt", append=T, quote=F, sep=" ", row.names=F, col.names=F) }))


a <- matrix(runif(20, 0, 100), nrow=5)
a
apply(a, 1, FUN=function(A) {sum(A)})
