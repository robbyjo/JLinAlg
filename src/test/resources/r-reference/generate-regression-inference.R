# Independent R references. Run from the repository root with Rscript.
dir <- "src/test/resources/r-reference/regression-inference"
dir.create(dir, recursive=TRUE, showWarnings=FALSE)
.libPaths(c("build/r-library",.libPaths()))
set.seed(20260919)
n <- 240
c <- rnorm(n); z1 <- rnorm(n); z2 <- rnorm(n); u <- rnorm(n)
d <- .5*c + .9*z1 - .6*z2 + u
y <- 1 + .3*c + .8*d + .7*u + rnorm(n,sd=.6)
trials <- rep(c(3,7,11,5),length.out=n)
offset <- .15*sin(seq_len(n))
binary <- rbinom(n,trials,pnorm(-.4+.55*c-.25*z1+offset))/trials
cluster <- rep(seq_len(n/4),each=4)
data <- data.frame(y,c,z1,z2,d,binary,trials,offset,cluster)
con <- file(file.path(dir,"reference.properties"),"w")
put <- function(key,x) writeLines(paste0(key,"=",paste(sprintf("%.17g",x),collapse=",")),con)
writeLines(paste0("R.version=",getRversion()),con)
latent <- 1+.4*c-.6*z1+rnorm(n,sd=.8)
tobit <- pmax(-.4,pmin(2.5,latent)); cens <- ifelse(latent<=-.4,-1,ifelse(latent>=2.5,1,0))
time0 <- exp(.8+.3*c-.2*z1+.65*log(-log(runif(n))))
limit <- exp(rnorm(n,1.2,.7)); time <- pmin(time0,limit); censTime <- as.integer(time0>limit)
likelihoodData <- data.frame(tobit,c,z1,cens,time,censTime)
write.table(likelihoodData,file.path(dir,"censored.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
library(survival)
left <- ifelse(cens<0,NA,tobit);right<-ifelse(cens>0,NA,tobit)
fit <- survreg(Surv(left,right,type="interval2")~c+z1,dist="gaussian")
put("censored.gaussian.beta",coef(fit));put("censored.gaussian.scale",fit$scale);put("censored.gaussian.cov",t(vcov(fit)));put("censored.gaussian.logLik",as.numeric(logLik(fit)))
for(dist in c("lognormal","weibull","exponential")) {
 fit<-survreg(Surv(time,1-censTime)~c+z1,dist=dist)
 put(paste0("censored.",dist,".beta"),coef(fit));put(paste0("censored.",dist,".scale"),fit$scale)
 put(paste0("censored.",dist,".cov"),t(vcov(fit)));put(paste0("censored.",dist,".logLik"),as.numeric(logLik(fit)))
}
ordinal <- as.integer(cut(.7*c-.3*z1+rlogis(n),breaks=c(-Inf,-.6,.4,1.5,Inf)))-1
write.table(data.frame(ordinal,c,z1),file.path(dir,"ordinal.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
for(link in c("logistic","probit")) {
 fit<-MASS::polr(ordered(ordinal)~c+z1,method=link,Hess=TRUE,control=list(reltol=1e-12,maxit=2000))
 key<-if(link=="logistic")"logit" else link
 put(paste0("ordinal.",key,".beta"),coef(fit));put(paste0("ordinal.",key,".thresholds"),fit$zeta)
 put(paste0("ordinal.",key,".cov"),t(vcov(fit)));put(paste0("ordinal.",key,".logLik"),as.numeric(logLik(fit)))
}
rare <- rbinom(n,1,plogis(-2+.6*c-.2*z1))
weight <- exp(.3*c);stratum<-ceiling(cluster/12); count<-rpois(n,exp(.4+.2*c))
surveyData<-data.frame(y,c,z1,rare,count,weight,stratum,cluster)
write.table(surveyData,file.path(dir,"survey.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
fit<-glm(rare~c+z1,binomial(),control=glm.control(epsilon=1e-12))
put("rare.ml.beta",coef(fit));put("rare.ml.cov",t(vcov(fit)))
X<-model.matrix(fit);mu<-fitted(fit);ybar<-mean(rare)
for(tau in c(ybar,.04)) {
 # King-Zeng ISQ equation 11: evaluate independently with dense cross-products.
 w1<-tau/ybar;w0<-(1-tau)/(1-ybar);W<-mu*(1-mu)*ifelse(rare==1,w1,w0)
 Q<-solve(crossprod(X,X*W));qdiag<-rowSums((X%*%Q)*X)
 xi<-.5*qdiag*((1+w1)*mu-w1);bias<-Q%*%crossprod(X,W*xi)
 beta<-coef(fit)-drop(bias);beta[1]<-beta[1]+qlogis(tau)-qlogis(ybar)
 put(if(tau==ybar)"rare.bias.beta" else "rare.prior.beta",beta)
}
library(survey)
sampling<-svydesign(ids=~cluster,strata=~stratum,weights=~weight,data=surveyData,nest=TRUE)
for(name in c("gaussian","logit","probit","poisson")) {
 family<-switch(name,gaussian=gaussian(),logit=quasibinomial(),probit=quasibinomial("probit"),poisson=quasipoisson())
 form<-if(name=="gaussian")y~c+z1 else if(name=="poisson")count~c+z1 else rare~c+z1
 fit<-svyglm(form,sampling,family=family,control=glm.control(epsilon=1e-12))
 put(paste0("survey.",name,".beta"),coef(fit));put(paste0("survey.",name,".cov"),t(vcov(fit)))
}
writeLines(paste0("survey.version=",packageVersion("survey")),con)
writeLines(paste0("survival.version=",packageVersion("survival")),con)
writeLines(paste0("MASS.version=",packageVersion("MASS")),con)
close(con)
