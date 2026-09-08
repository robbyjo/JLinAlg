.libPaths(c('build/r-library', .libPaths()))
suppressPackageStartupMessages(library(lme4))
options(digits=17)
i <- 0:39
d <- data.frame(y=2+2*sin(i/5)+cos(2.1*i),g=factor(i%/%5))
fit <- lmer(y~1+(1|g),d,control=lmerControl(optimizer='bobyqa'))
values <- c(randomVariance=as.numeric(VarCorr(fit)[[1]]),residualVariance=sigma(fit)^2,
  logLik=as.numeric(logLik(fit)))
Z <- model.matrix(~0+g,d); X <- matrix(1,40,1)
nll <- function(v) {
  V <- v[1]*tcrossprod(Z)+diag(v[2],40)
  inv <- chol2inv(chol(V)); information <- as.numeric(crossprod(X,inv%*%X))
  beta <- as.numeric(crossprod(X,inv%*%d$y))/information
  residual <- d$y-beta
  .5*(39*log(2*pi)+as.numeric(determinant(V,logarithm=TRUE)$modulus)+log(information)+sum(residual*(inv%*%residual)))
}
bounded <- optim(c(.05,.05),nll,method='L-BFGS-B',lower=c(.01,.01),upper=c(.1,.1))
values <- c(values,boundedRandom=bounded$par[1],boundedResidual=bounded$par[2],boundedLogLik=-bounded$value)
for(sign in c(1,-1)) {
  i <- 0:79; j <- i%%8; g <- i%/%8; x <- j-3.5
  dd <- data.frame(y=2+.3*x+2*sin(g)*(1+sign*.5*x)+ifelse(j%%4==0|j%%4==3,.3,-.3),x=x,g=factor(g))
  fm <- lmer(y~x+(1+x|g),dd,REML=FALSE,control=lmerControl(optimizer='bobyqa',optCtrl=list(rhoend=1e-10)))
  devfun <- lmer(y~x+(1+x|g),dd,REML=FALSE,devFunOnly=TRUE)
  theta <- getME(fm,'theta'); initial <- c(theta[1],sqrt(sum(theta[2:3]^2)))
  profiled <- function(rho) nlminb(initial,function(t) devfun(c(t[1],rho*t[2],sqrt(1-rho*rho)*t[2])),lower=c(0,0),
    control=list(rel.tol=1e-12,eval.max=5000,iter.max=3000))$objective
  target <- profiled(sign)+qchisq(.95,1)
  crossing <- uniroot(function(rho) profiled(rho)-target,sort(c(sign,sign*.99)),tol=1e-10)$root
  values[paste0('correlation',if(sign==1)'Lower' else 'Upper')] <- crossing
}
writeLines(c(paste0('# R ',getRversion(),'; lme4 ',packageVersion('lme4'),'; bounded likelihood independently evaluated'),
  paste(names(values),format(values,digits=17,scientific=TRUE),sep='=')),
  'src/test/resources/r-reference/mixed-review.properties')
print(values)
