.libPaths(c('build/r-library', .libPaths()))
# Run from repository root. Fixtures use the FULL conditional probability mass,
# including binomial coefficients/Poisson factorials, integrated over N(0,sd^2).
options(digits = 17)
root <- 'src/test/resources/r-reference'
softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
group_integral <- function(y, eta, sd, family, trials = rep(1, length(y))) {
  logf <- function(u) {
    e <- eta + sd * u
    if (family == 'binomial') {
      s <- round(y * trials)
      sum(lchoose(trials, s) - s * softplus(-e) - (trials-s) * softplus(e)) + dnorm(u, log=TRUE)
    } else sum(dpois(y, exp(e), log=TRUE)) + dnorm(u, log=TRUE)
  }
  if (sd == 0) return(logf(0) - dnorm(0, log=TRUE))
  m <- optimize(function(u) -logf(u), c(-50, 50), tol=1e-12)$minimum
  shift <- logf(m)
  # Split at the mode so narrow posterior peaks cannot be missed by integrate.
  f <- function(v) vapply(v, function(z) exp(logf(z+m)-shift), 0.0)
  a <- integrate(f, -Inf, 0, rel.tol=2e-12, subdivisions=1000L)
  b <- integrate(f, 0, Inf, rel.tol=2e-12, subdivisions=1000L)
  shift + log(a$value + b$value)
}
cases <- list(
  list('balanced100', 'binomial', rep(c(0,1),50), 0, 1),
  list('balanced10000', 'binomial', rep(c(0,1),5000), 0, 1),
  list('rare1of100', 'binomial', c(1,rep(0,99)), -5, 2),
  list('allzero100', 'binomial', rep(0,100), -4, 3),
  list('rare1of1000', 'binomial', c(1,rep(0,999)), -7, 2.5),
  list('countlarge', 'poisson', rep(c(90,110),50), 4, 1.5),
  list('countzero', 'poisson', rep(0,30), -2, 2),
  list('countmixed', 'poisson', c(0,0,0,1,2), -2, 1.2))
fixed <- do.call(rbind, lapply(cases, function(z) data.frame(
  name=z[[1]], family=z[[2]], n=length(z[[3]]), successes=sum(z[[3]]),
  beta=z[[4]], sd=z[[5]], logLikelihood=group_integral(z[[3]], rep(z[[4]],length(z[[3]])),z[[5]],z[[2]]))))
write.table(fixed, file.path(root,'quadrature-integrals.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
print(fixed)
if ('integrals-only' %in% commandArgs(TRUE)) quit(status=0)
suppressPackageStartupMessages(library(lme4))

set.seed(182753)
properties <- c(paste0('R.version=', getRversion()), paste0('lme4.version=',packageVersion('lme4')))
for (kind in c('binary','rare','binomial','poisson')) {
  ng <- 32L; per <- 12L; n <- ng*per
  group <- factor(rep(seq_len(ng), each=per))
  x <- rep(seq(-1.5,1.5,length.out=per), ng)
  offset <- .15 * sin(seq_len(n))
  eta <- (if(kind=='rare') -3.3 else if(kind=='poisson') .4 else -.3) + .7*x +
    rep(rnorm(ng,sd=if(kind=='rare') 1.4 else .9),each=per) + offset
  trials <- if(kind=='binomial') rep(c(5,10,20),length.out=n) else rep(1,n)
  y <- if(kind=='poisson') rpois(n,exp(eta)) else rbinom(n,trials,plogis(eta))/trials
  d <- data.frame(y,x,group,offset,trials)
  family <- if(kind=='poisson') poisson() else binomial()
  elapsed <- system.time(fit <- glmer(y~x+offset(offset)+(1|group), data=d, family=family,
    weights=trials, nAGQ=25, control=glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=200000))))[['elapsed']]
  b <- unname(fixef(fit)); sd <- sqrt(as.numeric(VarCorr(fit)$group))
  # lme4's Hessian is for deviance in (SD, beta0, beta1). Transform the
  # inverse observed information into (beta0, beta1, variance).
  joint <- solve(fit@optinfo$derivs$Hessian/2)[c(2,3,1),c(2,3,1)] * outer(c(1,1,2*sd),c(1,1,2*sd))
  ll <- sum(vapply(split(seq_len(n),group), function(i)
    group_integral(y[i],b[1]+b[2]*x[i]+offset[i],sd,if(kind=='poisson') 'poisson' else 'binomial',trials[i]),0.0))
  # glmer nAGQ>1 reports a deviance-residual-based objective for counts;
  # freeze the independently integrated FULL LL, not that shifted logLik.
  prefix <- function(key,value) paste0(kind,'.',key,'=',format(value,digits=17,scientific=FALSE,trim=TRUE))
  selected <- c(1L,193L,384L)
  marginal <- vapply(selected,function(i) if(kind=='poisson') exp(b[1]+b[2]*x[i]+offset[i]+sd^2/2)
    else exp(group_integral(1,b[1]+b[2]*x[i]+offset[i],sd,'binomial')),0.0)
  properties <- c(properties,prefix('intercept',b[1]),prefix('slope',b[2]),prefix('sd',sd),
    prefix('logLikelihood',ll),prefix('se0',sqrt(vcov(fit)[1,1])),prefix('se1',sqrt(vcov(fit)[2,2])),
    prefix('cov01',vcov(fit)[1,2]),prefix('cov02',joint[1,3]),prefix('cov12',joint[2,3]),
    prefix('varVariance',joint[3,3]),prefix('mean0',marginal[1]),prefix('mean192',marginal[2]),prefix('mean383',marginal[3]),
    prefix('lme4Seconds',elapsed),prefix('lme4ReportedLogLik',as.numeric(logLik(fit))))
  write.table(d,file.path(root,paste0('quadrature-',kind,'.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
}
writeLines(properties,file.path(root,'quadrature-fits.properties'))
print(fixed)
cat(paste(properties,collapse='\n'),'\n')
