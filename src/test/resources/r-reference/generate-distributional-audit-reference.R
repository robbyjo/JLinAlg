.libPaths(c('build/r-library', .libPaths()))
suppressPackageStartupMessages(library(betareg))
options(digits=17)
root <- 'src/test/resources/r-reference'
densities <- list()
for (phi in c(.1,10,1e8,1e14)) for (mu in c(.0001,.4,.9999)) {
  y <- mu
  densities[[length(densities)+1]] <- data.frame(family='beta',y=y,mu=mu,parameter=phi,ll=dbeta(y,mu*phi,(1-mu)*phi,log=TRUE))
}
for (size in c(.01,10,1e8,1e14)) for (mu in c(.01,2,1e5)) for(y in c(0,2,1e5))
  densities[[length(densities)+1]] <- data.frame(family='nb',y=y,mu=mu,parameter=size,ll=dnbinom(y,mu=mu,size=size,log=TRUE))
write.table(do.call(rbind,densities),file.path(root,'distributional-audit-densities.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
set.seed(90317)
n <- 1500; x <- seq(-1,1,length.out=n); z <- sin(seq_len(n)*.31)
mu <- plogis(-.7+2.4*x); phi <- exp(2.1+.8*z)
y <- rbeta(n,mu*phi,(1-mu)*phi)
d <- data.frame(y=y,x=x,z=z)
write.table(d,file.path(root,'distributional-audit-beta.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
run <- function() betareg(y~x|z,d,link='logit',link.phi='log',type='ML',control=betareg.control(maxit=1000,reltol=1e-10))
fit <- run(); stopifnot(fit$converged)
times <- replicate(5,system.time({checked <- run();stopifnot(checked$converged)})['elapsed'])
values <- c(logLik=as.numeric(logLik(fit)),setNames(coef(fit),paste0('coefficient',0:3)),seconds=median(times),checksum=as.numeric(logLik(fit))+sum(coef(fit)))
writeLines(c(paste0('# R ',getRversion(),'; betareg ',packageVersion('betareg')),paste(names(values),format(values,digits=17,scientific=TRUE),sep='=')),file.path(root,'distributional-audit-beta.properties'))
print(values)
suppressPackageStartupMessages(library(glmmTMB))
nd <- d[1:400,]; nd$x <- seq(-1,1,length.out=400)
set.seed(90318); nd$y <- rnbinom(400,mu=exp(.4+.8*nd$x),size=exp(.8+.5*nd$z))
nbfit <- glmmTMB(y~x,dispformula=~z,family=nbinom2,data=nd,
  control=glmmTMBControl(optimizer=optim,optArgs=list(method='BFGS'),optCtrl=list(maxit=2000,reltol=1e-14)))
stopifnot(nbfit$fit$convergence==0,nbfit$sdr$pdHess)
write.table(nd,file.path(root,'distributional-audit-nb.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
nbv <- c(logLik=as.numeric(logLik(nbfit)),setNames(c(fixef(nbfit)$cond,fixef(nbfit)$disp),paste0('coefficient',0:3)),setNames(as.vector(t(vcov(nbfit,full=TRUE))),paste0('covariance',0:15)))
writeLines(c(paste0('# R ',getRversion(),'; glmmTMB ',packageVersion('glmmTMB')),paste(names(nbv),format(nbv,digits=17,scientific=TRUE),sep='=')),file.path(root,'distributional-audit-nb.properties'))
# Independent likelihood maximization for a large-mean hurdle; information is
# the zero-truncated-Poisson variance times positive mass, not Poisson variance.
hy <- c(rep(0,40),990+rep(0:20,length.out=160)); p0 <- mean(hy==0)
objective <- function(eta) {mu <- exp(eta); -sum(ifelse(hy==0,log(p0),log1p(-p0)+dpois(hy,mu,log=TRUE)-log(-expm1(-mu))))}
hf <- optimize(objective,c(log(900),log(1100)),tol=1e-12)
write.table(data.frame(y=hy),file.path(root,'distributional-audit-hurdle.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(paste(c('mu','zero','logLik'),format(c(exp(hf$minimum),p0,-hf$objective),digits=17,scientific=TRUE),sep='='),file.path(root,'distributional-audit-hurdle.properties'))
