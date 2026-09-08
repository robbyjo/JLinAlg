.libPaths(c('build/r-library',.libPaths()))
# Run from the repository root. Uses only the existing n=160 design and a fixed
# lambda: glmnet objective RSS/(2*n)+lambda*||beta||_1, fixedLassoInf penalty n*lambda.
options(digits=17)
n<-160L;i<-0:(n-1L);z<--3+6*i/(n-1)
X<-cbind(cos(i*.7)+.2*z,sin(i*.23)-.1*z,cos(i*.33),sin(i*.41))
y<-.5+.9*X[,1]-.6*X[,2]+.5*sin(i*1.13)
lambda<-.12;sigma<-.5
fit <- function() {
  g<-glmnet::glmnet(X,y,lambda=lambda,alpha=1,standardize=FALSE,intercept=TRUE,thresh=1e-14)
  beta<-as.numeric(coef(g))[-1]
  f<-selectiveInference::fixedLassoInf(X,y,beta,lambda*n,sigma=sigma,intercept=TRUE,type='partial',alpha=.05)
  list(beta=beta,f=f)
}
loginterval <- function(a,b) {
  if(a>=b)return(-Inf)
  if(a>=0){high<-pnorm(-a,log.p=TRUE);low<-pnorm(-b,log.p=TRUE)}
  else {high<-pnorm(b,log.p=TRUE);low<-pnorm(a,log.p=TRUE)}
  high+log(-expm1(low-high))
}
extract <- function(run,tight=TRUE) {
  f<-run$f;count<-length(f$coef0)
  out<-data.frame(predictor_index=which(run$beta!=0)-1L,beta=f$coef0,sd=f$sd,
    lower=rep(0,count),upper=rep(0,count),p=rep(0,count),
    ci_lower=f$ci[,1],ci_upper=f$ci[,2],grid_lower=f$ci[,1],grid_upper=f$ci[,2])
  for(j in seq_len(count)) {
    sg<-f$sign[j];t<-sg*f$coef0[j];sd<-f$sd[j];lo<-f$vlo[j];hi<-f$vup[j]
    out[j,c('lower','upper')]<-sort(sg*c(lo,hi))
    denom<-loginterval(lo/sd,hi/sd)
    out$p[j]<-2*exp(min(loginterval(lo/sd,t/sd),loginterval(t/sd,hi/sd))-denom)
    if(tight) {
      cdf<-function(mu)exp(loginterval((lo-mu)/sd,(t-mu)/sd)-loginterval((lo-mu)/sd,(hi-mu)/sd))
      roots<-sapply(c(.975,.025),function(target) {
        a<-t-sd;b<-t+sd;step<-sd
        while(cdf(a)<target){step<-step*2;a<-t-step};step<-sd
        while(cdf(b)>target){step<-step*2;b<-t+step}
        root<-uniroot(function(mu)cdf(mu)-target,c(a,b),tol=1e-12)$root
        stopifnot(abs(cdf(root)-target)<1e-8)
        root
      })
      out[j,c('ci_lower','ci_upper')]<-sort(sg*roots)
    }
  }
  out
}
run<-fit();stopifnot(identical(which(run$beta!=0),c(1L,2L)))
ref<-extract(run)
# Independent OLS projection and KKT normalization gate, outside timings.
xc<-scale(X,center=TRUE,scale=FALSE);yc<-y-mean(y);active<-which(run$beta!=0)
stopifnot(max(abs(drop(solve(crossprod(xc[,active]),crossprod(xc[,active],yc)))-ref$beta))<1e-8)
score<-drop(crossprod(xc,yc-xc%*%run$beta))/n
stopifnot(max(abs(score[active]-lambda*sign(run$beta[active])))<1e-7,
  max(abs(score[-active]))<=lambda+1e-7,all(is.finite(as.matrix(ref))),all(ref$p>0))
path<-'src/benchmark/resources/r-reference/selective-n160.csv'
dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
write.csv(ref,path,row.names=FALSE,quote=FALSE)
checksum<-function(table)sum(as.matrix(table[,c('beta','sd','lower','upper','p','ci_lower','ci_upper')]))
repetitions<-as.integer(Sys.getenv('SELECTIVE_REPETITIONS','200'))
stopifnot(repetitions>0)
for(tight in c(FALSE,TRUE)) {
  for(j in 1:3)invisible(extract(fit(),tight))
  times<-numeric(3);checks<-numeric(3)
  for(batch in 1:3)times[batch]<-system.time(for(j in seq_len(repetitions))checks[batch]<-checks[batch]+checksum(extract(fit(),tight)))[['elapsed']]
  cat(sprintf('R_%s repetitions=%d median_seconds=%.9f checksum=%.17g\n',
    if(tight)'native_grid_plus_tight_inversion' else 'native_grid',repetitions,median(times),checks[1]))
}
cat('max_native_grid_vs_tight_CI_difference=',max(abs(as.matrix(ref[,c('ci_lower','ci_upper')])-as.matrix(ref[,c('grid_lower','grid_upper')]))),'\n')
cat(R.version.string,'; glmnet=',as.character(packageVersion('glmnet')),
  '; selectiveInference=',as.character(packageVersion('selectiveInference')),'\n')
