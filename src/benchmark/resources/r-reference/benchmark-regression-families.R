# Same deterministic input formulas as RegressionFamiliesBenchmark.
.libPaths(c("build/r-library",.libPaths()))
args<-commandArgs(trailingOnly=TRUE)
n<-if(length(args))as.integer(args[1]) else 512L
reps<-if(length(args)>1)as.integer(args[2]) else 7L
batch<-if(length(args)>2)as.integer(args[3]) else 20L
i<-0:(n-1L);z<--3+6*i/(n-1);x<-cos(i*1.7)+.2*z;X<-cbind(1,x)
y<-1+.4*x+sin(z)+.15*cos(i*.37);Y<-cbind(y,-.5+.8*x+.1*sin(i*.8))
l1<-exp(.3+.5*x);l2<-exp(-.2-.3*x);u<-((i*37+11)%%101+.5)/101
cls<-factor(ifelse(u<1/(1+l1+l2),0,ifelse(u<(1+l1)/(1+l1+l2),1,2)))
measure<-function(name,fn) {
  for(j in 1:5)checksum<-fn()
  seconds<-numeric(reps)
  for(j in seq_len(reps)) {start<-as.numeric(Sys.time());for(k in seq_len(batch))checksum<-fn();seconds[j]<-(as.numeric(Sys.time())-start)/batch}
  cat(sprintf("%s,%d,%.9f,%.12f\n",name,n,median(seconds),checksum))
}
cat("operation,rows,median_seconds,checksum\n")
measure("multivariate_ols",function()sum(coef(lm(Y~x))))
measure("multinomial_logit",function(){m<-nnet::multinom(cls~x,trace=FALSE,maxit=2000,reltol=1e-12,abstol=1e-12);stopifnot(m$convergence==0);as.numeric(logLik(m))})
fn<-function(b){r<-y-X%*%b;a<--r/.001;mean(.5*r+.001*(pmax(a,0)+log1p(exp(-abs(a)))))}
gr<-function(b)as.numeric(crossprod(X,plogis(-(y-X%*%b)/.001)-.5)/n)
measure("quantile_smoothed",function(){m<-optim(c(0,0),fn,gr,method="BFGS",control=list(reltol=1e-14,maxit=20000));stopifnot(max(abs(gr(m$par)))<1e-6);sum(m$par)})
measure("quantile_nonsmooth_rq",function()sum(quantreg::rq(y~x,tau=.5)$coef))
bw<-np::npregbw(xdat=z,ydat=y,bws=.4,bandwidth.compute=FALSE,regtype="lc",ckertype="gaussian")
measure("kernel_np",function()sum(fitted(np::npreg(bws=bw,txdat=z,tydat=y))))
measure("supersmoother",function()sum(supsmu(z,y)$y))
# Pure R implementation of the same fixed-bandwidth Robinson/HC3 operations.
robinson<-function() {
  W<-exp(-.5*(outer(z,z,"-")/.4)^2);W<-W/rowSums(W)
  xp<-x-as.numeric(W%*%x);yp<-y-as.numeric(W%*%y);bread<-1/sum(xp*xp);b<-sum(xp*yp)*bread
  g<-as.numeric(W%*%(y-x*b));intercept<-mean(g)
  A<-diag(n)-W;L<-crossprod(xp,A)*bread;B<-A-xp%*%L;res<-as.numeric(B%*%y)
  variance<-sum(L^2*(res/diag(B))^2);vh<-sum(res^2)/sum(B^2)*sum(L^2)
  stopifnot(is.finite(variance),is.finite(vh));intercept+b
}
measure("partially_linear_R_matrix",robinson)
cat("R:",R.version.string,"\n")
print(sapply(c("nnet","quantreg","np"),function(p)as.character(packageVersion(p))))
