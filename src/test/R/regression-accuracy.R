# Regenerate independent R fixtures; run from repository root.
.libPaths(c("build/r-library", .libPaths()))
options(digits=17)
values <- list()
put <- function(name, x) values[[name]] <<- as.numeric(x)
n <- 160L; i <- 0:(n-1L); z <- -3+6*i/(n-1)
x <- cos(i*1.7)+.2*z; X <- cbind(1,x)
y <- 1+.4*x+sin(z)+.15*cos(i*.37)
Y <- cbind(y,-.5+.8*x+.1*sin(i*.8))
m <- lm(Y ~ x)
put("ols_beta",coef(m)); put("ols_cov",crossprod(resid(m))/df.residual(m))
l1 <- exp(.3+.5*x); l2 <- exp(-.2-.3*x); u <- ((i*37+11)%%101+.5)/101
cls <- ifelse(u<1/(1+l1+l2),0,ifelse(u<(1+l1)/(1+l1+l2),1,2))
m <- nnet::multinom(factor(cls)~x,trace=FALSE,maxit=2000,reltol=1e-12,abstol=1e-12)
put("multinomial_beta",t(coef(m)));put("multinomial_loglik",as.numeric(logLik(m)))
for(tau in c(.25,.5,.75))put(paste0("rq_",tau),quantreg::rq(y~x,tau=tau)$coefficients)
# Smoothed pinball reference matches the exact objective exposed by Java.
for(tau in c(.25,.5,.75)) {
  eps <- 1e-3
  fn <- function(beta){r<-y-X%*%beta;a<--r/eps;mean(tau*r+eps*(pmax(a,0)+log1p(exp(-abs(a)))))}
  gr <- function(beta)as.numeric(crossprod(X,plogis(-(y-X%*%beta)/eps)-tau)/n)
  m<-optim(quantreg::rq(y~x,tau=tau)$coef,fn,gr,method="BFGS",control=list(reltol=1e-14,maxit=20000))
  stopifnot(max(abs(gr(m$par)))<1e-6)
  put(paste0("smooth_rq_",tau),m$par)
}
kernel <- function(v) {sapply(z,function(q){a<--.5*((q-z)/.4)^2;w<-exp(a-max(a));sum(w*v)/sum(w)})}
put("kernel",kernel(y))
put("supsmu",supsmu(z,y)$y)
put("supsmu_weighted",supsmu(z,y,wt=1+(i%%3),bass=5)$y)
put("supsmu_fixed",supsmu(z,y,span=.3)$y)
put("supsmu_periodic",supsmu((z+3)/6,y,periodic=TRUE)$y)
xp<-x-kernel(x);yp<-y-kernel(y);b<-sum(xp*yp)/sum(xp*xp)
g<-kernel(y-x*b);intercept<-mean(g);g<-g-intercept
S<-exp(-.5*(outer(z,z,"-")/.4)^2);S<-S/rowSums(S);A<-diag(n)-S
L<-crossprod(xp,A)/sum(xp*xp);B<-A-xp%*%L
res<-as.numeric(B%*%y)
v<-sum(L^2*(res/diag(B))^2)
vh<-sum(res^2)/sum(B^2)*sum(L^2)
put("robinson_beta",c(intercept,b));put("robinson_slope_cov",v);put("robinson_smooth",g)
put("robinson_homoskedastic_cov",vh);put("robinson_noise_df",sum(B^2))
# A fixed-lambda LASSO example with non-orthogonal predictors and known sigma.
XS<-cbind(cos(i*.7)+.2*z,sin(i*.23)-.1*z,cos(i*.33),sin(i*.41))
ys<-.5+XS[,1]*.9-XS[,2]*.6+.5*sin(i*1.13)
lambda<-.12
gfit<-glmnet::glmnet(XS,ys,lambda=lambda,alpha=1,standardize=FALSE,intercept=TRUE,thresh=1e-14)
beta<-as.numeric(coef(gfit))[-1]
sfit<-selectiveInference::fixedLassoInf(XS,ys,beta,lambda*n,sigma=.5,intercept=TRUE,type="partial",alpha=.05)
put("selective_beta",sfit$coef0);put("selective_p",2*pmin(sfit$pv,1-sfit$pv));put("selective_ci_grid",t(sfit$ci))
# fixedLassoInf uses a coarse interval grid (attained tails need not equal .025).
# Independently invert its reported truncation distribution with tight root tolerance.
log_interval <- function(a,b) {
  if(a>=0){hi<-pnorm(-a,log.p=TRUE);lo<-pnorm(-b,log.p=TRUE)}
  else {hi<-pnorm(b,log.p=TRUE);lo<-pnorm(a,log.p=TRUE)}
  hi+log(-expm1(lo-hi))
}
ci<-matrix(0,length(sfit$coef0),2);limits<-ci
for(j in seq_along(sfit$coef0)) {
  sg<-sfit$sign[j];t<-sg*sfit$coef0[j];s<-sfit$sd[j];lo<-sfit$vlo[j];hi<-sfit$vup[j]
  cdf<-function(mu)exp(log_interval((lo-mu)/s,(t-mu)/s)-log_interval((lo-mu)/s,(hi-mu)/s))
  roots<-vapply(c(.975,.025),function(target){a<-t-s;b<-t+s;step<-s
    while(cdf(a)<target){step<-step*2;a<-t-step};step<-s
    while(cdf(b)>target){step<-step*2;b<-t+step}
    uniroot(function(mu)cdf(mu)-target,c(a,b),tol=1e-12)$root},numeric(1))
  ci[j,]<-sort(sg*roots);limits[j,]<-sort(sg*c(lo,hi))
}
put("selective_ci",t(ci));put("selective_limits",t(limits))
out<-do.call(rbind,lapply(names(values),function(name)data.frame(name=name,index=seq_along(values[[name]])-1,value=values[[name]])))
dir.create("src/test/resources/regression",recursive=TRUE,showWarnings=FALSE)
write.csv(out,"src/test/resources/regression/r-accuracy.csv",row.names=FALSE,quote=FALSE)
cat("Reference fixtures:",nrow(out),"values;",R.version.string,"\n")
print(sapply(c("nnet","quantreg","glmnet","selectiveInference"),function(p)as.character(packageVersion(p))))
