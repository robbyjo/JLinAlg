.libPaths(c('build/r-library',.libPaths()))
suppressPackageStartupMessages(library(lme4))
options(digits=17)
root='src/test/resources/r-reference/remaining-mixed'
dir.create(root,recursive=TRUE,showWarnings=FALSE)
out=paste0('R.version=',getRversion())
put=function(key,value) out <<- c(out,paste0(key,'=',paste(format(value,digits=17,scientific=TRUE,trim=TRUE),collapse=',')))
softplus=function(x) pmax(x,0)+log1p(exp(-abs(x)))
for(kind in c('rare','poisson')) {
 d=read.delim(paste0('src/test/resources/r-reference/quadrature-',kind,'.tsv'))
 write.table(d,file.path(root,paste0(kind,'.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
 fam=if(kind=='rare')binomial() else poisson()
 fit=glmer(y~x+offset(offset)+(1|group),d,family=fam,nAGQ=1,
  control=glmerControl(optimizer='bobyqa',tolPwrss=1e-12,optCtrl=list(rhoend=1e-10)))
 # Independently evaluate EXACTLY the first-order Laplace approximation, not AGQ.
 X=model.matrix(~x,d);groups=split(seq_len(nrow(d)),d$group)
 objective=function(par) {
  v=exp(par[3]);eta=drop(X%*%par[1:2])+d$offset
  sum(vapply(groups,function(rows) {
   e=eta[rows];y=d$y[rows]
   mu=function(b)if(kind=='rare')plogis(e+b) else exp(e+b)
   score=function(b)sum(y-mu(b))-b/v
   b=uniroot(score,c(-100,100),tol=1e-12)$root
   m=mu(b);curvature=1/v+sum(if(kind=='rare')m*(1-m) else m)
   mass=if(kind=='rare')sum(y*(e+b)-softplus(e+b)) else sum(dpois(y,m,log=TRUE))
   -(mass-b*b/(2*v)-.5*log(v)-.5*log(curvature))
  },0.))
 }
 start=c(fixef(fit),log(as.data.frame(VarCorr(fit))$vcov[1]))
 exact=optim(start,objective,method='BFGS',control=list(reltol=1e-13,ndeps=rep(1e-5,3)),hessian=TRUE)
 stopifnot(exact$convergence==0)
 put(paste0(kind,'.beta'),exact$par[1:2]);put(paste0(kind,'.variance'),exp(exact$par[3]))
 put(paste0(kind,'.LL'),-exact$value);put(paste0(kind,'.covariance'),as.vector(t(solve(exact$hessian)[1:2,1:2])))
 put(paste0(kind,'.lme4LL'),as.numeric(logLik(fit)));put(paste0(kind,'.lme4Beta'),fixef(fit))
 cat(kind,'independent LL',-exact$value,'lme4 LL',as.numeric(logLik(fit)),'\n')
}
# Unbalanced Gaussian REML, nonzero fixed means, independently optimized dense likelihood.
set.seed(2931);group=rep(1:12,c(2,3,4,5,6,2,7,3,4,5,6,7));n=length(group)
x=sin(seq_len(n)*.37);y=1+.6*x+rnorm(12,sd=.8)[group]+rnorm(n,sd=.4)
d=data.frame(y=y,x=x,group=group);write.table(d,file.path(root,'reml.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
X=model.matrix(~x,d);K=outer(group,group,'==')*1
for(method in c('REML','ML')) {
 likelihood=function(par,details=FALSE) {
  v=exp(par);V=v[1]*K+diag(v[2],n);C=chol(V);Vi=chol2inv(C);H=crossprod(X,Vi%*%X);beta=solve(H,crossprod(X,Vi%*%y));r=y-drop(X%*%beta)
  value=.5*((n-if(method=='REML')ncol(X) else 0)*log(2*pi)+2*sum(log(diag(C)))+
   if(method=='REML')as.numeric(determinant(H,logarithm=TRUE)$modulus) else 0)
  value=value+.5*sum(r*(Vi%*%r))
  if(details)list(beta=beta,covariance=solve(H),LL=-value) else value
 }
 # Parentheses avoid R's if/else expression binding into the final quadratic.
 fitted=optim(log(c(.7,.2)),likelihood,method='BFGS',control=list(reltol=1e-13))
 check=likelihood(fitted$par,TRUE);put(paste0('reml.',method,'.variance'),exp(fitted$par));
 put(paste0('reml.',method,'.beta'),check$beta);put(paste0('reml.',method,'.covariance'),as.vector(t(check$covariance)));put(paste0('reml.',method,'.LL'),check$LL)
 lmerFit=lmer(y~x+(1|group),d,REML=method=='REML')
 stopifnot(abs(as.numeric(logLik(lmerFit))-check$LL)<1e-7)
}
# Fixed-dispersion PQL: independent dense R working REML with residual scale ONE.
# MASS::glmmPQL estimates an additional residual scale, so it is not this estimand.
y=c(0,1,1,1,2,3,3,4,5,7,8,9);n=length(y);X=matrix(1,n,1);K=kronecker(diag(4),matrix(1,3,3))
eta=rep(log(mean(y)),n);v=.5
for(iteration in 1:100) {
 mu=exp(eta);z=eta+(y-mu)/mu;R=diag(1/mu)
 working=function(logv,details=FALSE) {
  V=exp(logv)*K+R;Vi=solve(V);H=crossprod(X,Vi%*%X);beta=solve(H,crossprod(X,Vi%*%z));r=z-drop(X%*%beta)
  nll=.5*((n-1)*log(2*pi)+as.numeric(determinant(V,logarithm=TRUE)$modulus)+log(H[1,1])+sum(r*(Vi%*%r)))
  if(details)list(beta=beta,random=drop(exp(logv)*K%*%Vi%*%r),LL=-nll,covariance=solve(H)) else nll
 }
 optimum=optimize(working,c(-20,5),tol=1e-10);solved=working(optimum$minimum,TRUE)
 nextEta=drop(X%*%solved$beta)+solved$random;nextV=exp(optimum$minimum)
 change=max(abs(nextEta-eta)/(1+abs(eta)),abs(nextV-v)/(1+v));eta=nextEta;v=nextV
 if(change<1e-9)break
}
stopifnot(change<1e-9)
put('pql.beta',solved$beta);put('pql.variance',v);put('pql.workingLL',solved$LL)
put('pql.random',solved$random);put('pql.covariance',solved$covariance)
# Known numerator matrix: unrelated founders a,b; c,d full siblings;
# e=c x d is inbred; f has known sire e and an unrelated unknown dam.
parents=list(c(0,0),c(0,0),c(1,2),c(1,2),c(3,4),c(5,0));A=matrix(0,6,6)
for(i in 1:6) {
 pa=parents[[i]];known=pa[pa>0]
 if(i>1)for(j in seq_len(i-1))A[i,j]=A[j,i]=sum(A[known,j])/2
 A[i,i]=1+if(all(pa>0))A[pa[1],pa[2]]/2 else 0
}
put('pedigree.A',as.vector(t(A)));put('pedigree.precision',as.vector(t(solve(A))))
# Finite grouped-binomial optimum with a positive tail rounded to one in double.
tailScore=function(b) 1e22*plogis(-b-50)-3*plogis(b)
tailBeta=uniroot(tailScore,c(-10,10),tol=1e-13)$root
put('tail.beta',tailBeta)
put('tail.LL',-1e22*log1p(exp(-tailBeta-50))-3*log1p(exp(tailBeta)))
put('tail.pql.response',51);put('tail.pql.covariance',exp(50)/1e22)
writeLines(out,file.path(root,'reference.properties'))
