.libPaths(c('build/r-library', .libPaths()))
# Independent fitting/inference references and warmed CPU timing. No installs.
# Run from repository root. Optional arguments: fixture output, benchmark rows.
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else 'src/test/resources/fitting-audit-v030/reference.properties'
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
library(glmnet)
options(digits=17)
make_data <- function(n) {
  i <- 0:(n-1)
  x <- cbind(sin(i*.37), cos(i*.19))
  w <- 1+(i%%5)
  offset <- .2*sin(i*.11)
  y <- 1+.7*x[,1]-.4*x[,2]+.3*sin(i*1.7)+offset
  list(x=x, X=cbind(1,x), w=w, offset=offset, y=y, id=i%/%3,
       counts=as.numeric((i*7+3)%%11), gamma=exp(.3+.5*x[,1]+.7*sin(i)))
}
d <- make_data(60)
values <- list()
save_value <- function(key,value) values[[key]] <<- paste(format(as.numeric(value), digits=17, scientific=TRUE, trim=TRUE),collapse=',')
save_fit <- function(key,f) {
  save_value(paste0(key,'.beta'),coef(f));save_value(paste0(key,'.covariance'),t(vcov(f)))
  save_value(paste0(key,'.p'),coef(summary(f))[,4]);save_value(paste0(key,'.ll'),logLik(f))
}
ols <- lm(d$y ~ d$x, weights=d$w, offset=d$offset)
save_fit('ols',ols)
gaussian <- glm(d$y ~ d$x,weights=d$w,offset=d$offset,family=gaussian())
save_fit('gaussian',gaussian)
poisson <- glm(d$counts ~ d$x,weights=d$w,offset=d$offset,family=poisson(),control=glm.control(epsilon=1e-12,maxit=200))
stopifnot(poisson$converged);save_fit('poisson',poisson)
gamma <- glm(d$gamma ~ d$x,family=Gamma(link='log'),control=glm.control(epsilon=1e-12,maxit=200))
stopifnot(gamma$converged);save_value('gamma.beta',coef(gamma));save_value('gamma.dispersion',summary(gamma)$dispersion)
# Unlike stats::Gamma()'s AIC approximation, maximize the actual Gamma density.
profile <- optimize(function(lp) -sum(dgamma(d$gamma,shape=exp(-lp),scale=fitted(gamma)*exp(lp),log=TRUE)),c(-10,5),tol=1e-12)
save_value('gamma.profile.ll',-profile$objective)
ig <- glm(d$gamma ~ d$x, family=inverse.gaussian(link='log'),control=glm.control(epsilon=1e-12,maxit=200))
stopifnot(ig$converged);save_value('ig.beta',coef(ig));save_value('ig.covariance',t(vcov(ig)))
phi <- deviance(ig)/length(d$gamma)
save_value('ig.profile.ll',sum(-.5*log(2*pi*phi*d$gamma^3)-(d$gamma-fitted(ig))^2/(2*phi*d$gamma*fitted(ig)^2)))
save_value('binomial.rare.ll',sum(dbinom(1:3,1e14,2e-14,log=TRUE)))
tail_score <- function(b) 1e22*plogis(-b-50)-3*plogis(b)
tail_beta <- uniroot(tail_score,c(-10,10),tol=1e-13)$root
tail_ll <- -1e22*log1p(exp(-tail_beta-50))-3*log1p(exp(tail_beta))
tail_info <- 1e22*dlogis(tail_beta+50)+3*dlogis(tail_beta)
save_value('binomial.tail',c(tail_beta,tail_ll,1/sqrt(tail_info)))

# Fixed-correlation Gaussian GEE, phi=1: exact GLS and cluster sandwich.
fixed_gee <- function(d,rho=.2) {
  p <- ncol(d$X); B <- matrix(0,p,p); rhs <- rep(0,p)
  groups <- split(seq_along(d$y),d$id)
  inv <- lapply(groups,function(j) {
    R <- matrix(rho,length(j),length(j));diag(R)<-1
    solve(R / sqrt(outer(d$w[j],d$w[j])))
  })
  for(k in seq_along(groups)){j<-groups[[k]]; X<-d$X[j,,drop=FALSE];B<-B+crossprod(X,inv[[k]]%*%X);rhs<-rhs+drop(crossprod(X,inv[[k]]%*%(d$y[j]-d$offset[j])))}
  bread <- solve(B);beta <- drop(bread%*%rhs);meat <- matrix(0,p,p)
  for(k in seq_along(groups)){j<-groups[[k]];u<-crossprod(d$X[j,,drop=FALSE],inv[[k]]%*%(d$y[j]-d$offset[j]-d$X[j,,drop=FALSE]%*%beta));meat<-meat+tcrossprod(u)}
  list(beta=beta,cov=bread%*%meat%*%bread)
}
gee <- fixed_gee(d);save_value('gee.beta',gee$beta);save_value('gee.covariance',t(gee$cov))
lambda <- c(.2,.08,.02)
lasso <- glmnet(d$x,d$y,weights=d$w,alpha=1,lambda=lambda,thresh=1e-14,maxit=100000)
stopifnot(lasso$jerr==0);save_value('lasso.beta',as.matrix(coef(lasso)))
# Exact small-p active-sign enumeration for Java's stated elastic-net objective,
# including an unpenalized covariate (glmnet rescales supplied penalty factors).
enet_exact <- function(d,lambda,alpha) {
  w<-d$w/sum(d$w);xm<-colSums(d$x*w);ym<-sum(d$y*w)
  z<-sweep(d$x,2,xm);s<-sqrt(colSums(z*z*w));z<-sweep(z,2,s,'/')
  G<-crossprod(z,z*w);r<-drop(crossprod(z,(d$y-ym)*w));pf<-c(0,2)
  H<-G+diag(lambda*(1-alpha)*pf);best<-Inf;answer<-NULL
  signs<-as.matrix(expand.grid(-1:1,-1:1))
  for(i in seq_len(nrow(signs))) {
    sign<-signs[i,];active<-which(sign!=0);beta<-numeric(2)
    if(length(active))beta[active]<-solve(H[active,active,drop=FALSE],r[active]-lambda*alpha*pf[active]*sign[active])
    grad<-drop(H%*%beta-r)
    if(any(beta*sign < -1e-10)||any(abs(grad[sign==0])>lambda*alpha*pf[sign==0]+1e-10))next
    obj<-.5*sum(w*(d$y-ym-z%*%beta)^2)+lambda*(alpha*sum(pf*abs(beta))+.5*(1-alpha)*sum(pf*beta^2))
    if(obj<best){best<-obj;answer<-c(ym-sum(xm*beta/s),beta/s)}
  }
  stopifnot(!is.null(answer));answer
}
for(a in c(0,.35))save_value(paste0('enet.',a),unlist(lapply(lambda,function(l)enet_exact(d,l,a))))
# Weighted held-out risk, independent foldwise glmnet fits. Java emits its
# deterministic seed-17 assignment as a separate reproducible input below.
fold_file <- 'src/test/resources/fitting-audit-v030/folds.txt'
if(file.exists(fold_file)) {
  fold <- scan(fold_file,quiet=TRUE);stopifnot(length(fold)==length(d$y))
  fold_mse <- matrix(0,5,length(lambda));fold_w <- numeric(5)
  for(k in 0:4) {
    train<-fold!=k;test<-!train;f<-glmnet(d$x[train,],d$y[train],weights=d$w[train],alpha=1,lambda=lambda,thresh=1e-14,maxit=100000)
    stopifnot(f$jerr==0);pred<-predict(f,d$x[test,,drop=FALSE]);fold_w[k+1]<-sum(d$w[test])
    fold_mse[k+1,]<-colSums((d$y[test]-pred)^2*d$w[test])/fold_w[k+1]
  }
  mse<-colSums(fold_mse*fold_w)/sum(fold_w)
  se<-sqrt(colSums(sweep(fold_mse,2,mse)^2*fold_w)/sum(fold_w)/4)
  save_value('cv.mse',mse);save_value('cv.se',se)
}
writeLines(paste0(names(values),'=',unlist(values)),out)
cat('Reference fixture:',out,' R=',R.version.string,' glmnet=',as.character(packageVersion('glmnet')),'\n')

n <- if(length(args)>1) as.integer(args[2]) else 3000L
b <- make_data(n)
bench <- function(name,fun,checksum) {
  batch <- 20L
  for(i in 1:3) checksum(fun())
  timings<-numeric(7);total<-0
  for(i in 1:7){start<-proc.time()[['elapsed']];for(j in seq_len(batch)) {f<-fun();total<-total+checksum(f)};timings[i]<-1000*(proc.time()[['elapsed']]-start)/batch}
  cat(sprintf('%s n=%d median_ms=%.6f checksum=%.17g\n',name,n,median(timings),total))
}
bench('ols',function()lm(b$y~b$x,weights=b$w,offset=b$offset),function(f)sum(coef(f)))
bench('poisson',function()glm(b$counts~b$x,weights=b$w,offset=b$offset,family=poisson(),control=glm.control(epsilon=1e-10,maxit=200)),function(f){stopifnot(f$converged);sum(coef(f))})
bench('gee_fixed',function()fixed_gee(b),function(f)sum(f$beta))
bench('lasso_path',function()glmnet(b$x,b$y,weights=b$w,alpha=1,lambda=lambda,thresh=1e-14,maxit=100000),function(f){stopifnot(f$jerr==0);sum(as.matrix(coef(f)))})
