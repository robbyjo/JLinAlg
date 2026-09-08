.libPaths(c('build/r-library', .libPaths()))
suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(clubSandwich))
options(digits=17)
out <- 'src/test/resources/meta'
dir.create(out, showWarnings=FALSE, recursive=TRUE)
reference <- data.frame(key=character(), value=double())
save_values <- function(key, value) {
  reference <<- rbind(reference, data.frame(key=paste0(key, '.', seq_along(value)-1L), value=as.numeric(value)))
}
save_fit <- function(key, fit) {
  save_values(paste0(key,'.beta'), coef(fit))
  save_values(paste0(key,'.cov'), as.vector(t(vcov(fit))))
  save_values(paste0(key,'.p'), fit$pval)
  save_values(paste0(key,'.ll'), logLik(fit))
}

bias <- data.frame(yi=c(.1,.12,.08,.15,.2,.28,.34,.45,.55,.7,.9,1.1),
                   sei=c(.05,.07,.09,.11,.14,.18,.22,.28,.35,.44,.55,.7))
write.table(bias, file.path(out,'bias.csv'), sep=',', row.names=FALSE, quote=FALSE)
for (method in c('FE','REML','DL','PM')) {
  f <- rma(yi, sei=sei, data=bias, method=method, control=list(threshold=1e-12,tol=1e-12,maxiter=1000))
  for (estimator in c('L0','R0','Q0')) for (side in c('left','right')) {
    key <- paste('trim',method,estimator,side,sep='.')
    tf <- tryCatch(suppressWarnings(trimfill(f, estimator=estimator, side=side)),error=function(e) NULL)
    if (is.null(tf)) { save_values(paste0(key,'.undefined'),1); next }
    # trimfill's final augmented-data refit drops the input fit's control list.
    # Refit the identical augmented data at the declared comparison tolerance.
    refined <- rma(tf$yi,tf$vi,method=method,control=list(threshold=1e-12,tol=1e-12,maxiter=1000))
    save_values(paste0(key,'.k0'), tf$k0)
    save_values(paste0(key,'.beta'), coef(refined))
    save_values(paste0(key,'.se'), refined$se)
  }
  tf <- trimfill(f)
  refined <- rma(tf$yi,tf$vi,method=method,control=list(threshold=1e-12,tol=1e-12,maxiter=1000))
  save_values(paste0('trim.',method,'.auto.beta'),coef(refined))
  save_values(paste0('trim.',method,'.auto.k0'),tf$k0)
}
pet <- lm(yi ~ sei, weights=1/sei^2, data=bias)
peese <- lm(yi ~ I(sei^2), weights=1/sei^2, data=bias)
for (key in c('pet','peese')) {
  fit <- get(key)
  save_values(paste0(key,'.beta'),coef(fit))
  save_values(paste0(key,'.cov'),as.vector(t(vcov(fit))))
  save_values(paste0(key,'.p'),coef(summary(fit))[,4])
  save_values(paste0(key,'.ci'),as.vector(t(confint(fit,level=.90))))
}
eg <- regtest(rma(yi,sei=sei,data=bias,method='FE'),model='lm',predictor='sei',ret.fit=TRUE)
save_values('egger', c(eg$zval,eg$pval,coef(eg$fit)[2]))

j <- 1:20; moderate <- data.frame(yi=.2+.45*sin(j*2.1)+.6*(.08+.025*j),sei=.08+.025*j)
write.table(moderate,file.path(out,'trim-moderate.csv'),sep=',',row.names=FALSE,quote=FALSE)
for (estimator in c('L0','R0','Q0')) {
  f <- rma(yi,sei=sei,data=moderate,method='REML',control=list(threshold=1e-12,tol=1e-12,maxiter=1000))
  tf <- trimfill(f,side='left',estimator=estimator)
  refined <- rma(tf$yi,tf$vi,method='REML',control=list(threshold=1e-12,tol=1e-12,maxiter=1000))
  save_values(paste0('moderate.',estimator),c(tf$k0,coef(refined),refined$se,refined$tau2))
}

g <- rep(1:12,each=3); x <- rep(c(-1,0,1),12); i <- seq_along(g)
vi <- .035 + .009*(i %% 4)
y <- .4 + .3*x + .7*sin(g*1.7) + .45*cos(g*.9)*x + .22*sin(i*2.3)
dat <- data.frame(yi=y, vi=vi, x=x, cluster=g, effect=i)
V <- diag(vi)
for (a in i) for (b in i) if (g[a]==g[b] && a!=b) V[a,b] <- .25*sqrt(vi[a]*vi[b])
write.table(dat,file.path(out,'hierarchical.csv'),sep=',',row.names=FALSE,quote=FALSE)
write.table(V,file.path(out,'sampling.csv'),sep=',',row.names=FALSE,col.names=FALSE,quote=FALSE)
for (method in c('ML','REML')) {
  cat('Generating random coefficients',method,'\n')
  fit <- rma.mv(yi,V,mods=~x,random=~1+x|cluster,struct='GEN',data=dat,method=method,
                control=list(rel.tol=1e-10,iter.max=2000,eval.max=4000))
  save_fit(paste0('hier.',method),fit)
  save_values(paste0('hier.',method,'.G'),as.vector(t(fit$G)))
  cat('Generating nested intercepts',method,'\n')
  nested <- rma.mv(yi,V,mods=~x,random=list(~1|cluster,~1|effect),data=dat,method=method,
                   control=list(rel.tol=1e-10,iter.max=2000,eval.max=4000))
  save_fit(paste0('nested.',method),nested)
  save_values(paste0('nested.',method,'.variance'),nested$sigma2)
  if (method=='REML') {
    diagonal <- rma.mv(yi,V,mods=~x,random=~1+x|cluster,struct='GDIAG',data=dat,method=method,
                      control=list(rel.tol=1e-10,iter.max=2000,eval.max=4000))
    save_fit('diagonal.REML',diagonal)
    save_values('diagonal.REML.G',as.vector(t(diagonal$G)))
    robust_hier <- vcovCR(fit,cluster=dat$cluster,type='CR2')
    save_values('hier.CR2.cov',as.vector(t(robust_hier)))
    save_values('hier.CR2.df',coef_test(fit,vcov=robust_hier)$df_Satt)
    save_values('hier.CR2.p',coef_test(fit,vcov=robust_hier)$p_Satt)
  }
}
known <- rma.mv(yi,V,mods=~x,data=dat,method='ML')
save_fit('known',known)
for (correction in c('CR0','CR1','CR1p','CR1S','CR2')) {
  cr <- vcovCR(known,cluster=dat$cluster,type=correction)
  save_values(paste0(correction,'.cov'),as.vector(t(cr)))
  test <- coef_test(known,vcov=cr,test=if (correction=='CR2') 'Satterthwaite' else 'naive-tp')
  save_values(paste0(correction,'.p'),test[[ncol(test)]])
  if (correction=='CR2') save_values('CR2.df',test$df_Satt)
}
# rma.mv uses REML; DL and PM below are independently evaluated generalized estimating equations.
X <- matrix(1,nrow(dat),1)
gls <- function(t) {
  W <- solve(V+diag(t,nrow(V))); B <- solve(t(X)%*%W%*%X)
  beta <- B%*%t(X)%*%W%*%y; e <- y-X%*%beta
  P <- W-W%*%X%*%B%*%t(X)%*%W
  list(beta=beta,cov=B,q=drop(t(e)%*%W%*%e),P=P)
}
q0 <- gls(0)
for (method in c('REML','DL','PM')) {
  tau <- switch(method,
    REML=rma.mv(yi,V,random=~1|effect,data=dat,method='REML',control=list(rel.tol=1e-12))$sigma2,
    DL=max(0,(q0$q-nrow(dat)+1)/sum(diag(q0$P))),
    PM=uniroot(function(t) gls(t)$q-nrow(dat)+1,c(0,10),tol=1e-12)$root)
  f <- gls(tau)
  save_values(paste0('correlated.',method),c(tau,f$beta,f$cov,f$q,q0$q))
}
write.table(reference,file.path(out,'metafor-reference.csv'),sep=',',row.names=FALSE,quote=FALSE)
cat('Generated',nrow(reference),'reference scalars with metafor',as.character(packageVersion('metafor')),
    'and clubSandwich',as.character(packageVersion('clubSandwich')),'\n')
