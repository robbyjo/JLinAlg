# Reproduce independent fixtures with R stats, lavaan, and direct likelihoods.
.libPaths(c('build/r-library', .libPaths()))
library(lavaan)
out <- 'src/test/resources/estimator-extensions'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
write_matrix <- function(x, name) write.table(x, file.path(out,name), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE, na='NaN')
options(digits=17)
t <- 1:80
x <- cbind(1,sin(t*.35))
e <- sin(t*1.37)+.3*cos(t*.71)
for(i in 2:80) e[i] <- .4*e[i-1]+e[i]
y <- 1.2+.7*x[,2]+e
y[c(4,17,18,55)] <- NA
f <- arima(y,order=c(1,0,0),xreg=x,include.mean=FALSE,method='ML',optim.control=list(reltol=1e-12))
write_matrix(cbind(y,x),'arima-data.tsv')
write_matrix(c(coef(f),f$sigma2,as.numeric(logLik(f))),'arima-reference.tsv')
sm <- KalmanSmooth(y,makeARIMA(.4,.2,numeric()))
write_matrix(sm$smooth,'smooth-stationary.tsv')
sm <- KalmanSmooth(y,makeARIMA(.4,.2,1))
write_matrix(sm$smooth,'smooth-integrated.tsv')

set.seed(81703)
n <- 600
latent <- rnorm(n)
z <- sapply(c(.55,.8,.7,.9),function(a) a*latent+rnorm(n))
dat <- as.data.frame(apply(z,2,function(v) as.integer(cut(v,c(-Inf,-.4,.5,Inf)))-1))
names(dat) <- paste0('y',1:4)
f <- cfa('F =~ y1+y2+y3+y4',data=dat,ordered=names(dat),std.lv=TRUE,parameterization='theta',estimator='WLSMV')
stopifnot(lavInspect(f,'converged'))
obs <- lavInspect(f,'wls.obs'); gamma <- lavInspect(f,'gamma')
keys <- c(unlist(lapply(names(dat),function(s)paste0(s,'|t',1:2))), 'y1~~y2','y1~~y3','y2~~y3','y1~~y4','y2~~y4','y3~~y4')
ix <- match(keys,names(obs));stopifnot(!anyNA(ix))
write_matrix(obs[ix],'dwls-moments.tsv');write_matrix(gamma[ix,ix],'dwls-gamma.tsv')
test <- lavInspect(f,'test')$scaled.shifted
write_matrix(c(coef(f)[1:4],test$stat,test$scaling.factor,test$shift.parameter),'dwls-reference.tsv')
write_matrix(vcov(f)[1:4,1:4],'dwls-covariance.tsv')
write_matrix(as.matrix(dat),'ordinal-data.tsv')

# Direct mixed likelihood: continuous X, binary Y, MAR missingness depending on X.
set.seed(820)
n <- 140
x <- rnorm(n,.4,1.1);latent <- .65*x+rnorm(n)
y <- as.numeric(latent>.2)
y[x>.7 & (1:n)%%3==0] <- NA
x[(1:n)%%11==0] <- NA
dat <- cbind(x,y)
nll <- function(p) {
    mu<-p[1];sd<-exp(p[2]);b<-p[3];threshold<-p[4]
    sum(vapply(1:n,function(i) {
        if(is.na(x[i])) {
            if(is.na(y[i]))return(0)
            prob<-pnorm((threshold-b*mu)/sqrt(1+b*b*sd*sd))
            return(-log(if(y[i]==0)prob else 1-prob))
        }
        val<- -dnorm(x[i],mu,sd,log=TRUE)
        if(!is.na(y[i]))val<-val-pnorm(threshold-b*x[i],lower.tail=y[i]==0,log.p=TRUE)
        val
    },numeric(1)))
}
fit<-optim(c(.3,0,.4,0),nll,method='BFGS',control=list(reltol=1e-12),hessian=TRUE)
stopifnot(fit$convergence==0)
write_matrix(dat,'mixed-data.tsv')
write_matrix(c(fit$par,-fit$value),'mixed-reference.tsv')
write_matrix(solve(fit$hessian),'mixed-covariance.tsv')

# Independently construct observed-information LR scaling in saturated moment
# coordinates. Center case scores at the fitted structured model, matching the
# explicitly documented convention; this is not the Mplus YB shortcut.
library(numDeriv)
d <- read.csv('src/test/resources/r-reference/sem-joint/missing.csv')
model <- 'f =~ 1*x1+l2*x2+l3*x3
          g =~ 1*y1+l5*y2+l6*y3
          g ~ b*f
          f ~~ vf*f
          g ~~ vg*g'
f <- sem(model,data=d,missing='ml',meanstructure=TRUE,information='observed')
stopifnot(lavInspect(f,'converged'))
implied<-fitted(f);mu<-as.vector(implied$mean);sigma<-implied$cov
p<-length(mu);lower<-which(lower.tri(sigma,diag=TRUE),arr.ind=TRUE)
s<-c(mu,sigma[lower]);n<-nrow(d);d<-as.matrix(d)
scores<-function(moment) {
    m<-moment[1:p];v<-matrix(0,p,p)
    v[lower]<-moment[-(1:p)];v<-v+t(v)-diag(diag(v))
    vapply(1:n,function(i){keep<-which(is.finite(d[i,]));if(!length(keep))return(0)
        delta<-d[i,keep]-m[keep];cv<-v[keep,keep,drop=FALSE]
        -.5*(length(keep)*log(2*pi)+as.numeric(determinant(cv,logarithm=TRUE)$modulus)+sum(delta*solve(cv,delta)))
    },numeric(1))
}
A1<-hessian(function(z)-sum(scores(z)),s)
B1<-jacobian(scores,s);B1<-scale(B1,center=TRUE,scale=FALSE)
delta<-lavInspect(f,'delta')
# lavaan moment order is means then column-wise lower covariance triangle.
A0inv<-solve(lavInspect(f,'information')*n)
U<-solve(A1)-delta%*%A0inv%*%t(delta)
scaling<-sum((B1%*%U)*B1)/fitMeasures(f,'df')
write_matrix(c(fitMeasures(f,'chisq'),scaling,fitMeasures(f,'chisq')/scaling),'fiml-robust-reference.tsv')

# Local linear mixed-kernel HC3 from explicit base-R weighted matrix algebra.
i<-0:119;x<-(i%%60)/59;group<-i%/%60;y<-1+2*x+3*group+.3*sin(i*1.73)
keep<-group==1;xx<-x[keep];yy<-y[keep];h<-min(sd(xx),diff(quantile(xx,c(.25,.75)))/1.349)*length(xx)^(-1/3)
weights<-function(query) {
    w<-exp(-.5*((xx-query)/h)^2);design<-cbind(1,(xx-query)/h)
    drop(solve(crossprod(design,w*design))[1,,drop=FALSE]%*%t(w*design))
}
influence<-weights(.5);estimate<-sum(influence*yy);variance<-0
for(j in seq_along(xx)){w<-weights(xx[j]);variance<-variance+(influence[j]*(yy[j]-sum(w*yy))/(1-w[j]))^2}
write_matrix(cbind(x,group,y),'kernel-data.tsv')
write_matrix(c(estimate,sqrt(variance),h),'kernel-reference.tsv')
writeLines(trimws(capture.output(sessionInfo()),which='right'),file.path(out,'R-session.txt'))
