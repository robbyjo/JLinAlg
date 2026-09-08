.libPaths(c('build/r-library', .libPaths()))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(lme4))
options(digits=17)
root <- 'src/test/resources/remaining-model-audit'
dir.create(root,recursive=TRUE,showWarnings=FALSE)
refs <- data.frame(model=character(),quantity=character(),index=integer(),value=double())
save <- function(model,quantity,value) {
 refs <<- rbind(refs,data.frame(model=model,quantity=quantity,index=seq_along(value)-1L,value=as.numeric(value)))
}
write_data <- function(name,x) write.table(x,file.path(root,paste0(name,'.csv')),sep=',',quote=FALSE,row.names=FALSE,col.names=FALSE)
timings <- data.frame(model=character(),repeats=integer(),ms=double(),checksum=double())
time_fit <- function(name,run,repeats=50L) {
 for(k in 1:10)run()
 checksum<-0;started<-Sys.time()
 for(k in seq_len(repeats))checksum<-checksum+sum(run())
 elapsed<-as.numeric(difftime(Sys.time(),started,units='secs'))*1000/repeats
 timings <<- rbind(timings,data.frame(model=name,repeats=repeats,ms=elapsed,checksum=checksum))
}

i<-0:239;x<-sin(i*.71)+.2*cos(i*.13);z<-cos(i*.39)
start<-(i%%11)*.17;stop<-ceiling((start+1+(i%%31)*.23)*2)/2
event<-i%%4!=0;stratum<-i%%3;offset<-.2*sin(i*.29)
write_data('cox',cbind(start,stop,as.integer(event),stratum,x,z,offset))
for(ties in c('breslow','efron')) {
 run<-function(){f<-coxph(Surv(start,stop,event)~x+z+offset(offset)+strata(stratum),ties=ties,control=coxph.control(eps=1e-10,iter.max=50));stopifnot(f$iter<50);f}
 fit<-run();name<-paste0('cox_',ties)
 save(name,'beta',coef(fit));save(name,'se',sqrt(diag(vcov(fit))));save(name,'ll',fit$loglik[2])
 time_fit(name,function()coef(run()))
}

i<-0:35;time<-i+1;event<-i%%4!=0;x<-sin(i*1.3);g<-i%%4;theta<-.4
write_data('gamma_frailty',cbind(time,as.integer(event),x,g))
partial<-function(v){eta<-v[1]*x+v[2:5][g+1];sum(sapply(which(event),function(i){r<-eta[time>=time[i]];a<-max(r);eta[i]-a-log(sum(exp(r-a)))}))}
gamma_objective<-function(v)-partial(v)-sum(dgamma(exp(v[2:5]),shape=1/theta,scale=theta,log=TRUE)+v[2:5])
gammafit<-function()optim(rep(0,5),gamma_objective,method='BFGS',control=list(reltol=1e-12,maxit=1000),hessian=TRUE)
fit<-gammafit();stopifnot(fit$convergence==0)
save('gamma_frailty','beta',fit$par[1]);save('gamma_frailty','mode',fit$par[2:5]);save('gamma_frailty','partial',partial(fit$par))
save('gamma_frailty','laplace',-fit$value+2*log(2*pi)-sum(log(diag(chol(fit$hessian[2:5,2:5])))))
time_fit('gamma_frailty',function(){f<-gammafit();stopifnot(f$convergence==0);f$par},20L)

i<-0:119;t<-(i%%13-6)/3;c<-cos(i*.23);g<-i%/%10
m<-1+.7*t-.4*c+1.3*sin(g*1.7)+.5*sin(i*2.3)
y<-2+.3*t+.8*m+.2*c+.7*cos(g*2.1)+.4*cos(i*1.3)
write_data('mediation',cbind(y,t,m,c,g))
med<-function(mixed=FALSE){
 if(mixed){a<-lmer(m~t+c+(1|g),REML=TRUE);b<-lmer(y~t+m+c+(1|g),REML=TRUE);total<-lmer(y~t+c+(1|g),REML=TRUE);coefn<-fixef
 stopifnot(all(vapply(list(a,b,total),function(f)is.null(f@optinfo$conv$lme4$messages)&&f@optinfo$conv$opt==0,logical(1))))}
 else {a<-lm(m~t+c);b<-lm(y~t+m+c);total<-lm(y~t+c);coefn<-coef}
 av<-coefn(a)[2];bv<-coefn(b)[3];se<-sqrt(bv^2*vcov(a)[2,2]+av^2*vcov(b)[3,3])
 c(a=av,b=bv,direct=coefn(b)[2],total=coefn(total)[2],indirect=av*bv,se=se)
}
save('mediation','effects',med());save('mediation_mixed','effects',med(TRUE))
time_fit('mediation',function()med());time_fit('mediation_mixed',function()med(TRUE),20L)

i<-0:239;x<-i/120;y<-1.7*exp(.35*x)+.15*sin(i*1.7)+.1*cos(i*.31)
write_data('nonlinear',cbind(y,x))
nlsfit<-function()nls(y~a*exp(b*x),start=list(a=1,b=.1),control=nls.control(tol=1e-9,maxiter=100))
fit<-nlsfit();save('nonlinear','beta',coef(fit));save('nonlinear','se',sqrt(diag(vcov(fit))));save('nonlinear','sse',deviance(fit))
time_fit('nonlinear',function()coef(nlsfit()))

i<-0:47;x<-i/24;id<-i%%6
A<-diag(6);A[1,c(3,4)]<-.5;A[2,c(3,4)]<-.5;A[3,4]<-.5;A[lower.tri(A)]<-t(A)[lower.tri(A)]
Z<-diag(6)[id+1,,drop=FALSE];K<-Z%*%A%*%t(Z)
y<-1.7*exp(.35*x)+c(.6,-.5,.3,-.2,.8,-.7)[id+1]+.2*sin(i*1.7)
write_data('nonlinear_pedigree',cbind(y,x,id))
objective<-function(v){
 V<-exp(v[3])*K+diag(exp(v[4]),length(y));R<-chol(V)
 e<-y-v[1]*exp(v[2]*x);u<-backsolve(R,e,transpose=TRUE)
 .5*(length(y)*log(2*pi)+2*sum(log(diag(R)))+sum(u*u))
}
pedfit<-function()optim(c(1.7,.35,log(.3),log(.03)),objective,method='L-BFGS-B',lower=c(.1,-2,log(1e-8),log(1e-8)),upper=c(5,2,log(100),log(100)),control=list(factr=1e3,pgtol=1e-7,ndeps=rep(1e-5,4),maxit=1000))
fit<-pedfit();stopifnot(fit$convergence==0)
save('nonlinear_pedigree','beta',fit$par[1:2]);save('nonlinear_pedigree','variance',exp(fit$par[3:4]));save('nonlinear_pedigree','ll',-fit$value)
time_fit('nonlinear_pedigree',function(){f<-pedfit();stopifnot(f$convergence==0);f$par},20L)
write.table(refs,file.path(root,'reference.csv'),sep=',',row.names=FALSE,quote=FALSE)
print(timings,row.names=FALSE)
evidence <- 'src/benchmark/resources/remaining-model-audit-benchmark'
dir.create(evidence,recursive=TRUE,showWarnings=FALSE)
write.table(timings,file.path(evidence,'measured-r.csv'),sep=',',row.names=FALSE,quote=FALSE)
cat('R',as.character(getRversion()),'survival',as.character(packageVersion('survival')),'lme4',as.character(packageVersion('lme4')),'\n')
