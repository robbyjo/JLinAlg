.libPaths(c('build/r-library', .libPaths()))
suppressPackageStartupMessages(library(lavaan))
options(digits=17)
out <- 'src/test/resources/r-reference/sem-joint'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
set.seed(260908)
n <- 600L
f <- rnorm(n)
g <- .55*f + .8*rnorm(n)
# Non-normal disturbances make robust corrections meaningfully different.
noise <- function(s) s*rt(n,7)/sqrt(7/5)
d <- data.frame(x1=1.2+f+noise(.65), x2=-.3+.8*f+noise(.7),
                x3=.4+1.1*f+noise(.6), y1=.8+g+noise(.55),
                y2=-.7+.9*g+noise(.65), y3=.3+.7*g+noise(.7))
syntax <- 'f =~ 1*x1 + l2*x2 + l3*x3
g =~ 1*y1 + l5*y2 + l6*y3
g ~ b*f
f ~~ vf*f
g ~~ vg*g'
common <- list(model=syntax, data=d, fixed.x=FALSE, auto.cov.y=FALSE,
               information='expected')
fit_case <- function(args) do.call(sem,args)
cases <- list(latent=c(common,list(meanstructure=FALSE)),
              means=c(common,list(meanstructure=TRUE)),
              robust=c(common,list(meanstructure=TRUE, estimator='MLM')))
fixedMean<-common;fixedMean$model<-paste(syntax,'x1 ~ 1.0*1',sep='\n')
cases$fixedMean<-c(fixedMean,list(meanstructure=TRUE))
missing <- d
missing$y2[seq_len(n) %% 5==0 | d$x1>2.2] <- NA
missing$x2[seq_len(n) %% 7==0] <- NA
missing$y1[seq_len(n) %% 9==0] <- NA
missing$x3[seq_len(n) %% 11==0] <- NA
fiml <- common; fiml$data <- missing; fiml$information <- 'observed'
cases$fiml <- c(fiml,list(meanstructure=TRUE,missing='ml'))
med <- data.frame(x=rnorm(n));med$m <- .6+.5*med$x+rnorm(n);med$y <- -.4+.7*med$m+.2*med$x+rnorm(n)
cases$mediation <- list(model='m ~ a*x; y ~ b*m + c*x; indirect := a*b',data=med,
                       meanstructure=TRUE,fixed.x=FALSE,information='expected')
of <- rnorm(n)
latent_response <- cbind(.9*of+rnorm(n), .7*of+rnorm(n),1.1*of+rnorm(n),.8*of+rnorm(n))
ordinal <- as.data.frame(apply(latent_response,2,function(x)as.integer(cut(x,c(-Inf,-.4,.6,Inf)))-1L))
names(ordinal)<-paste0('z',1:4)
cases$ordinal <- list(model='f =~ a*z1 + b*z2 + c*z3 + d*z4',data=ordinal,
                     ordered=names(ordinal),std.lv=TRUE,parameterization='theta',estimator='PML')
write.csv(d,file.path(out,'continuous.csv'),row.names=FALSE,na='NaN')
write.csv(missing,file.path(out,'missing.csv'),row.names=FALSE,na='NaN')
write.csv(med,file.path(out,'mediation.csv'),row.names=FALSE)
write.csv(ordinal,file.path(out,'ordinal.csv'),row.names=FALSE)
emit <- function(con,key,value) writeLines(paste0(key,'=',paste(format(value,digits=17,trim=TRUE,scientific=TRUE),collapse=',')),con)
for(name in names(cases)) {
  args <- cases[[name]]
  fit <- fit_case(args)
  stopifnot(lavInspect(fit,'converged'))
  pt <- parTable(fit); pt <- pt[pt$free>0,];pt <- pt[order(pt$free),]
  labels <- ifelse(nzchar(pt$label),pt$label,paste0(pt$lhs,pt$op,pt$rhs))
  con <- file(file.path(out,paste0(name,'.properties')),'w')
  writeLines(paste0('# R ',getRversion(),'; lavaan ',packageVersion('lavaan')),con)
  writeLines(paste0('labels=',paste(labels,collapse=',')),con)
  emit(con,'estimates',coef(fit));emit(con,'covariance',as.vector(t(vcov(fit))))
  emit(con,'se',sqrt(diag(vcov(fit))))
  if(name!='ordinal') {
    emit(con,'logLikelihood',as.numeric(logLik(fit)))
    if(name=='fiml')for(key in c('chisq','cfi','tli','rmsea'))emit(con,key,fitMeasures(fit,key))
    mi <- modindices(fit); mi <- mi[mi$lhs=='x1' & mi$op=='~~' & mi$rhs=='x2',]
    if(nrow(mi)) {emit(con,'mi.x1.x2',mi$mi);emit(con,'epc.x1.x2',mi$epc)}
    if(name=='robust')emit(con,'scalingFactor',fitMeasures(fit,'chisq.scaling.factor'))
    if(name=='mediation'){pe<-parameterEstimates(fit);emit(con,'indirect.se',pe$se[pe$op==':='])}
    if(name=='means') {
      # Reproducible CR0 with expected information, on lavaan's parameter scale.
      scores <- lavScores(fit)
      groups <- rep(seq_len(n/5),each=5)
      sums <- rowsum(scores,groups)
      bread <- solve(lavInspect(fit,'information')*n)
      emit(con,'clusterCovariance',as.vector(t(bread %*% crossprod(sums) %*% bread)))
    }
  } else {
    # Evaluate the actual sum of pairwise cell log probabilities: optim$fx has
    # lavaan-internal scaling and is not a log likelihood.
    beta<-coef(fit)[1:4];threshold<-matrix(coef(fit)[5:12],nrow=4,byrow=TRUE)
    sd<-sqrt(1+beta^2);threshold<-threshold/sd
    Phi2<-function(a,b,r) {
      if(a==-Inf || b==-Inf)return(0)
      if(a==Inf)return(pnorm(b));if(b==Inf)return(pnorm(a))
      pbivnorm::pbivnorm(a,b,r)
    }
    ll<-0
    for(i in 2:4)for(j in seq_len(i-1)) {
      rho<-beta[i]*beta[j]/(sd[i]*sd[j]);tab<-table(factor(ordinal[[i]],levels=0:2),factor(ordinal[[j]],levels=0:2))
      ti<-c(-Inf,threshold[i,],Inf);tj<-c(-Inf,threshold[j,],Inf)
      for(a in 1:3)for(b in 1:3) {
        prob<-Phi2(ti[a+1],tj[b+1],rho)-Phi2(ti[a],tj[b+1],rho)-Phi2(ti[a+1],tj[b],rho)+Phi2(ti[a],tj[b],rho)
        ll<-ll+tab[a,b]*log(prob)
      }
    }
    emit(con,'pairwiseLogLikelihood',ll)
  }
  elapsed <- replicate(5,system.time(for(batch in 1:5){again<-fit_case(args);stopifnot(lavInspect(again,'converged'))})[['elapsed']]/5)
  emit(con,'R.medianSeconds',median(elapsed))
  close(con)
  cat(name, 'converged; median seconds',median(elapsed),'\n')
}
