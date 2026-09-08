.libPaths(c('build/r-library',.libPaths()))
suppressPackageStartupMessages(library(lme4))
root='src/test/resources/r-reference/remaining-mixed'
raw=readLines(file.path(root,'reference.properties'));parts=strsplit(raw[!startsWith(raw,'R.version=')],'=',fixed=TRUE)
references=setNames(lapply(parts,function(p)as.numeric(strsplit(p[2],',',fixed=TRUE)[[1]])),vapply(parts,`[`,'',1))
check=function(actual,want,tolerance)stopifnot(all(is.finite(actual)),max(abs(actual-want))<tolerance)
cat('engine\tcase\trun\tseconds\tchecksum\tLLerror\n')
for(kind in c('rare','poisson','REML','ML')) {
 gaussian=kind%in%c('REML','ML');d=read.delim(file.path(root,paste0(if(gaussian)'reml' else kind,'.tsv')))
 prefix=if(gaussian)paste0('reml.',kind) else kind;times=numeric(7);checksum=0
 for(run in -2:7) {
  batch=if(gaussian)20L else 1L
  time=system.time({
   fits=lapply(seq_len(batch),function(unused) {
    fit=if(gaussian)lmer(y~x+(1|group),d,REML=kind=='REML') else
    glmer(y~x+offset(offset)+(1|group),d,family=if(kind=='rare')binomial() else poisson(),nAGQ=1,
     control=glmerControl(optimizer='bobyqa',tolPwrss=1e-12,optCtrl=list(rhoend=1e-10)))
    list(fit=fit,covariance=as.matrix(vcov(fit)))
   })
  })[['elapsed']]/batch
  for(result in fits) {
  fit=result$fit;covariance=result$covariance
  stopifnot(is.null(fit@optinfo$conv$lme4$messages),fit@optinfo$conv$opt==0)
  check(fixef(fit),references[[paste0(prefix,'.beta')]],2e-6)
  check(as.data.frame(VarCorr(fit))$vcov,references[[paste0(prefix,'.variance')]],3e-6)
  check(as.vector(t(covariance)),references[[paste0(prefix,'.covariance')]],3e-6)
  error=abs(as.numeric(logLik(fit))-references[[paste0(prefix,'.LL')]]);stopifnot(error<2e-8)
  checksum=checksum+as.numeric(logLik(fit))+sum(fixef(fit))+sum(sqrt(diag(covariance)))
  }
  if(run>0){times[run]=time;cat(sprintf('R-lme4\t%s\t%d\t%.9f\t%.14g\t%.3g\n',kind,run,time,checksum,error))}
 }
 cat(sprintf('R-lme4\t%s\tmedian\t%.9f\t%.14g\t0\n',kind,median(times),checksum))
}
