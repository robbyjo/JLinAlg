.libPaths(c('build/r-library', .libPaths()))
suppressPackageStartupMessages(library(glmmTMB))
options(digits=17)
root <- 'src/test/resources/r-reference'
out <- 'src/benchmark/resources/distributional-audit'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
b <- read.delim(file.path(root,'beta-mixed-glmmtmb.tsv')); b$group <- factor(b$group)
z <- read.delim(file.path(root,'zero-inflated-mixed-data.tsv')); z$group <- factor(z$group)
runs <- list(
  betaMixed=function() glmmTMB(y~x+(1|group),data=b,family=beta_family()),
  zip=function() glmmTMB(y_zip~x+(1|group),ziformula=~x+(1|group),data=z,family=poisson),
  zinb=function() glmmTMB(y_zinb~x+(1|group),ziformula=~x,data=z,family=nbinom2))
values <- c()
for(name in names(runs)) {
  run <- runs[[name]]
  for(i in 1:2) {f<-run();stopifnot(f$fit$convergence==0,f$sdr$pdHess)}
  times <- numeric(5); checksum<-0
  for(i in 1:5) {
    times[i]<-system.time(f<-run())[['elapsed']]
    stopifnot(f$fit$convergence==0,f$sdr$pdHess)
    checksum<-checksum+as.numeric(logLik(f))+sum(fixef(f)$cond)+sum(fixef(f)$zi)
  }
  v<-c(seconds=median(times),logLik=as.numeric(logLik(f)),checksum=checksum)
  values<-c(values,setNames(v,paste(name,names(v),sep='.')))
  cat(name,'seconds=',paste(times,collapse=','),'logLik=',as.numeric(logLik(f)),'checksum=',checksum,'\n')
}
writeLines(c(paste0('# R ',getRversion(),'; glmmTMB ',packageVersion('glmmTMB')),
  paste(names(values),format(values,digits=17,scientific=TRUE),sep='=')),file.path(out,'r-timings.properties'))
