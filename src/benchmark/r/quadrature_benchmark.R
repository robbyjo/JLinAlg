.libPaths(c('build/r-library', .libPaths()))
suppressPackageStartupMessages(library(lme4))
# Same four datasets, full ML fit plus fixed-effect inference, 3 warmups then
# median of 7 timed fits. Parsing and accuracy checks are outside the timer.
root <- 'src/test/resources/r-reference'
reference <- readLines(file.path(root,'quadrature-fits.properties'))
getref <- function(kind,key) as.numeric(sub('^[^=]*=','',reference[startsWith(reference,paste0(kind,'.',key,'='))]))
cat('family\tmedianSeconds\tchecksum\n')
for(kind in c('binary','rare','binomial','poisson')) {
  d <- read.delim(file.path(root,paste0('quadrature-',kind,'.tsv')))
  d$group <- factor(d$group)
  family <- if(kind=='poisson') poisson() else binomial()
  times <- numeric(7); checksum <- 0
  for(iteration in -2:7) {
    seconds <- system.time({
      fit <- glmer(y~x+offset(offset)+(1|group),data=d,family=family,weights=trials,nAGQ=25,
        control=glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=200000)))
      se <- sqrt(diag(vcov(fit)))
    })[['elapsed']]
    stopifnot(is.null(fit@optinfo$conv$lme4$messages),fit@optinfo$conv$opt==0,
      max(abs(fixef(fit)-c(getref(kind,'intercept'),getref(kind,'slope'))))<1e-5,
      max(abs(se-c(getref(kind,'se0'),getref(kind,'se1'))))<1e-5,
      abs(sqrt(as.numeric(VarCorr(fit)$group))-getref(kind,'sd'))<1e-5)
    if(iteration>0) {
      times[iteration] <- seconds
      # Correct the count-specific constant in lme4's reported nAGQ>1 LL.
      full <- as.numeric(logLik(fit))-getref(kind,'lme4ReportedLogLik')+getref(kind,'logLikelihood')
      stopifnot(abs(full-getref(kind,'logLikelihood'))<2e-6)
      checksum <- checksum+full+fixef(fit)[1]+se[1]
    }
  }
  cat(sprintf('%s\t%.9f\t%.12f\n',kind,median(times),checksum))
}
