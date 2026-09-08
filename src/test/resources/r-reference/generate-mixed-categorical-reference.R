.libPaths(c('build/r-library', .libPaths()))
suppressPackageStartupMessages(library(lme4))
options(digits=17)
out <- 'src/test/resources/r-reference'
levels_f <- c('c','a','b')
design_data <- expand.grid(f=factor(levels_f,levels=levels_f),h=factor(c('v','u'),levels=c('v','u')),x=c(-1,1))
design_data$g <- 'g0'; design_data$y <- seq_len(nrow(design_data))
write.table(design_data,file.path(out,'mixed-categorical-design-data.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
expressions <- c('1+f','0+f','f*x','f:x','x+f:x','f:h','f+h+f:h','0+f+h',
  '0+f*x','0+f:h','0+x+f:x','0+x:f+h','0+x:f+h:x','f+h:f','f:h+x:f:h','x*f*h',
  'f*x+f:h','f:h+f:x+h:x','f*x-1','x:x+f')
records <- list()
for (coding in c('TREATMENT','SUM')) {
  dd <- design_data
  for (name in c('f','h')) {
    contrast <- if (coding=='SUM') contr.sum(nlevels(dd[[name]])) else contr.treatment(nlevels(dd[[name]]))
    colnames(contrast) <- if (coding=='SUM') levels(dd[[name]])[-nlevels(dd[[name]])] else levels(dd[[name]])[-1]
    contrasts(dd[[name]]) <- contrast
  }
  for (expression in expressions) {
    mm <- model.matrix(as.formula(paste('~',expression)),dd)
    for (column in seq_len(ncol(mm))) records[[length(records)+1]] <- data.frame(expression=expression,
      coding=coding,column=colnames(mm)[column],values=paste(format(mm[,column],digits=17,trim=TRUE),collapse=','))
  }
}
write.table(do.call(rbind,records),file.path(out,'mixed-categorical-designs.tsv'),sep='\t',quote=FALSE,row.names=FALSE)

set.seed(81427)
d <- expand.grid(j=0:11,group=0:49)
d$f <- factor(levels_f[d$j%%3+1],levels=levels_f)
d$x <- ifelse(d$j%%6<3,-1,1)
d$g <- sprintf('g%02d',d$group)
d$w <- .7+.2*(d$j%%4); d$o <- .05*sin(d$j)
contrasts(d$f) <- contr.treatment(3)
z <- model.matrix(~f*x,d)
covariance <- diag(c(.9,.5,.7,.4,.3,.35))+.12
b <- matrix(rnorm(50*6),50,6)%*%chol(covariance)
d$y <- 2+.8*d$x+rowSums(z*b[d$group+1,])+d$o+rnorm(nrow(d),sd=.45/sqrt(d$w))
write.table(d[c('y','x','f','g','w','o')],file.path(out,'mixed-categorical-fit-data.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
fits <- list(treatment='1+f',indicators='0+f',sum='1+f',interaction='0+f:x',factorial='1+f*x')
for (name in names(fits)) {
  dd <- d
  contrasts(dd$f) <- if(name=='sum') contr.sum(3) else contr.treatment(3)
  f <- as.formula(paste0('y~x+offset(o)+(',fits[[name]],'|g)'))
  elapsed <- system.time(fit <- lmer(f,dd,weights=w,REML=TRUE,
    control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000,rhoend=1e-9))))['elapsed']
  covariance <- as.vector(t(as.matrix(VarCorr(fit)[[1]])))
  modes <- as.numeric(ranef(fit)[[1]][1,])
  values <- c(fixef(fit),sigma(fit)^2,as.numeric(logLik(fit)),covariance,modes,head(fitted(fit),12),elapsed)
  names(values) <- c('beta0','beta1','residual','logLik',paste0('cov',seq_along(covariance)-1),
    paste0('mode',seq_along(modes)-1),paste0('fitted',0:11),'fitSeconds')
  writeLines(c(paste0('# R ',getRversion(),'; lme4 ',packageVersion('lme4'),'; random=',fits[[name]]),
    paste(names(values),format(values,digits=17,scientific=TRUE),sep='=')),
    file.path(out,paste0('mixed-categorical-',name,'.properties')))
  print(c(case=name,logLik=as.numeric(logLik(fit)),singular=isSingular(fit),seconds=elapsed))
}
