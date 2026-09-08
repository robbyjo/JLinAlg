.libPaths(c('build/r-library', .libPaths()))
suppressPackageStartupMessages({ library(lme4); library(lmerTest); library(pbkrtest) })
options(digits=17)
out <- 'src/test/resources/r-reference'
save_case <- function(name, d, multiple=FALSE, profile_ci=FALSE) {
  write.table(d, file.path(out,paste0('mixed-',name,'.tsv')), sep='\t', row.names=FALSE, quote=FALSE)
  two <- name == 'two-correlated'
  f <- if (two) y ~ x + offset(o) + (1+x|g) + (1+x|h) else if (multiple) y ~ x + offset(o) + (1+x|g) + (1|h) else y ~ x + (1+x|g)
  ctrl <- lmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000, rhobeg=.1, rhoend=1e-9))
  fit <- lmerTest::lmer(f,d,weights=w,REML=TRUE,control=ctrl)
  elapsed <- replicate(3,system.time(lmerTest::lmer(f,d,weights=w,REML=TRUE,control=ctrl))['elapsed'])
  s <- coef(summary(fit,ddf='Satterthwaite'))
  # pbkrtest's residual basis is I; explicit whitening preserves the weighted
  # estimand and permits a correct KR comparison even for nonconstant weights.
  dd <- transform(d, ys=(y-o)*sqrt(w), ws=sqrt(w), xs=x*sqrt(w))
  fk <- if (two) ys~0+ws+xs+(0+ws+xs|g)+(0+ws+xs|h) else if (multiple) ys~0+ws+xs+(0+ws|h)+(0+ws+xs|g) else ys~0+ws+xs+(0+ws+xs|g)
  white <- lmerTest::lmer(fk,dd,REML=TRUE,control=ctrl)
  krtime <- system.time(kr <- coef(summary(white,ddf='Kenward-Roger')))['elapsed']
  values <- c(beta=fixef(fit), covariance=as.vector(t(as.matrix(VarCorr(fit)[[1]]))),
    residual=sigma(fit)^2, logLik=as.numeric(logLik(fit)), se=s[,2], df=s[,3],
    krse=kr[,2], krdf=kr[,3], fitSeconds=median(elapsed), krSeconds=krtime)
  if (multiple) values <- c(values, secondVariance=as.matrix(VarCorr(fit)[[2]])[1,1])
  if (two) values <- c(values, secondCovariance=as.matrix(VarCorr(fit)[[2]])[1,2],secondSlopeVariance=as.matrix(VarCorr(fit)[[2]])[2,2])
  if (profile_ci) {
    ci <- confint(fit,parm='x',method='profile',quiet=TRUE)
    values <- c(values, profileLower=ci[1,1],profileUpper=ci[1,2])
    vc <- confint(fit,parm='theta_',method='profile',quiet=TRUE,oldNames=FALSE)
    intercept <- grep('^sd_.*Intercept',rownames(vc)); slope <- grep('^sd_x',rownames(vc))
    correlation <- grep('^cor_',rownames(vc)); residual <- grep('sigma',rownames(vc))
    values <- c(values, randomSdLower=vc[intercept,1],randomSdUpper=vc[intercept,2],
      slopeSdLower=vc[slope,1],slopeSdUpper=vc[slope,2],
      correlationLower=vc[correlation,1],correlationUpper=vc[correlation,2],
      residualSdLower=vc[residual,1],residualSdUpper=vc[residual,2])
  }
  # Stable names independent of R coefficient-name concatenation.
  names(values)[1:14] <- c('beta0','beta1','cov00','cov01','cov10','cov11','residual',
    'logLik','se0','se1','df0','df1','krse0','krse1')
  names(values)[15:18] <- c('krdf0','krdf1','fitSeconds','krSeconds')
  writeLines(c(paste0('# R ',getRversion(),'; lme4 ',packageVersion('lme4'),
    '; lmerTest ',packageVersion('lmerTest'),'; pbkrtest ',packageVersion('pbkrtest')),
    paste(names(values),format(values,digits=17,scientific=TRUE),sep='=')),
    file.path(out,paste0('mixed-',name,'.properties')))
  print(c(case=name, values))
}
d <- data.frame(y=sleepstudy$Reaction,x=sleepstudy$Days,g=as.character(sleepstudy$Subject),h='h0',w=1,o=0)
save_case('sleepstudy',d,profile_ci=TRUE)
save_case('unbalanced',d[-c(2:5,13,21:26,41,53:57,81:83,95,112:115,143:145,162:165),])
set.seed(81724)
d <- expand.grid(j=0:7,group=0:29)
d$x <- (d$j-3.5)/2; d$g <- sprintf('g%02d',d$group)
d$h <- sprintf('h%02d',(d$group*3+d$j)%%12)
b <- matrix(rnorm(60),30,2) %*% chol(matrix(c(1.4,.38,.38,.7),2))
h <- rnorm(12,sd=.7)
d$w <- .6+(d$j%%4)*.5; d$o <- .2*sin(d$j)
d$y <- 2+.8*d$x+b[d$group+1,1]+b[d$group+1,2]*d$x+h[as.integer(sub('h','',d$h))+1]+d$o+rnorm(nrow(d),sd=.6/sqrt(d$w))
save_case('weighted-crossed',d[c('y','x','g','h','w','o')],multiple=TRUE)
d$y <- d$y + .5*sin(as.integer(sub('h','',d$h))*1.13)*d$x
save_case('two-correlated',d[c('y','x','g','h','w','o')],multiple=TRUE)
