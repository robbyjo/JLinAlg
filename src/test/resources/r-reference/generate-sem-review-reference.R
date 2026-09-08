.libPaths(c('build/r-library',.libPaths()))
suppressPackageStartupMessages(library(lavaan))
options(digits=17)
destination<-'src/test/resources/r-reference/sem-review.properties'
con<-file(destination,'w')
emit<-function(key,value)writeLines(paste0(key,'=',paste(format(value,digits=17,scientific=TRUE,trim=TRUE),collapse=',')),con)
writeLines(paste0('# R ',getRversion(),'; lavaan ',packageVersion('lavaan'),'; analytic/base integrate checks'),con)
rho<-.9998;p<-.25+asin(rho)/(2*pi)
emit('orthant',p);emit('orthantDerivative',1/(2*pi*sqrt(1-rho^2)))
emit('nearBoundaryLL',9936*log(p)+64*log(.5-p))
normalInterval<-function(a,b)ifelse(a>=0,pnorm(a,lower.tail=FALSE)-pnorm(b,lower.tail=FALSE),pnorm(b)-pnorm(a))
rectangle<-function(a,b,c,d,rho) {
  # Positive conditional-normal integrand, independent of the Java Plackett CDF.
  integrate(function(x)dnorm(x)*normalInterval((c-rho*x)/sqrt(1-rho*rho),(d-rho*x)/sqrt(1-rho*rho)),
    max(-40,a),min(40,b),abs.tol=1e-300,rel.tol=2e-11,subdivisions=2000L)$value
}
cases<-rbind(c(4,5,-5,-4,.9),c(2.3,Inf,-Inf,-2.3,.95),
             c(-Inf,.25,-Inf,.3,.9998),c(-Inf,.25,-Inf,-.3,-.9998),
             c(-Inf,-4,-Inf,-4,-.8))
for(i in seq_len(nrow(cases))) {
 emit(paste0('rectangle.',i,'.bounds'),cases[i,])
 emit(paste0('rectangle.',i,'.probability'),do.call(rectangle,as.list(cases[i,])))
}
d<-data.frame(x=c(0,0,2,2),y=c(0,2,0,2))
fit<-sem('x ~~ 1*x; y ~~ 1*y; x ~~ 0*y; x ~ 0*1; y ~ 0*1',data=d,fixed.x=FALSE,meanstructure=TRUE)
emit('srmr.fixedMean',unlist(lavaan:::lav_fit_srmr_lavobject(fit,'srmr')))
# lavaan suppresses global tests for a zero-free-parameter model; use its
# Gaussian likelihood-ratio identity for this analytic fixed-moment case.
emit('srmr.fixedMean.chi',nrow(d)*sum(colMeans(d)^2))
emit('largeLocationLL',-6*(log(2*pi)+1))
emit('singularSampleLL',-4*(log(2*pi)+1+log(2)))
close(con)
cat('Wrote',destination,'\n')
