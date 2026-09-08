.libPaths(c('build/r-library', .libPaths()))
# Base-R references for the four independent quadrature review counterexamples.
# No fitted Java parameters are used to generate these expectations.
options(digits=17)
root <- 'src/test/resources/r-reference'
sizes <- c(64,65,1e4,1e6,1e8,1e10,1e12,1e14,1e16)
mass <- data.frame(n=sizes,binomial=vapply(sizes,function(n) dbinom(floor(n/2),n,.5,log=TRUE),0.),
  poisson=vapply(sizes,function(n) dpois(n,n,log=TRUE),0.))
write.table(mass,file.path(root,'quadrature-large-counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
properties <- paste0('R.version=',getRversion())
put <- function(key,value) properties <<- c(properties,paste0(key,'=',format(value,digits=17,scientific=FALSE,trim=TRUE)))
n <- 3000000; k <- 1225; v0 <- 4/n
group <- function(s,beta,v) {
  if(v==0) return(dbinom(s,n,plogis(beta),log=TRUE))
  sd <- sqrt(v)
  score <- function(u) sd*(s-n*plogis(beta+sd*u))-u
  mode <- uniroot(score,c(-50,50),tol=1e-12)$root
  curvature <- 1+v*n*plogis(beta+sd*mode)*plogis(-beta-sd*mode)
  scale <- 1/sqrt(curvature)
  logf <- function(z) dbinom(s,n,plogis(beta+sd*(mode+scale*z)),log=TRUE)+dnorm(mode+scale*z,log=TRUE)+log(scale)
  peak <- logf(0)
  peak+log(integrate(function(z) exp(logf(z)-peak),-Inf,Inf,rel.tol=1e-12,abs.tol=0,subdivisions=2000)$value)
}
lik <- function(beta,t) group(n/2-k,beta,v0*expm1(t))+group(n/2+k,beta,v0*expm1(t))
best <- optimize(function(t) lik(0,t),c(0,3),maximum=TRUE,tol=1e-10)
t <- best$maximum; v <- v0*expm1(t); h <- 1e-5; ht <- 1e-3
put('boundary.n',n);put('boundary.k',k);put('boundary.variance',v)
put('boundary.logLikelihood',best$objective);put('boundary.zeroLogLikelihood',lik(0,0))
put('boundary.rightNllDerivative',n/4-k*k)
put('boundary.betaSE',sqrt(1/(-(lik(h,t)-2*lik(0,t)+lik(-h,t))/h^2)))
put('boundary.varianceCovariance',(v+v0)^2/(-(lik(0,t+ht)-2*lik(0,t)+lik(0,t-ht))/ht^2))
# A shallower interior optimum also fooled the old boundary derivative.
n <- 300000; k <- 287; v0 <- 4/n
best <- optimize(function(t) lik(0,t),c(0,3),maximum=TRUE,tol=1e-10)
put('shallow.n',n);put('shallow.k',k)
put('shallow.variance',v0*expm1(best$maximum))
put('shallow.logLikelihood',best$objective);put('shallow.zeroLogLikelihood',lik(0,0))
for(offset in c(100,300,500,1000)) {
  mode <- uniroot(function(u) 1-exp(offset+u)-u,c(-2*offset,0),tol=1e-12)$root
  scale <- 1/sqrt(1+exp(offset+mode))
  logf <- function(z) dpois(1,exp(offset+mode+scale*z),log=TRUE)+dnorm(mode+scale*z,log=TRUE)+log(scale)
  peak <- logf(0)
  value <- peak+log(integrate(function(z) exp(logf(z)-peak),-40,40,rel.tol=1e-12,abs.tol=0,subdivisions=2000)$value)
  put(paste0('offset.',offset),value)
}
put('proportion.trials',354343171);put('proportion.successes',186369615)
put('proportion.logLikelihood',dbinom(186369615,354343171,186369615/354343171,log=TRUE))
writeLines(properties,file.path(root,'quadrature-stress.properties'))
print(mass);cat(paste(properties,collapse='\n'),'\n')
