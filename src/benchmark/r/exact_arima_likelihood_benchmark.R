.libPaths(c('build/r-library', .libPaths()))
# Same fixed parameters and full state initialization as ExactArimaLikelihoodBenchmark.
# R's approximate diffuse kappa=1e6; Java evaluates the symbolic diffuse limit.
args <- commandArgs(trailingOnly=TRUE)
batch <- if(length(args)>0) as.integer(args[1]) else 1000L
repetitions <- if(length(args)>1) as.integer(args[2]) else 11L
root <- 'src/test/resources/timeseries'
cases <- c('stationary','stationary-missing','integrated-ar','integrated-ma','integrated-missing',
           'seasonal','seasonal-missing','seasonal-only','twice-integrated','drift',
           'seasonal-drift','seasonal-ar-ma','monthly-missing')
conv <- function(a,b) as.numeric(convolve(a,rev(b),type='open'))
cat('case,r_likelihood_ms,loglik_error,checksum\n')
checksum <- 0
for(name in cases) {
  y <- read.csv(file.path(root,paste0(name,'-series.csv')))$value
  tab <- read.csv(file.path(root,paste0(name,'-reference.csv')))
  ref <- setNames(tab$value,tab$key)
  coefs <- ref[grepl('^coef',names(ref))]
  offset <- 0L
  take <- function(n) { x <- coefs[offset+seq_len(n)]; offset <<- offset+n; unname(x) }
  ar <- take(ref['p']); ma <- take(ref['q']); sar <- take(ref['P']); sma <- take(ref['Q'])
  if(ref['drift']==1) y <- y-take(1)*seq_along(y)
  s <- ref['period']
  sp <- function(x,sign) { z <- numeric(length(x)*s+1); z[1] <- 1; if(length(x)) z[1+seq_along(x)*s] <- sign*x; z }
  phi <- -conv(c(1,-ar),sp(sar,-1))[-1]; theta <- conv(c(1,ma),sp(sma,1))[-1]
  delta <- 1
  for(i in seq_len(ref['d'])) delta <- conv(delta,c(1,-1))
  for(i in seq_len(ref['D'])) { z <- numeric(s+1); z[c(1,s+1)] <- c(1,-1); delta <- conv(delta,z) }
  delta <- -delta[-1]
  run <- function() {
    model <- stats:::makeARIMA(phi,theta,delta,SSinit='Rossignol2011')
    value <- .Call(stats:::C_ARIMA_Like,y,model,0L,FALSE)
    -.5*(value[3]*(log(2*pi)+1+log(value[1]/value[3]))+value[2])
  }
  times <- numeric(repetitions)
  for(repetition in seq.int(-2L,repetitions)) {
    start <- as.numeric(Sys.time())
    for(i in seq_len(batch)) { ll <- run(); checksum <- checksum+ll }
    if(repetition>0) times[repetition] <- 1000*(as.numeric(Sys.time())-start)/batch
  }
  stopifnot(abs(ll-ref['loglik'])<1e-7)
  cat(sprintf('%s,%.6f,%.9g,%.9g\n',name,median(times),ll-ref['loglik'],checksum))
}
