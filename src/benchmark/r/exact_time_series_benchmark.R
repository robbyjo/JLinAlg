.libPaths(c('build/r-library', .libPaths()))
# Full stats::arima ML fits matching ExactTimeSeriesBenchmark. No CSS comparisons.
args <- commandArgs(trailingOnly=TRUE)
warmups <- if(length(args)>0) as.integer(args[1]) else 3L
repetitions <- if(length(args)>1) as.integer(args[2]) else 7L
root <- 'src/test/resources/timeseries'
cases <- c('stationary','stationary-missing','integrated-ar','integrated-ma','integrated-missing',
           'seasonal','seasonal-missing','seasonal-only','twice-integrated','drift',
           'seasonal-drift','seasonal-ar-ma','monthly-missing')
cat('case,r_median_ms,loglik_error,converged,checksum\n')
checksum <- 0
for(name in cases) {
  y <- read.csv(file.path(root,paste0(name,'-series.csv')))$value
  tab <- read.csv(file.path(root,paste0(name,'-reference.csv')))
  ref <- setNames(tab$value,tab$key)
  run <- function() arima(y, order=unname(ref[c('p','d','q')]),
    seasonal=list(order=unname(ref[c('P','D','Q')]),period=unname(ref['period'])),
    include.mean=FALSE,method='ML',SSinit='Rossignol2011',
    xreg=if(ref['drift']==1) seq_along(y) else NULL,
    optim.control=list(reltol=1e-12,maxit=2000))
  consume <- function(f) { stopifnot(f$code==0,abs(f$loglik-ref['loglik'])<1e-7); f$loglik+f$sigma2+sum(f$coef) }
  for(i in seq_len(warmups)) checksum <- checksum+consume(run())
  times <- numeric(repetitions)
  for(i in seq_len(repetitions)) {
    t0 <- as.numeric(Sys.time()); f <- run()
    if(ref['d']+ref['D']>0) predict(f,n.ahead=8,newxreg=if(ref['drift']==1) length(y)+1:8 else NULL)
    times[i] <- 1000*(as.numeric(Sys.time())-t0)
    checksum <- checksum+consume(f)
  }
  cat(sprintf('%s,%.6f,%.9g,%s,%.9g\n',name,median(times),f$loglik-ref['loglik'],f$code==0,checksum))
}
