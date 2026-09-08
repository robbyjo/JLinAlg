.libPaths(c('build/r-library', .libPaths()))
# Run from the repository root. Frozen stats::arima ML and predict.Arima references.
options(digits=17)
set.seed(620260908)
out <- 'src/test/resources/timeseries'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
conv <- function(a,b) as.numeric(convolve(a,rev(b),type='open'))
make <- function(name, order, seasonal=c(0,0,0,1), ar=numeric(), ma=numeric(),
                 sar=numeric(), sma=numeric(), missing=integer(), drift=0) {
  s <- seasonal[4]
  sp <- function(x,sign) { z <- numeric(length(x)*s+1); z[1] <- 1; if(length(x)) z[1+seq_along(x)*s] <- sign*x; z }
  phi <- -conv(c(1,-ar),sp(sar,-1))[-1]
  theta <- conv(c(1,ma),sp(sma,1))[-1]
  x <- as.numeric(arima.sim(list(ar=phi,ma=theta),n=240))
  delta <- 1
  for(i in seq_len(order[2])) delta <- conv(delta,c(1,-1))
  for(i in seq_len(seasonal[2])) { z <- numeric(s+1); z[c(1,s+1)] <- c(1,-1); delta <- conv(delta,z) }
  if(length(delta)>1) x <- as.numeric(filter(x,-delta[-1],method='recursive'))
  x <- x + drift*seq_along(x)
  x[missing] <- NA_real_
  fitfun <- function() arima(x, order=order, seasonal=list(order=seasonal[1:3],period=s),
    include.mean=FALSE, method='ML', SSinit='Rossignol2011',
    xreg=if(drift!=0) seq_along(x) else NULL,
    optim.control=list(reltol=1e-12,maxit=2000))
  f <- fitfun()
  stopifnot(f$code==0)
  pred <- predict(f,n.ahead=8,newxreg=if(drift!=0) length(x)+1:8 else NULL)
  write.table(data.frame(value=x),file.path(out,paste0(name,'-series.csv')),row.names=FALSE,sep=',',na='NaN',quote=FALSE)
  ref <- c(p=order[1],d=order[2],q=order[3],P=seasonal[1],D=seasonal[2],Q=seasonal[3],period=s,
           drift=as.integer(drift!=0),coef=unname(f$coef),variance=f$sigma2,loglik=f$loglik,nobs=f$nobs,
           forecast=as.numeric(pred$pred),se=as.numeric(pred$se))
  if(length(delta)>1) for(k in c(1e4,1e5,1e6,1e7)) {
    fixed <- arima(x,order=order,seasonal=list(order=seasonal[1:3],period=s),include.mean=FALSE,
      fixed=unname(f$coef),method='ML',SSinit='Rossignol2011',kappa=k,
      xreg=if(drift!=0) seq_along(x) else NULL)
    ref[paste0('kappa',format(k,scientific=FALSE))] <- fixed$loglik
  }
  write.table(data.frame(key=names(ref),value=unname(ref)),file.path(out,paste0(name,'-reference.csv')),row.names=FALSE,sep=',',quote=FALSE)
  invisible(f)
}
make('stationary',c(1,0,1),ar=.55,ma=.35)
make('stationary-missing',c(1,0,1),ar=.55,ma=.35,missing=c(1:3,20:24,77,180,239:240))
make('integrated-ar',c(1,1,0),ar=.65)
make('integrated-ma',c(0,1,1),ma=.45)
make('integrated-missing',c(1,1,1),ar=.55,ma=.35,missing=c(1:3,20:24,77,180,239:240))
make('seasonal',c(0,1,1),c(0,1,1,4),ma=.4,sma=.3)
make('seasonal-missing',c(0,1,1),c(0,1,1,4),ma=.4,sma=.3,missing=c(1,2,7,20:24,77,180,239:240))
make('seasonal-only',c(0,0,0),c(1,1,0,4),sar=.55)
make('twice-integrated',c(0,2,1),ma=.3,missing=c(1,20:24,239:240))
make('drift',c(1,1,0),ar=.55,drift=.2,missing=c(1,20:24,239:240))
make('seasonal-drift',c(0,0,0),c(1,1,0,4),sar=.55,drift=.2,missing=c(1,20:24,239:240))
make('seasonal-ar-ma',c(1,1,1),c(1,1,1,4),ar=.35,ma=.2,sar=.3,sma=.15,missing=c(1,20:24,239:240))
make('monthly-missing',c(0,1,1),c(0,1,1,12),ma=.4,sma=.3,missing=c(1:5,20:24,77,180,239:240))
writeLines(c(R.version.string, 'stats::arima(method="ML", SSinit="Rossignol2011", kappa=1e6)',
             'Independent seeded Gaussian simulations; forecasts include missing terminal observations.'),file.path(out,'reference-version.txt'))
