.libPaths(c('build/r-library', .libPaths()))
suppressPackageStartupMessages(library(lavaan))
suppressPackageStartupMessages(library(quantreg))
options(digits=17)
out <- 'src/test/resources/feasible-extensions'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
emit <- function(con,key,value) writeLines(paste0(key,'=',paste(format(value,digits=17,trim=TRUE,scientific=TRUE),collapse=',')),con)
writeLines(c(R.version.string,paste('lavaan',packageVersion('lavaan')),
             paste('quantreg',packageVersion('quantreg'))),file.path(out,'versions.txt'))

# Same frozen ordinal sample, response-independent missingness. Every row has
# at least two observed entries. This targets available-pair, not MAR FIML.
d <- read.csv('src/test/resources/r-reference/sem-joint/ordinal.csv')
for(j in 1:4) d[seq_len(nrow(d)) %% (4+j)==0,j] <- NA
d <- d[rowSums(!is.na(d))>=2,]
fit <- cfa('f =~ a*z1 + b*z2 + c*z3 + d*z4',data=d,ordered=names(d),
           std.lv=TRUE,parameterization='theta',estimator='PML',missing='pairwise')
stopifnot(lavInspect(fit,'converged'))
write.table(d,file.path(out,'ordinal-missing.csv'),row.names=FALSE,col.names=FALSE,sep=',',na='-1',quote=FALSE)
con <- file(file.path(out,'ordinal-missing.properties'),'w')
emit(con,'estimates',coef(fit));emit(con,'covariance',as.vector(t(vcov(fit))))
close(con)

# rq coefficients are independent of Java's LP; the density sandwich is
# independently evaluated with R QR/solve and both varying and pooled density.
d <- as.matrix(read.csv('src/test/resources/exact-quantile/regular.csv',header=FALSE))
y <- d[,1]; x <- d[,-1]; tau <- .35
fit <- rq(y~x-1,tau=tau,method='fn')
f <- .15 + seq_along(y)/(4*length(y))
A <- crossprod(x,f*x); inv <- solve(A)
cov <- tau*(1-tau)*inv %*% crossprod(x) %*% inv
h <- .45
f0 <- mean(dnorm(residuals(fit)/h))/h
con <- file(file.path(out,'quantile.properties'),'w')
emit(con,'tau',tau);emit(con,'estimates',coef(fit));emit(con,'densities',f)
emit(con,'covariance',as.vector(t(cov)));emit(con,'bandwidth',h);emit(con,'kernelDensity',f0)
emit(con,'kernelCovariance',as.vector(t(tau*(1-tau)/f0^2*solve(crossprod(x)))))
close(con)

# Refit existing frozen ARIMA samples; stats::arima uses finite diffuse
# initialization. Compare information away from unresolved cancellation modes.
for(name in c('integrated-ar','integrated-ma','integrated-missing','seasonal',
              'seasonal-only','twice-integrated','drift','seasonal-drift','monthly-missing')) {
  y <- read.csv(paste0('src/test/resources/timeseries/',name,'-series.csv'))$value
  ref <- read.csv(paste0('src/test/resources/timeseries/',name,'-reference.csv'))
  v <- setNames(ref$value,ref$key)
  s <- v[['period']]; drift <- v[['drift']]==1
  fit <- arima(y,order=v[c('p','d','q')],seasonal=list(order=v[c('P','D','Q')],period=s),
               include.mean=FALSE,method='ML',SSinit='Rossignol2011',
               xreg=if(drift)seq_along(y) else NULL,
               optim.control=list(reltol=1e-12,maxit=2000))
  stopifnot(fit$code==0)
  transform <- rep(1,length(fit$coef))
  if(drift && v[['D']]>0) transform[length(transform)] <- s
  cov <- fit$var.coef * outer(transform,transform)
  con <- file(file.path(out,paste0(name,'.properties')),'w')
  emit(con,'estimates',unname(fit$coef)*transform);emit(con,'covariance',as.vector(t(cov)))
  close(con)
}
