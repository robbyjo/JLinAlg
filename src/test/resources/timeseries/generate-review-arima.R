.libPaths(c('build/r-library', .libPaths()))
options(digits=17)
root <- 'src/test/resources/timeseries'
ma <- c(-1.5699824440879728,1.7418872190064827,1.0916093731684835,
        .2963556882787251,.8255498200370827,-.9177455894496253)
mafit <- arima(ma,order=c(0,0,1),include.mean=FALSE,method='ML',init=-.15,
               optim.control=list(reltol=1e-13,maxit=2000))
arfit <- arima(c(1,-2,3,-1),order=c(2,0,0),include.mean=FALSE,method='ML',
               init=c(-1,-.5),optim.control=list(reltol=1e-13,maxit=2000))
small <- arima(1e-8*c(1,-2,3,-1,2),order=c(0,0,0),include.mean=FALSE,method='ML')
y <- rep(NA_real_,101); y[c(1,2,51,101)] <- c(0,1,5000,10000)
drift <- arima(y,order=c(0,1,0),xreg=seq_along(y),fixed=100,method='ML')
z <- rep(NA_real_,104); z[1:4] <- 0; z[5:8] <- 1; z[101:104] <- 10000
seasonal <- arima(z,seasonal=list(order=c(0,1,0),period=4),xreg=seq_along(z)/4,
                  fixed=400,method='ML')
ref <- c(ma1=unname(mafit$coef),ma_loglik=mafit$loglik,ar1=unname(arfit$coef[1]),
         ar2=unname(arfit$coef[2]),ar_loglik=arfit$loglik,small_variance=small$sigma2,
         small_loglik=small$loglik,drift_loglik=drift$loglik,seasonal_drift_loglik=seasonal$loglik)
write.table(data.frame(key=names(ref),value=unname(ref)),file.path(root,'review-reference.csv'),
            row.names=FALSE,sep=',',quote=FALSE)
print(ref)
