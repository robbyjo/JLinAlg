.libPaths(c('build/r-library', .libPaths()))
# Final independent-review regressions. R polyroot enumerates the stationary
# phases; Java instead brackets real roots using derivative-root partitions.
options(digits=17)
cases <- list(slow=matrix(c(4,0,2,3,21,4,3,3,1),3,byrow=TRUE),
              twin=matrix(c(2,2,1,1,13,3,1,3,0),3,byrow=TRUE),
              independent=matrix(c(4,2,0,3,14,2,4,2,0),3,byrow=TRUE))
result <- c(susie.mean=1e-200*1e-100, susie.variance=1e-200)
phase_fit <- function(count) {
  n <- sum(count); h <- count[2,2]
  known <- c(2*count[1,1]+count[1,2]+count[2,1],count[1,2]+2*count[1,3]+count[2,3],
             count[2,1]+2*count[3,1]+count[3,2],count[2,3]+count[3,2]+2*count[3,3])
  a<-known[1];b<-known[2];c<-known[3];d<-known[4]
  roots <- polyroot(c(-a*d,a*d+(b+h)*(c+h)-h*(a+d),h*(a+d-b-c)-3*h*h,2*h*h))
  phases <- c(0,1,Re(roots[abs(Im(roots))<1e-8 & Re(roots)>=0 & Re(roots)<=1]))
  evaluate <- function(t) {
    f <- (known+h*c(t,1-t,1-t,t))/(2*n)
    p <- matrix(c(f[1]^2,2*f[1]*f[2],f[2]^2,2*f[1]*f[3],2*(f[1]*f[4]+f[2]*f[3]),
                  2*f[2]*f[4],f[3]^2,2*f[3]*f[4],f[4]^2),3,byrow=TRUE)
    pa<-f[3]+f[4];pb<-f[2]+f[4]
    c(phase=t,LL=sum(count[count>0]*log(p[count>0])),rSquared=(f[4]-pa*pb)^2/(pa*(1-pa)*pb*(1-pb)))
  }
  values <- t(vapply(phases,evaluate,numeric(3))); best <- values[which.max(values[,2]),]
  best
}
for (name in names(cases)) {
  best <- phase_fit(cases[[name]])
  names(best) <- paste0('ld.',name,'.',names(best)); result<-c(result,best)
}
# Normalize already-weighted predictors: no overflowing outcome/exposure ratio.
x <- c(1e-200,1e-200); y <- c(1e200,1); se <- c(1e200,1)
wx<-x/se; wy<-y/se; scale<-max(abs(wx)); z<-wx/scale
result<-c(result,ivw.beta=sum(z*wy)/sum(z*z)/scale,ivw.se=1/sqrt(sum(z*z))/scale)
writeLines(c('# Base R conjugacy, weighted regression, and complete phase-root enumeration',
             paste0(names(result),'=',format(result,digits=17,scientific=TRUE,trim=TRUE))),
           'src/test/resources/r-reference/genetic-final-fixes.properties')
print(result)
# Warm paired solver timing, using the same 41-founder table. No BED I/O or
# package comparison is implied. Include complete root and likelihood selection.
for (i in seq_len(100)) phase_fit(cases$slow)
times <- numeric(5); checksum <- 0; repeats <- 1000L
for (batch in seq_len(5)) {
  elapsed <- system.time(for(i in seq_len(repeats)) {
    value <- phase_fit(cases$slow)[['rSquared']]
    stopifnot(is.finite(value),abs(value/result[['ld.slow.rSquared']]-1)<2e-10)
    checksum <- checksum+value*value
  })[['elapsed']]
  times[batch] <- elapsed*1000/repeats
}
cat(sprintf('LD41_baseR_polyroot median_ms=%.9g checksum=%.17g calls=%d\n',median(times),checksum,5L*repeats))
