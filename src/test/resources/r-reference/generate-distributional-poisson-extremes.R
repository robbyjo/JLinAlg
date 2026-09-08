.libPaths(c('build/r-library', .libPaths()))
options(digits=17)
mu <- c(1e-10,.1,1,20,1e6,1e9,1e12,1e15)
d <- do.call(rbind,lapply(mu,function(m) {
  y <- unique(c(1,round(m),round(m+sqrt(m))))
  data.frame(y=y,mu=m,zip=ifelse(y==0,log(.2+.8*exp(-m)),log(.8)+dpois(y,m,log=TRUE)),
    hurdle=ifelse(y==0,log(.2),log(.8)+dpois(y,m,log=TRUE)-log(-expm1(-m))))
}))
write.table(d,'src/test/resources/r-reference/distributional-poisson-extremes.tsv',sep='\t',quote=FALSE,row.names=FALSE)
