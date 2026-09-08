# Independent integer-count NB finite product, avoiding dnbinom's asymptotic
# small-count branch. Retain R's result alongside it to expose disagreement.
options(digits=17)
grid <- expand.grid(y=c(0,1,2,10,100), mu=c(.01,2,1e5), size=c(1e8,1e12,1e14))
grid$r_log_density <- with(grid, dnbinom(y, mu=mu, size=size, log=TRUE))
finite_product <- function(y, mu, size) {
  correction <- if (y == 0) 0 else sum(log1p((seq_len(y)-1)/size))
  correction + y*log(mu) - lgamma(y+1) - (size+y)*log1p(mu/size)
}
grid$finite_product_log_density <- mapply(finite_product, grid$y, grid$mu, grid$size)
stopifnot(abs(finite_product(2,1e5,1e14) - (-99977.66724625262)) < 2e-11)
write.table(grid, 'src/test/resources/r-reference/jdistlib-upgrade-nb.tsv',
            sep='\t', quote=FALSE, row.names=FALSE)
print(c(R_version=as.character(getRversion()), rows=nrow(grid),
        maximum_R_discrepancy=max(abs(grid$r_log_density-grid$finite_product_log_density))))
