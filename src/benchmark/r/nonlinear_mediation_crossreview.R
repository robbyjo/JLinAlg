.libPaths(c('build/r-library', .libPaths()))
# Independent edge references; no fitting tolerance is enlarged for quantized data.
options(digits=17)
levels <- c(.95, .999999999999, 1 - .Machine$double.eps/2)
tails <- do.call(rbind, lapply(levels, function(level) {
  p <- (1-level)/2
  c(level, qnorm(p, lower.tail=FALSE), qt(p,114:117,lower.tail=FALSE))
}))
write.table(apply(tails,2,function(x)sprintf('%.17g',x)),
            'src/test/resources/remaining-model-audit/mediation-tail.csv',
            sep=',',row.names=FALSE,col.names=FALSE,quote=FALSE)
print(tails)
i <- 0:99; x <- i%%9-4
for (shift in c(0,1e12,1e15)) {
  y <- shift + 2*x + .25*ifelse(i%%2==0,1,-1)
  # Solve the centered linear problem, then evaluate the original mean.
  b <- sum(x*(y-shift))/sum(x*x)
  cat('shift',shift,'beta',b,'evaluated_SSE',sum((y-(shift+b*x))^2),'\n')
}
