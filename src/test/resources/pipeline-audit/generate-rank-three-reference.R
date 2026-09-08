.libPaths(c('build/r-library',.libPaths()))
# Positive gamma-mixture identity; independently sum fixed 1000 terms and check
# mass/mean, unlike Java's dynamically bounded truncation. Base R only.
lambda <- c(1,.5,.25); beta <- min(lambda); terms <- 1000L
coeff <- numeric(terms+1L); coeff[1] <- prod(sqrt(beta/lambda))
g <- vapply(1:terms,function(j) sum((1-beta/lambda)^j)/2,numeric(1))
for(k in 1:terms) coeff[k+1] <- sum(g[1:k]*rev(coeff[1:k]))/k
stopifnot(abs(sum(coeff)-1)<1e-14,
          abs(sum(coeff*(length(lambda)+2*(0:terms))*beta)-sum(lambda))<1e-13)
q <- c(10,30,50,100,200)
tail <- vapply(q,function(x) sum(coeff*pchisq(x/beta,df=length(lambda)+2*(0:terms),lower.tail=FALSE)),0)
fixture <- data.frame(q=q,tail=tail,upper=pchisq(q,3,lower.tail=FALSE))
writeLines(c('q\ttail\tupper',apply(fixture,1,function(row) paste(sprintf('%.17g',row),collapse='\t'))),
           'src/test/resources/pipeline-audit/rank-three-reference.tsv')
