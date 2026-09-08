# Independent base-R fixtures and paired timing for the v0.3.0 core audit.
options(digits=17)
out <- "src/test/resources/core-audit"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
rows <- list()
record <- function(id, z) {
  rows[[length(rows)+1L]] <<- data.frame(id=id, statistic=unname(z$statistic),
    p=unname(z$p.value), estimate=if(is.null(z$estimate)) NA_real_ else unname(z$estimate)[1])
}
x <- 1:6; y <- c(2,1,4,3,6,5)
record("pearson", cor.test(x,y))
record("t-one",t.test(c(1,2,4,3)))
record("t-welch",t.test(c(1,2,4,3),c(3,5,4,6)))
record("variance",var.test(c(1,2,4,3),c(3,5,4,6)))
record("kendall-zero", suppressWarnings(cor.test(c(-0,0,1,2,3,4),c(1,2,3,2,5,4),method="kendall")))
record("spearman-zero", suppressWarnings(cor.test(c(-0,0,1,2,3,4),c(1,2,3,2,5,4),method="spearman")))
for(n in c(12,30,49,50,180,1300,5000)) {
  x <- seq_len(n)-1; y <- (x*37)%%(n+1)
  if(n <= 50) record(paste0("kendall-exact-",n),cor.test(x,y,method="kendall",exact=TRUE))
  record(paste0("kendall-asym-",n),cor.test(x,y,method="kendall",exact=FALSE))
  record(paste0("spearman-",n),suppressWarnings(cor.test(x,y,method="spearman")))
}
g <- list(c(1,1,1),c(2,3,4),c(3,5,4))
record("anova-constant-group",oneway.test(unlist(g)~rep(seq_along(g),lengths(g)),var.equal=TRUE))
g <- list(c(1,2,1.5,3),c(2,3,4,4.1),c(3,5,4,6))
for(equal in c(TRUE,FALSE)) record(paste0("anova-",equal),oneway.test(unlist(g)~rep(seq_along(g),lengths(g)),var.equal=equal))
write.csv(do.call(rbind,rows),file.path(out,"reference.csv"),row.names=FALSE,na="NaN")
if ("--fixtures-only" %in% commandArgs(TRUE)) quit(save="no")
bench <- list()
for(n in c(5000,20000)) {
 x <- (0:(n-1))%%101; y <- ((0:(n-1))*37)%%97
 f <- function() cor.test(x,y,method="kendall",exact=FALSE)
 for(i in 1:3) f()
 times <- replicate(5,system.time(for(i in 1:3) z<-f())[["elapsed"]]/3)
 z <- f()
 bench[[length(bench)+1]]<-data.frame(n=n,median_seconds=median(times),checksum=unname(z$estimate))
}
dir.create("src/benchmark/resources/core-audit",recursive=TRUE,showWarnings=FALSE)
write.csv(do.call(rbind,bench),"src/benchmark/resources/core-audit/r-timing.csv",row.names=FALSE)
print(do.call(rbind,bench)); print(sessionInfo())
