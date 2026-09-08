.libPaths(c('build/r-library', .libPaths()))
suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(clubSandwich))
args <- commandArgs(trailingOnly=TRUE)
repeats <- if (length(args)) as.integer(args[1]) else 100L
stopifnot(repeats>0)
root <- 'src/test/resources/meta'
dat <- read.csv(file.path(root,'hierarchical.csv'))
V <- as.matrix(read.csv(file.path(root,'sampling.csv'),header=FALSE))
bias <- read.csv(file.path(root,'bias.csv'))
reference <- read.csv(file.path(root,'metafor-reference.csv'))
ref <- setNames(reference$value,reference$key)
expected <- function(a,na,b,nb,c,nc) unname(ref[c(paste0(a,'.',seq_len(na)-1),paste0(b,'.',seq_len(nb)-1),paste0(c,'.',seq_len(nc)-1))])
control <- list(rel.tol=1e-10,iter.max=2000,eval.max=4000)
cases <- list(
  known_gls=function() {
    f <- rma.mv(yi,V,mods=~x,data=dat,method='ML')
    c(coef(f),as.vector(t(vcov(f))),logLik(f))
  },
  random_moderators_reml=function() {
    f <- rma.mv(yi,V,mods=~x,random=~1+x|cluster,struct='GEN',data=dat,method='REML',control=control)
    c(coef(f),as.vector(t(f$G)),logLik(f))
  },
  nested_reml=function() {
    f <- rma.mv(yi,V,mods=~x,random=list(~1|cluster,~1|effect),data=dat,method='REML',control=control)
    c(coef(f),f$sigma2,logLik(f))
  },
  cr2_correlated=function() {
    f <- rma.mv(yi,V,mods=~x,data=dat,method='ML')
    cr <- vcovCR(f,cluster=dat$cluster,type='CR2'); ct <- coef_test(f,vcov=cr)
    c(as.vector(t(cr)),ct$df_Satt,ct$p_Satt)
  },
  trimfill_fe_l0=function() {
    f <- trimfill(rma(yi,sei=sei,data=bias,method='FE'),side='left',estimator='L0')
    c(coef(f),f$se,f$k0)
  },
  pet=function() {
    f <- lm(yi~sei,weights=1/sei^2,data=bias)
    c(coef(f),as.vector(t(vcov(f))),coef(summary(f))[,4])
  })
targets <- list(
  known_gls=expected('known.beta',2,'known.cov',4,'known.ll',1),
  random_moderators_reml=expected('hier.REML.beta',2,'hier.REML.G',4,'hier.REML.ll',1),
  nested_reml=expected('nested.REML.beta',2,'nested.REML.variance',2,'nested.REML.ll',1),
  cr2_correlated=expected('CR2.cov',4,'CR2.df',2,'CR2.p',2),
  trimfill_fe_l0=expected('trim.FE.L0.left.beta',1,'trim.FE.L0.left.se',1,'trim.FE.L0.left.k0',1),
  pet=expected('pet.beta',2,'pet.cov',4,'pet.p',2))
cat('method,repeats,milliseconds_per_fit,max_absolute_reference_error,checksum\n')
for (name in names(cases)) {
  run <- cases[[name]]
  error <- max(abs(run()-targets[[name]])); stopifnot(is.finite(error),error<2e-5)
  for (i in 1:10) run()
  checksum <- 0
  elapsed <- system.time(for (i in seq_len(repeats)) checksum <- checksum + sum(run()))[['elapsed']]
  cat(sprintf('%s,%d,%.9f,%.12g,%.12g\n',name,repeats,1000*elapsed/repeats,error,checksum))
}
