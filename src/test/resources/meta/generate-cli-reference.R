# Independent metafor fixtures for the streaming CLI and primitive-array path.
# Run from the repository root with Rscript src/test/resources/meta/generate-cli-reference.R
if (dir.exists("build/r-library")) .libPaths(c("build/r-library", .libPaths()))
library(metafor)
y <- c(.2,.5,.1,.7); se <- c(.1,.2,.15,.25); dose <- c(-1,0,1,2)
lines <- c(paste0("# R ", getRversion(), "; metafor ", packageVersion("metafor")))
add <- function(key, values) {
  for (i in seq_along(values)) lines <<- c(lines, paste0(key, ".", i-1, "=", sprintf("%.17g",values[i])))
}
for (kind in c("pool", "regression")) for (method in c("FE","REML","DL","PM")) {
  fit <- if (kind=="regression")
    rma(yi=y, sei=se, mods=~dose, method=method, control=list(threshold=1e-12,tol=1e-12,maxiter=1000)) else
    rma(yi=y, sei=se, method=method, control=list(threshold=1e-12,tol=1e-12,maxiter=1000))
  key <- paste(kind,method,sep=".")
  add(paste0(key,".beta"),coef(fit)); add(paste0(key,".se"),fit$se)
  add(paste0(key,".p"),fit$pval); add(paste0(key,".tau2"),fit$tau2)
  add(paste0(key,".q"),fit$QE)
}
writeLines(lines,"src/test/resources/meta/cli-reference.properties")
