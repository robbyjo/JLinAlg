# Frozen SVA, ComBat, PCA, and AutoSVA reference generator.
# Usage:
#   Rscript confounder-reference.R SVA_SOURCE AUTO_SVA_SOURCE OUTPUT_DIRECTORY
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("expected SVA_SOURCE AUTO_SVA_SOURCE OUTPUT_DIRECTORY")
sva_source <- normalizePath(args[[1]], mustWork=TRUE)
auto_source <- normalizePath(args[[2]], mustWork=TRUE)
output <- args[[3]]
dir.create(output, recursive=TRUE, showWarnings=FALSE)

# Load the precise helpers used by sva 3.60.0 without requiring an installed
# Bioconductor package or its parallel executor. The native monotone helper is
# exactly a cumulative maximum in src/sva.c.
rowVars <- function(x, na.rm=FALSE) apply(x, 1, var, na.rm=na.rm)
bpparam <- function(...) NULL
bplapply <- function(X, FUN, ..., BPPARAM=NULL) lapply(X, FUN, ...)
helper_lines <- readLines(file.path(sva_source, "R", "helper.R"))
eval(parse(text=helper_lines[!grepl("^library\\(BiocParallel\\)", helper_lines)]))
mono <- function(x) cummax(x)
source(file.path(sva_source, "R", "f.pvalue.R"))
source(file.path(sva_source, "R", "irwsva.build.R"))
source(file.path(sva_source, "R", "ComBat.R"))

set.seed(1949)
n <- 16L
m <- 48L
case <- rep(c(0,1), length.out=n)
batch <- rep(c("A","B"), each=n/2)
hidden <- sin(seq_len(n)*0.61) + ifelse(batch=="B", 0.7, -0.7)
dat <- matrix(0, m, n)
for (i in seq_len(m)) {
  dat[i,] <- (i %% 4 - 1.5)*hidden + ifelse(i %% 5 == 0, 1.2*case, 0) + rnorm(n,0,.12)
}
rownames(dat) <- paste0("g",seq_len(m))
colnames(dat) <- paste0("s",seq_len(n))
mod <- model.matrix(~case)
mod0 <- model.matrix(~1, data=data.frame(case=case))
write.table(cbind(feature_id=rownames(dat), dat), file.path(output,"input-data.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(sample_id=colnames(dat), case=case, batch=batch), file.path(output,"input-pheno.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

set.seed(7201)
sva_fit <- irwsva.build(dat, mod, mod0, n.sv=2, B=5)
colnames(sva_fit$sv) <- paste0("SV", seq_len(ncol(sva_fit$sv)))
write.table(cbind(sample_id=colnames(dat), sva_fit$sv), file.path(output,"sva-factors.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(feature_id=rownames(dat), pprob_gam=sva_fit$pprob.gam, pprob_b=sva_fit$pprob.b), file.path(output,"sva-weights.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

pca_fit <- prcomp(t(dat), center=TRUE, scale.=FALSE, rank.=2)
write.table(cbind(sample_id=colnames(dat), pca_fit$x[,1:2]), file.path(output,"pca-factors.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

combat_fit <- ComBat(dat, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
write.table(cbind(feature_id=rownames(dat), combat_fit), file.path(output,"combat-adjusted.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# AutoSVA defines its own irwsva.build; fixed K avoids the external mclapply
# dependency and freezes the scientific core separately from the K search.
source(auto_source)
set.seed(7201)
auto_fit <- irwsva.build(dat, mod, mod0, n.sv=2, B=20)
write.table(cbind(sample_id=colnames(dat), auto_fit$sv), file.path(output,"autosva-factors.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(feature_id=rownames(dat), weight=auto_fit$wt, pprob_gam=auto_fit$pprob.gam, pprob_b=auto_fit$pprob.b), file.path(output,"autosva-weights.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

writeLines(c(
  "generator=src/test/R/confounder-reference.R",
  "sva.version=3.60.0",
  "sva.commit=87a4798db8134fd0952a516489dd555545a92d36",
  "peer.commit=40bc4b2cd92459ce42f44dfe279717436395f3f6",
  paste0("r.version=", R.version.string),
  "seed.data=1949",
  "seed.fit=7201",
  "orientation=features_by_samples"
), file.path(output,"reference.properties"))
