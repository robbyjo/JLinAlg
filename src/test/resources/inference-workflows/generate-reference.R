# Frozen limma/voom, edgeR and DESeq2 comparison fixture.
# Run from the repository root after installing the packages into build/r-lib.
options(digits=17)
.libPaths(c(normalizePath('build/r-lib'), .libPaths()))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DESeq2))

destination <- 'src/test/resources/inference-workflows'
dir.create(destination, recursive=TRUE, showWarnings=FALSE)
write_tab <- function(value, name) write.table(value,
  file.path(destination, name), sep='\t', quote=FALSE, row.names=FALSE)
write_matrix <- function(value, name) {
  write_tab(data.frame(feature_id=rownames(value), value,
    check.names=FALSE), name)
}

set.seed(17092026)
samples <- paste0('s', seq_len(8))
group <- factor(rep(c('control','case'), each=4),
  levels=c('control','case'))
design <- model.matrix(~group)
features <- paste0('feature_', sprintf('%02d', seq_len(48)))

baseline <- seq(2, 8, length.out=length(features))
effect <- c(rep(1.25, 8), rep(-0.75, 6), rep(0, 34))
residual_sd <- seq(0.18, 0.65, length.out=length(features))
continuous <- matrix(0, length(features), length(samples),
  dimnames=list(features, samples))
for (feature in seq_along(features)) {
  continuous[feature,] <- baseline[feature] +
    effect[feature] * as.numeric(group == 'case') +
    rnorm(length(samples), sd=residual_sd[feature])
}
write_matrix(continuous, 'continuous.tsv')
limma_fit <- eBayes(lmFit(continuous, design), trend=FALSE,
  robust=FALSE)
write_tab(data.frame(feature_id=features,
  effect=limma_fit$coefficients[,2],
  t=limma_fit$t[,2], p=limma_fit$p.value[,2],
  moderated_variance=limma_fit$s2.post,
  df_total=limma_fit$df.total), 'limma.tsv')

library_factor <- c(0.72, 0.89, 1.02, 1.17, 0.81, 0.96, 1.19, 1.31)
base_mean <- exp(seq(log(12), log(700), length.out=length(features)))
log2_fold_change <- c(rep(1.1, 8), rep(-0.8, 6), rep(0, 34))
dispersion <- seq(0.04, 0.28, length.out=length(features))
counts <- matrix(0L, length(features), length(samples),
  dimnames=list(features, samples))
for (feature in seq_along(features)) {
  mean <- base_mean[feature] * library_factor *
    2^(log2_fold_change[feature] * as.numeric(group == 'case'))
  counts[feature,] <- rnbinom(length(samples), mu=mean,
    size=1/dispersion[feature])
}
write_matrix(counts, 'counts.tsv')

dge <- DGEList(counts=counts, group=group)
dge <- normLibSizes(dge)
dge <- estimateDisp(dge, design, robust=FALSE)
edge_fit <- glmQLFit(dge, design, robust=FALSE)
edge_test <- glmQLFTest(edge_fit, coef=2)
write_tab(data.frame(feature_id=features,
  log2_fold_change=edge_test$table$logFC,
  p=edge_test$table$PValue,
  dispersion=edge_fit$dispersion), 'edger.tsv')

voom_data <- voom(dge, design, plot=FALSE)
voom_fit <- eBayes(lmFit(voom_data, design), trend=FALSE,
  robust=FALSE)
write_tab(data.frame(feature_id=features,
  log2_fold_change=voom_fit$coefficients[,2],
  t=voom_fit$t[,2], p=voom_fit$p.value[,2],
  moderated_variance=voom_fit$s2.post,
  df_total=voom_fit$df.total), 'voom.tsv')

dds <- DESeqDataSetFromMatrix(countData=counts,
  colData=data.frame(group=group, row.names=samples), design=~group)
dds <- DESeq(dds, quiet=TRUE)
deseq <- results(dds, contrast=c('group','case','control'),
  independentFiltering=FALSE, cooksCutoff=FALSE)
write_tab(data.frame(feature_id=features,
  log2_fold_change=deseq$log2FoldChange,
  p=deseq$pvalue,
  dispersion=dispersions(dds)), 'deseq2.tsv')

versions <- data.frame(component=c('R','Bioconductor','limma','edgeR','DESeq2'),
  version=c(as.character(getRversion()), as.character(BiocManager::version()),
    as.character(packageVersion('limma')),
    as.character(packageVersion('edgeR')),
    as.character(packageVersion('DESeq2'))))
write_tab(versions, 'versions.tsv')
cat('Generated empirical-Bayes differential reference fixture\n')
