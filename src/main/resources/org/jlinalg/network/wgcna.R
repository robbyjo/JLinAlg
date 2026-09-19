# Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later
settings <- read.delim("settings.tsv", colClasses="character", check.names=FALSE)
opt <- setNames(settings$value, settings$key)
if ("r-library" %in% names(opt)) .libPaths(c(opt[["r-library"]], .libPaths()))
if (!requireNamespace("WGCNA", quietly=TRUE)) stop("Install WGCNA and its dependencies before running this adapter")
library(WGCNA)
set.seed(as.integer(opt[["seed"]]))
write_tsv <- function(x, path) write.table(x, path, sep="\t", quote=TRUE, row.names=FALSE, na="NA")
read_matrix <- function(path) {
  x <- read.delim(path, row.names=1, check.names=FALSE)
  if (!all(vapply(x, is.numeric, logical(1))) || any(!is.finite(as.matrix(x)))) stop("Non-numeric matrix")
  x
}
x <- read_matrix("matrix.tsv")
qc <- goodSamplesGenes(x, verbose=0)
if (!qc$allOK) stop("WGCNA quality checks failed; filter explicitly upstream")
net <- blockwiseModules(x, power=as.numeric(opt[["power"]]), networkType="signed",
  TOMType="signed", corType="pearson", minModuleSize=as.integer(opt[["min-module-size"]]),
  mergeCutHeight=as.numeric(opt[["merge-cut-height"]]), numericLabels=TRUE,
  randomSeed=as.integer(opt[["seed"]]), maxBlockSize=ncol(x), nThreads=1, verbose=0)
colors <- labels2colors(net$colors)
mes <- net$MEs
membership <- cor(x, mes, use="all.obs")
modules <- data.frame(gene=colnames(x), module=colors, module_label=net$colors, check.names=FALSE)
own <- match(paste0("ME",net$colors), colnames(mes))
modules$membership <- vapply(seq_len(ncol(x)), function(i) if(is.na(own[i])) NA_real_ else membership[i,own[i]], numeric(1))
write_tsv(data.frame(sample_id=rownames(x),mes,check.names=FALSE),"eigengenes.tsv")
adj <- adjacency(x, power=as.numeric(opt[["power"]]), type="signed", corFnc="cor", corOptions=list(use="all.obs"))
modules$within_module_connectivity <- vapply(seq_len(ncol(x)), function(i) {
  if(net$colors[i]==0) return(NA_real_)
  sum(adj[i,which(net$colors==net$colors[i] & seq_len(ncol(x))!=i)])
},numeric(1))
write_tsv(modules, "modules.tsv")
pick <- which(upper.tri(adj) & adj > as.numeric(opt[["edge-threshold"]]), arr.ind=TRUE)
write_tsv(data.frame(source=colnames(x)[pick[,1]], target=colnames(x)[pick[,2]], weight=adj[pick], type=rep("signed_coexpression_adjacency",nrow(pick))),"edges.tsv")
if (file.exists("pheno.tsv")) {
  pheno <- read.delim("pheno.tsv",check.names=FALSE)
  trait <- pheno$trait[match(rownames(x),pheno$sample_id)]
  if(anyNA(trait)||sd(trait)==0) stop("Trait must be complete and variable")
  correlations <- as.numeric(cor(mes,trait,use="all.obs"))
  probabilities <- corPvalueStudent(correlations,nrow(x))
  write_tsv(data.frame(module=colnames(mes),correlation=correlations,p=probabilities,bh=p.adjust(probabilities,"BH")),"module-trait.tsv")
}
if (file.exists("reference.tsv")) {
  ref <- read_matrix("reference.tsv")
  mp <- modulePreservation(list(discovery=list(data=x),replication=list(data=ref)),
    list(discovery=colors), referenceNetworks=1, nPermutations=as.integer(opt[["permutations"]]),
    randomSeed=as.integer(opt[["seed"]]), networkType="signed", quickCor=0, verbose=0)
  saveRDS(mp,"preservation.rds")
  z <- mp$preservation$Z[[1]][[2]]
  write_tsv(data.frame(module=rownames(z),z,check.names=FALSE),"preservation.tsv")
}
writeLines(capture.output(sessionInfo()),"sessionInfo.txt")
