# Copyright (C) 2026 JLinAlg contributors
# SPDX-License-Identifier: GPL-2.0-or-later
args <- commandArgs(trailingOnly=TRUE)
mode <- args[1]
if(args[2] != "-") .libPaths(c(args[2], .libPaths()))
tab <- function(path) read.delim(path,check.names=FALSE,stringsAsFactors=FALSE,
  colClasses="character",na.strings=character(),quote='"',comment.char="",fileEncoding="UTF-8")
write <- function(x,path) write.table(x,path,sep="\t",row.names=FALSE,quote=TRUE,
  qmethod="double",na="NA",fileEncoding="UTF-8")
if(mode == "state") {
  if(!requireNamespace("limma",quietly=TRUE)) stop("Install limma in the chosen R library before using state")
  raw <- tab("counts.tsv"); counts <- as.matrix(raw[,-1,drop=FALSE]); storage.mode(counts) <- "double"; rownames(counts) <- raw[[1]]
  design <- tab("design.tsv"); stopifnot(identical(colnames(counts),design[[1]]))
  x <- as.matrix(design[,-1,drop=FALSE]); storage.mode(x) <- "double"
  stopifnot(qr(x)$rank==ncol(x),nrow(x)>ncol(x),all(colSums(counts)>0))
  eligible <- rowSums(counts>=as.integer(args[3]))>=as.integer(args[4])
  answer <- data.frame(feature_id=rownames(counts),effect=NA_real_,se=NA_real_,lower=NA_real_,upper=NA_real_,df=NA_real_,statistic=NA_real_,p=NA_real_,status="low_counts")
  write(data.frame(sample_id=colnames(counts),library_count=colSums(counts),normalization_factor=1),"libraries.tsv")
  if(sum(eligible)>=2) {
    # Preserve full assayed-library sizes after feature filtering. No TMM or global-shift correction.
    v <- limma::voom(counts[eligible,,drop=FALSE],design=x,lib.size=colSums(counts),normalize.method="none",plot=FALSE)
    fit <- limma::eBayes(limma::lmFit(v,x))
    se <- fit$stdev.unscaled[,2]*sqrt(fit$s2.post); df <- fit$df.total; critical <- qt(.975,df)
    answer$effect[eligible] <- fit$coefficients[,2]; answer$se[eligible] <- se
    answer$lower[eligible] <- fit$coefficients[,2]-critical*se
    answer$upper[eligible] <- fit$coefficients[,2]+critical*se
    answer$df[eligible] <- df; answer$statistic[eligible] <- fit$t[,2]; answer$p[eligible] <- fit$p.value[,2]
    answer$status[eligible] <- "tested"
    if(any(!is.finite(as.matrix(answer[eligible,2:8])))) stop("Nonfinite limma inference")
    write(data.frame(feature_id=rownames(v),v$E,check.names=FALSE),"log2-cpm.tsv")
    write(data.frame(feature_id=rownames(v),v$weights,check.names=FALSE),"precision-weights.tsv")
  } else answer$status[eligible] <- "insufficient_features_for_voom"
  write(answer,"fit.tsv")
} else if(mode == "representation") {
  set.seed(as.integer(args[3])); raw <- tab("normalized.tsv"); x <- as.matrix(raw[,-1,drop=FALSE]); storage.mode(x) <- "double"
  variance <- apply(x,2,var); candidates <- which(variance>0)
  selected <- candidates[order(-variance[candidates],candidates)]
  selected <- head(selected,as.integer(args[4])); if(length(selected)<2) stop("Need two variable features")
  if(nrow(x)*length(selected)*min(nrow(x),length(selected))>1e8)
    stop("PCA exceeds 100 million SVD work units; reduce variable-features or observations")
  components <- min(as.integer(args[5]),nrow(x)-1,length(selected))
  pca <- prcomp(x[,selected,drop=FALSE],center=TRUE,scale.=FALSE,rank.=components)
  if(sum(pca$sdev>max(pca$sdev)*1e-10)<components) stop("Requested PCA components exceed numerical rank")
  clusters <- kmeans(pca$x,centers=as.integer(args[6]),nstart=20,iter.max=100)
  if(clusters$ifault!=0) stop("k-means did not converge")
  write(data.frame(obs_id=raw[[1]],cluster=clusters$cluster,pca$x,check.names=FALSE),"embedding.tsv")
  write(data.frame(feature_id=colnames(x),variance=variance,selected=seq_along(variance)%in%selected),"variable-features.tsv")
  write(data.frame(feature_id=rownames(pca$rotation),pca$rotation,check.names=FALSE),"loadings.tsv")
  write(data.frame(component=seq_along(pca$sdev),variance=pca$sdev^2),"pca-variance.tsv")
} else stop("Unknown adapter mode")
capture.output(sessionInfo(),file="session.txt")
