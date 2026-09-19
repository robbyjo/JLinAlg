# Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later
settings <- read.delim("settings.tsv",colClasses="character",check.names=FALSE)
opt <- setNames(settings$value,settings$key)
if ("r-library" %in% names(opt)) .libPaths(c(opt[["r-library"]],.libPaths()))
if (!requireNamespace("CARNIVAL",quietly=TRUE)) stop("Install Bioconductor CARNIVAL and the selected solver before running this adapter")
prior <- read.delim("prior.tsv",check.names=FALSE)
activities <- read.delim("activities.tsv",check.names=FALSE)
measurements <- setNames(activities$activity,activities$gene)
set.seed(as.integer(opt[["seed"]]))
solver <- opt[["solver"]]
options <- switch(solver,lpSolve=CARNIVAL::defaultLpSolveCarnivalOptions(),
  cbc=CARNIVAL::defaultCbcSolveCarnivalOptions(),cplex=CARNIVAL::defaultCplexCarnivalOptions())
dir.create("solver-output")
options$outputFolder <- normalizePath("solver-output",winslash="/")
if ("solver-path" %in% names(opt)) options$solverPath <- opt[["solver-path"]]
result <- CARNIVAL::runInverseCarnival(measurements=measurements,priorKnowledgeNetwork=prior,carnivalOptions=options)
saveRDS(result,"solution.rds")
writeLines(capture.output(str(result)),"solution-structure.txt")
edges <- result$weightedSIF
if (is.null(edges)) stop("CARNIVAL returned no weightedSIF solution; inspect solver feasibility/version")
write.table(edges,"weightedSIF.tsv",sep="\t",quote=TRUE,row.names=FALSE,na="NA")
write.table(result$nodesAttributes,"nodes.tsv",sep="\t",quote=TRUE,row.names=FALSE,na="NA")
edges <- edges[edges$Weight>0,,drop=FALSE]
if(nrow(edges)==0) stop("CARNIVAL returned no selected edges")
artificial <- edges$Node1=="Perturbation" | edges$Node2=="Perturbation"
write.table(edges[artificial,,drop=FALSE],"proposed-perturbations.tsv",sep="\t",quote=TRUE,row.names=FALSE,na="NA")
edges <- edges[!artificial,,drop=FALSE]
canonical <- data.frame(source=edges$Node1,target=edges$Node2,sign=edges$Sign,
  weight=edges$Weight/100,type=rep("contextual_signaling",nrow(edges)))
write.table(canonical,"edges.tsv",sep="\t",quote=TRUE,row.names=FALSE,na="NA")
writeLines(capture.output(sessionInfo()),"sessionInfo.txt")
