# Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later
settings <- read.delim("settings.tsv",colClasses="character",check.names=FALSE)
opt <- setNames(settings$value,settings$key)
if ("r-library" %in% names(opt)) .libPaths(c(opt[["r-library"]],.libPaths()))
if (!requireNamespace("GENIE3",quietly=TRUE)) stop("Install Bioconductor GENIE3 before running this adapter")
x <- t(as.matrix(read.delim("matrix.tsv",row.names=1,check.names=FALSE)))
regulators <- read.delim("regulators.tsv",colClasses="character")$gene
set.seed(as.integer(opt[["seed"]]))
weights <- GENIE3::GENIE3(x,regulators=regulators,nTrees=as.integer(opt[["trees"]]),nCores=1,verbose=FALSE)
links <- GENIE3::getLinkList(weights,threshold=as.numeric(opt[["edge-threshold"]]))
colnames(links) <- c("source","target","weight")
links$type <- rep("predicted_regulation",nrow(links))
# GENIE3 importance does not identify activation/inhibition.
links$sign <- rep("unknown",nrow(links))
write.table(links,"edges.tsv",sep="\t",quote=TRUE,row.names=FALSE,na="NA")
saveRDS(weights,"weights.rds")
writeLines(capture.output(sessionInfo()),"sessionInfo.txt")
