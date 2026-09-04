# Copyright (C) 2026 JLinAlg contributors
# SPDX-License-Identifier: GPL-2.0-or-later

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args) >= 1L) args[[1L]] else
  file.path("src", "benchmark", "resources", "susie",
            "N3finemapping.bin.gz")
local_library <- if (length(args) >= 2L) args[[2L]] else ".r-reference-lib"
if (dir.exists(local_library))
  .libPaths(c(normalizePath(local_library), .libPaths()))

if (!requireNamespace("susieR", quietly = TRUE))
  stop("install susieR before generating the benchmark fixture")

data("N3finemapping", package = "susieR", envir = environment())
X <- N3finemapping$X
y <- N3finemapping$Y[, 1L]
true_coefficient <- N3finemapping$true_coef[, 1L]

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
connection <- gzfile(output, open = "wb", compression = 9L)
writeBin(charToRaw("JLSUSIE1"), connection, endian = "big")
writeBin(as.integer(c(nrow(X), ncol(X))), connection,
         size = 4L, endian = "big")
writeBin(as.double(t(X)), connection, size = 8L, endian = "big")
writeBin(as.double(y), connection, size = 8L, endian = "big")
writeBin(as.double(true_coefficient), connection, size = 8L, endian = "big")
close(connection)
