# Generate the coloc.susie regression fixture from coloc's bundled example.
# Tested with coloc 5.2.3 and susieR 0.14.2.
local_library <- normalizePath(".r-reference-lib", mustWork = FALSE)
.libPaths(c(local_library, .libPaths()))

suppressPackageStartupMessages(library(coloc))
data(coloc_test_data, package = "coloc")
options(digits = 17)

fits <- lapply(c("D1", "D2", "D3", "D4"), function(name) {
  runsusie(coloc_test_data[[name]], maxit = 1000)
})
names(fits) <- c("D1", "D2", "D3", "D4")

selected <- list(
  D1.L1 = fits[["D1"]]$lbf_variable[fits[["D1"]]$sets$cs_index[1], ],
  D2.L1 = fits[["D2"]]$lbf_variable[fits[["D2"]]$sets$cs_index[1], ],
  D3.L1 = fits[["D3"]]$lbf_variable[fits[["D3"]]$sets$cs_index[1], ],
  D3.L2 = fits[["D3"]]$lbf_variable[fits[["D3"]]$sets$cs_index[2], ],
  D4.L1 = fits[["D4"]]$lbf_variable[fits[["D4"]]$sets$cs_index[1], ]
)

fixture <- data.frame(snp = coloc_test_data[["D1"]]$snp, selected,
                      check.names = FALSE)
destination <- file.path("src", "test", "resources", "r-reference",
                         "coloc-susie-example.tsv")
writeLines(c(
  "# Derived from coloc 5.2.3 coloc_test_data under GPL-3;",
  "# rows are variants and numeric columns are selected susieR 0.14.2 lbf_variable rows."
), destination)
write.table(fixture, destination, sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE, append = TRUE)

for (pair in list(c("D1", "D2"), c("D3", "D4"))) {
  result <- coloc.susie(fits[[pair[1]]], fits[[pair[2]]])
  cat(pair[1], pair[2], "\n")
  print(result$summary, digits = 16)
}
