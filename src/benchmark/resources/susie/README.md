# SuSiE vignette benchmark data

`N3finemapping.bin.gz` is a deterministic binary export of the first response
replicate in the `N3finemapping` data bundled with susieR 0.14.2. It contains
the 574-by-1,001 design matrix, response, and true coefficients used by the
official susieR fine-mapping vignette.

Source: <https://cran.r-project.org/package=susieR>

The source package is distributed under GPL (>= 3). The fixture is used only
for regression testing and benchmarking. Regenerate it with:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' `
  src/benchmark/r/generate_susie_vignette_data.R
```

Binary layout is gzip-compressed, big-endian `JLSUSIE1`, two 32-bit dimensions,
row-major IEEE-754 doubles for X, doubles for y, then the true coefficients.

Data SHA-256:
`0c00bec50cb8a9354e8e8fecdfedb42052afaf452d2afb2239b38b12694b2f42`.

`N3finemapping-reference.tsv` contains susieR's PIPs and coefficients plus
fit metadata. Its SHA-256 is
`2f428cd896e364aa3f1bc2e2e3807d6583106b85724459bd9e74a9d68140c2b7`.
