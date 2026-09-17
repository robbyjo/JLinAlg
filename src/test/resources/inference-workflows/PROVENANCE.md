# Inference-workflow fixture provenance

`generate-reference.R` creates the synthetic continuous and count matrices and
freezes outputs from limma/voom, edgeR, and DESeq2. The seed is `17092026`; no
participant data or downloaded analysis data are used. `versions.tsv` records
the exact R, Bioconductor, and package versions used for the checked-in files.

Reproduce from the repository root after installing the packages into a local
library:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' src/test/resources/inference-workflows/generate-reference.R
```

The package comparisons test shared effect direction/scale and cross-feature
moderation behavior. JLinAlg's compact negative-binomial and voom contracts are
not byte-for-byte reimplementations of edgeR, DESeq2, or limma's full option
sets. Independent hand calculations separately test the implemented likelihood,
spatial covariance, multiple-testing, and pooling equations.
