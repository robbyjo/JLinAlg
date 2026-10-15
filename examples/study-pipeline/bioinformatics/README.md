# Synthetic bioinformatics example

Copy this directory to a scratch project before running it. Every value is invented for software testing; the SNP IDs and coordinates do not represent a real biological locus. The deliberately small regional tables are not a realistic fine-mapping benchmark.

Run `powershell -File run.ps1 -Jar /path/to/jlinalg.jar` or `sh run.sh /absolute/path/to/jlinalg.jar`. The scripts execute the study offline and verify its audit chain. Repeating the command reuses completed, verified operations. Use `--rerun` with `study auto` to request new attempts without deleting history.

The fixture exercises annotation, allele harmonization (the second outcome SNP is reversed), fine-mapping, colocalization and its shared-prior sensitivity, clumped MR, MR estimator sensitivity, leave-one-instrument-out inference, candidate evidence, diagnostics, reporting and provenance. Variant scoring is skipped because the primary results measure an omics feature, not a GWAS SNP.

The contract uses manual goals to keep this numerical demonstration offline and focused. Real studies can use `mode: auto` to add public-source discovery, protein/sequence investigation, literature search, enrichment and other applicable branches. Remove offline/deny-all restrictions only after reviewing disclosure requirements through the amendment workflow for an already-resolved project.

Inspect `index/index.html`, `SESSION_HISTORY.md`, `runs/*/operations/*/genetic-evidence.yaml`, the individual command logs and source snapshots. All original inputs are named under `datasets.files` so cache invalidation includes their contents. These scripts are reproducibility entry points; retain them alongside the contract and source-build JAR hash.

LD correlations can overshoot a mathematical bound by floating-point roundoff at perfect correlation. The producer and reader canonicalize only values beyond -1 or 1 within `1e-12`; materially invalid correlations, asymmetry and non-unit diagonals still fail validation. This numerical tolerance does not change the instrument-selection or LD-clumping thresholds. Within a genetic-evidence operation, the unchanged eligible-feature index is loaded once for all pairs; it is never reused across operations or projects.

Generated dependency manifests and evidence use a bounded 128 MiB reader, including the shared prioritization, diagnostics and reporting paths. Ordinary study contracts retain their stricter 16-million-code-point budget, and duplicate keys and YAML aliases remain rejected.
