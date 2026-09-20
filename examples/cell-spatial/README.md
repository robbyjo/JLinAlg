# Synthetic cell and tissue example

Generated deterministically by `tools/generate_cell_spatial_examples.py` (seed 814).
Eight independent donors have control/case labels; the alternate paired sample
sheet treats the eight samples as four complete before/after donor pairs.
This is a schema/implementation fixture, not a biological cohort.

- `counts.tsv`: 9,920 sparse integer entries, 320 observations, 31 features.
- `cells.tsv`: manual T/B labels, specimen/section IDs, anatomical strata,
  physical coordinates and a supplied distance to a hypothetical boundary.
- `samples.tsv`, `paired-samples.tsv`: independent/paired designs.
- `features.tsv`: complete measured panel and mitochondrial flags.
- `selected-features.tsv`, `gene-sets.tsv`: prespecified example hypotheses.

The single-cell and spatial vignettes contain executable commands. Optional
`jlinalg.example.yaml` illustrates shared configuration; copy it to `jlinalg.yaml`
only after adjusting the explicit paths for the chosen project directory.
