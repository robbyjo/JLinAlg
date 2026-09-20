# Development TODO

Last reviewed: 2026-09-20.

This inventory tracks open work. Completed implementation and audit details are
recorded in the [advanced-method validation report](docs/advanced-validation.md),
the [v0.3.0 remaining-code audit](docs/release-0.3.0-audit.md), the
[release notes](RELEASE_NOTES.md), and the
[regression inference validation report](docs/regression-inference-validation.md).

## Omics pipeline foundations — shared design and data contracts

The following roadmap records planned work, not implemented capabilities. Build
reusable analytical operations before composing complete pipelines. Initial entry
points should accept quantified matrices, supplied annotations and segmentation;
raw-read alignment, raw mass-spectrum processing and image segmentation are
separate future scope. LLM integration remains deferred.

- [ ] **Pipeline families.** Define profiles for bulk single-source omics;
  specialized RNA (miRNA, lincRNA/other lncRNA and extracellular RNA);
  multiple-source/multimodal omics; single-cell RNA/protein; spatial omics; and
  microbiome taxonomic/functional profiles.
  Distinguish cell-resolved spatial assays from mixed-cell spots. Evaluate coverage
  against representative studies rather than promising a fixed coverage percentage.
- [ ] **Required source and study metadata.** Record tissue/biofluid, anatomical
  site, compartment, preparation, assay, species, donor, specimen, condition,
  visit/time, batch and covariates. Preserve donor-to-specimen-to-section-to-cell/
  spot hierarchy and distinguish extracellular collection source from inferred
  tissue of origin. Reject unidentifiable condition/batch designs for inference.
- [ ] **Replication and source relationships.** Declare independent same-source
  replication, cross-source generalization, orthogonal assay validation and
  cross-modality corroboration separately. Record matched donors and cohort
  overlap; do not pool different sources without an explicit common estimand.
  Support shared, source-specific and source-by-condition effects.
- [ ] **Data representations and identities.** Preserve original counts/intensities,
  normalized measurements, embeddings and aggregates as separate labeled objects.
  Record feature identifiers, units, detection limits, ambiguous mappings and
  measured feature universes. Distinguish same-cell, same-donor and unrelated
  reference integration; never silently fabricate matched measurements.
- [ ] **Scalable interchange.** Add sparse/chunked matrix processing and validated
  import/export adapters for AnnData/MuData, SingleCellExperiment/QFeatures and
  SpatialData or equivalent explicit tables. Preserve assay layers, sample IDs,
  coordinates and metadata; benchmark memory and avoid implicit dense conversion.
- [ ] **Composable execution.** Describe each operation's inputs, supported designs,
  assumptions, outputs, diagnostics and resource limits. Add resumable execution,
  versioned manifests and cached artifacts using project-over-local configuration.
  Support entry from raw quantified, preprocessed or already annotated matrices.

## Multivariate and phenotype-family analyses

Existing [multivariate Gaussian regression](docs/vignettes/regression-families.md)
provides shared-design coefficients and residual covariance;
[multivariate MR](docs/vignettes/multivariate-mr.md),
[genomic factors](docs/vignettes/genomic-factor.md), SEM and ordinary mixed/GEE
models provide specialized components. They do not yet constitute a general
phenotype-family inference and reporting workflow.

- [ ] **Phenotype-family schema.** Declare continuous, binary, ordinal, count and
  time-to-event traits; units, transformations, assay/time, missingness, covariates,
  outcome/exposure/mediator roles and planned contrasts. Provide editable lipid,
  cardiovascular, lung, kidney and hematology examples; require explicit panel
  membership rather than assuming a universal set of risk factors.
- [ ] **Derived and redundant phenotypes.** Record ratios, totals/components,
  calculated biomarkers and deterministic transformations. Diagnose rank deficiency
  and near-singular covariance; define supported estimable contrasts without
  silently dropping traits or treating derived measures as independent evidence.
- [ ] **Explicit analysis modes.** Separate marginal per-trait analyses, joint
  outcome tests, mutually adjusted conditional models and latent/composite traits.
  Declare the model direction and estimand. A joint association must not imply an
  association with every trait; conditional and marginal effects are distinct.
- [ ] **Family-wide separate-model runner.** Use an appropriate model for each
  outcome, explicit covariate policies and common versus per-trait sample sets.
  Report effects, SEs, CIs, participant/event counts, exclusions and comparable
  diagnostics. Reuse existing model APIs without overstating their supported scope.
- [ ] **Continuous multivariate inference.** Extend Gaussian regression with full
  coefficient sampling covariance, joint Wald/F tests, MANOVA-style tests
  (including Pillai/Wilks where assumptions permit), general contrasts and CLI
  reporting. Validate degrees of freedom, covariance rank and small-sample behavior.
- [ ] **Correlated-score infrastructure.** Export and combine score/effect vectors
  with labeled full cross-trait sampling covariance, retained participant sets and
  overlap metadata. Distinguish phenotypic, residual, genetic and sampling-error
  correlations; do not substitute one for another.
- [ ] **Mixed-outcome joint tests.** First combine appropriately fitted marginal
  continuous/binary/count/ordinal models through validated correlated-score tests.
  Add supported joint likelihood models subsequently; survival outcomes require
  censoring-aware contributions. Do not put mixed outcomes into Gaussian MANOVA
  or present a p-value combination as a complete joint response model.
- [ ] **Multivariate mixed and marginal models.** Model cross-trait residual and
  random-effect covariance, relatedness, donor clusters, repeated visits and paired
  tissues. Add trait/source/time interactions, robust covariance where justified,
  identification diagnostics and computational limits. Existing univariate mixed
  models alone do not supply these capabilities.
- [ ] **Missing outcomes.** Start with explicit complete-case/per-trait policies
  and missingness reports. Extend validated likelihood or multiple-imputation
  inference to phenotype families, preserving cross-trait covariance and pooling
  joint tests appropriately. Document missingness assumptions and sensitivity.
- [ ] **Latent and composite phenotypes.** Add validated family-level PCA/factor
  and prespecified composite workflows, with scale/orientation, model-fit checks,
  uncertainty and held-out projection. Keep phenotypic factors distinct from
  genetic factors and avoid outcome-driven score selection on validation samples.
- [ ] **Multiplicity and post-omnibus follow-up.** Define testing families across
  molecular features, phenotype groups, traits and cell types. Reuse existing
  correction primitives with validated hierarchical/selective rules; filtering
  traits after an omnibus test does not make unadjusted follow-up tests valid.
- [ ] **Joint replication and meta-analysis.** Support per-trait replication and
  effect-vector pooling with full within-study sampling covariance, overlapping
  cohorts, missing traits, harmonized scales and heterogeneous source effects.
  State a common estimand and distinguish replication from cross-tissue support.
- [ ] **Family reports.** Export effects and covariance, joint and trait-specific
  tests, missingness/sample-overlap matrices, forest plots and effect heatmaps.
  Retain original-unit interpretations alongside explicitly standardized effects.

## Single-cell analyses

The [single-cell RNA workflow](docs/vignettes/single-cell.md) now implements
quantified/annotated-input QC, sparse count aggregation, independent/paired donor
state analysis through limma voom, sample-level relative abundance, prespecified
gene-set scores and exploratory PCA/k-means. See the
[validation report](docs/cell-spatial-validation.md). The following broad roadmap
items remain unchecked where their full assay/design scope is not implemented.

- [x] **Annotated replicated RNA delivery.** CLI, strict sparse counts and measured
  panel, complete sample/donor metadata, QC/exclusion reports, pseudobulk with cell
  counts and library sizes, independent/paired designs, batch/numeric covariates,
  limma version/source provenance, joint BH, tutorials and independent numerical
  checks. This does not implement arbitrary longitudinal models or other assays.
- [x] **Bounded RNA exploration and sample summaries.** Library normalization,
  variance-selected PCA/k-means, relative population abundance and prespecified
  mean-expression gene-set scores, with explicit models and limits.

- [ ] **Cell/sample QC.** Add assay-aware library/feature coverage and mitochondrial
  metrics, empty-droplet/ambient-RNA assessment, doublet detection and explicit
  filtering reports. Check donor representation and low-information populations;
  support annotated-input workflows without rerunning irrelevant preprocessing.
  RNA library/feature/mitochondrial metrics and imported exclusions are delivered;
  automatic empty-droplet, ambient-RNA and doublet assessment remains open.
- [ ] **Representation and clustering.** Add normalization, variable-feature
  selection, PCA, nearest-neighbor graphs, clustering and visualization adapters.
  Keep integration/embedding values separate from measurements used for inference;
  assess preservation of biology and batch/condition confounding.
  Library normalization, variance selection, PCA and k-means are delivered;
  neighbor-based clustering, integration and their biological validation remain open.
- [ ] **Annotation and markers.** Support marker discovery, reference mapping,
  hierarchical cell labels, confidence/unknown labels and manual annotation import.
  Preserve reference versions and distinguish cluster markers from between-condition
  differential effects; do not force novel states into an existing reference label.
- [ ] **Pseudobulk and differential state.** Aggregate raw RNA counts by biological
  sample/visit and cell type, retaining donor structure, cell counts and offsets.
  Fit donor-aware case-control, paired, longitudinal and covariate-adjusted models
  within supported designs. Report insufficient cells/donors instead of creating
  artificial replicates; add validated cell-level hierarchical alternatives later.
  Independent donors and complete two-condition pairs with numeric covariates and
  batch are delivered. Longitudinal and cell-level hierarchical extensions remain open.
- [ ] **Differential cellular abundance.** Analyze sample-level population counts
  with compositional/denominator-aware inference and biological replication.
  Add neighborhood differential-abundance adapters (for example Milo) for continuous
  states; distinguish relative composition from absolute abundance.
  Equal-sample transformed-proportion inference is delivered; compositional joint
  models, rare-population calibration and Milo adapters remain open.
- [ ] **Cell-specific phenotype families.** Connect each cell type/state to the
  separate and joint phenotype-family operations above, preserving donor covariance
  and multiplicity across genes/proteins, cell types and phenotypes.
- [ ] **Functional analysis and replication.** Add cell-specific pathway/gene-set
  and regulator activity analyses with assay-appropriate backgrounds, independent
  cohort replication and effect meta-analysis. Assess sensitivity to annotation,
  filtering and donor composition; do not count cell numbers as independent donors.
- [ ] **CITE-seq/antibody-derived tags.** Add antibody/background/control QC,
  normalization, RNA-protein integration (for example totalVI), annotation and
  modality-specific differential analysis/concordance. Preserve measured versus
  inferred values; do not feed denoised estimates into arbitrary tests as raw data.
- [ ] **Mass-spectrometry single-cell proteomics.** Add identification-quality and
  run/channel QC, peptide-to-protein provenance/aggregation, missingness and
  detection-limit handling, normalization and differential protein/pathway analysis.
  Evaluate scp/QFeatures adapters; do not reuse RNA count aggregation or silently
  impose blanket imputation on protein intensities.
- [ ] **Flow/mass cytometry and imaging proteins.** Support compensation/unmixing
  metadata where applicable, assay transformations, control/gating import, population
  annotation, differential marker intensity and abundance. Route imaging protein
  measurements through the spatial pipeline after segmentation/measurement QC.
- [ ] **scATAC and paired RNA-ATAC.** Add accessibility QC, peak quantification,
  annotation/integration, donor-aware differential accessibility, motif enrichment,
  TF accessibility scores and peak-gene links through suitable adapters (for example
  Signac). Overlay fine-mapped variants with genome-build and uncertainty checks;
  do not equate accessibility or a peak-gene correlation with causal regulation.
- [ ] **Trajectories and perturbations.** Add trajectory/branch and condition-specific
  dynamics analyses with root/ordering uncertainty; evaluate RNA velocity only for
  appropriate inputs and assumptions. Add perturbation-aware models with guide/
  treatment assignment, controls and donor effects. Pseudotime is not elapsed time.
- [ ] **Regulatory and communication hypotheses.** Extend regulatory/network
  adapters to appropriate cell-specific data, including ligand-receptor evidence,
  downstream responses and multi-modal links. Preserve edge types, context and
  uncertainty; expression compatibility alone does not establish signaling or
  causality. Track SCENIC/GRNBoost extensions under the network backlog below.
- [ ] **Cell-specific genetic follow-up.** Connect validated cell-type/state
  findings to genetic enrichment, suitable QTLs, fine-mapping, colocalization and
  MR when inputs support them. Preserve ancestry, tissue/state, instrument strength
  and overlap constraints; imported atlas evidence is not a matched donor cohort.

## Spatial transcriptomics and proteomics

Existing regional EWAS concerns genomic distance. The new
[physical spatial RNA workflow](docs/vignettes/spatial-analysis.md) adds the
following quantified-input capabilities. Unchecked items retain their broader
unimplemented scope; none imply completed protein/imaging/deconvolution support.

- [x] **Bounded physical graphs and pattern tests.** Explicit units and supplied
  coordinates, section/compartment restrictions, radius/union-kNN/adjacency,
  2D/3D distances, SVG section maps, Moran I/Geary C, stratified random-label tests,
  neighborhood enrichment/depletion and global testing families.
- [x] **Replicated spatial RNA summaries.** Supplied domain expression/pathway/
  abundance via the shared RNA commands, sample-level edge composition and linear
  distance-slope comparisons with independent/paired donor designs. CLI, tutorials,
  numerical fixtures and citations are delivered within documented limits.

- [ ] **Spatial input/QC contract.** Import cell-resolved, subcellular or mixed-spot
  measurements with physical units, coordinate transforms, masks, image references,
  segmentation IDs/quality and donor/specimen/section hierarchy. Preserve measured
  gene-panel limits and flag uncertain cell assignments or tissue/background mixing.
  Explicit RNA tables, panel, hierarchy, supplied metadata and exclusions are
  delivered; native image/mask/transform import and subcellular assignment remain open.
- [ ] **Spatial graphs.** Build radius, k-nearest-neighbor and appropriate
  adjacency graphs with explicit scales, tissue boundaries/gaps and 2D/3D support.
  Do not connect separate specimens or unregistered sections merely because their
  numeric coordinates overlap. Validate scale and graph-construction sensitivity.
  Native radius/kNN and supplied adjacency are delivered; atlas-scale indexing,
  automatic mask intersection and biological scale-sensitivity studies remain open.
- [ ] **Spatially variable features.** Add Moran's I/Geary's C and suitable spatial
  expression/activity models, covariate adjustment, multiple-testing correction
  and declared null models. Distinguish within-specimen spatial pattern tests from
  reproducible between-condition effects across independent specimens.
  Moran/Geary with declared random-label nulls are delivered; covariate-adjusted
  spatial processes and real-cohort calibration remain open.
- [ ] **Neighborhood relationships.** Add neighborhood enrichment, distance-dependent
  co-occurrence and spatial distribution statistics using justified permutations
  within exchangeable anatomical/sample strata. Preserve population abundance and
  boundary effects appropriate to each null; assess multiple spatial scales.
  Stratified label-pair enrichment/depletion is delivered; additional distance-
  distribution statistics and multi-scale validation remain open.
- [ ] **Domains and niches.** Identify and characterize compartments and local
  cell communities; support pathology annotations and marker/pathway summaries.
  Assess stability and cross-sample correspondence; separate exploratory domain
  discovery from confirmatory tests performed on the same defining features.
- [ ] **Mixed-spot deconvolution.** Add reference-based mapping/deconvolution
  adapters (for example cell2location), reference compatibility/QC, uncertainty
  and missing-reference-population handling. Do not report inferred cell mixtures
  or imputed transcripts as directly observed single-cell measurements.
- [ ] **Replicated spatial comparisons.** Compare domain expression, composition,
  niche abundance and spatial relationships across conditions with donor/section
  dependence and spatial autocorrelation represented in the model. Support paired
  samples and replication; extra spots/sections cannot manufacture donor replication.
  Independent/paired donor summary comparisons are delivered; richer section-level
  dependence models and public replication studies remain open.
- [ ] **Spatial gradients and interfaces.** Model distances to vessels, lesions,
  tumor margins or other annotated structures, with nonlinear effects, covariates,
  sample variation and uncertainty in boundaries. Test condition-by-distance effects.
  Sample-specific linear slopes and donor-level slope contrasts are delivered;
  nonlinear covariate-adjusted trends and uncertain interfaces remain open.
- [ ] **Spatial signaling.** Combine proximity, ligand-receptor evidence and
  downstream activity into explicit hypotheses; compare against appropriate
  abundance/anatomy-aware nulls. Proximity alone is not evidence of causal signaling.
- [ ] **Histology and multi-section integration.** Add image-feature/pathology
  integration, section registration, cross-slice domain alignment and later 3D
  reconstruction, preserving registration uncertainty and sample identities.
- [ ] **Spatial multimodal reporting.** Link RNA/protein/epigenetic measurements,
  cell states and phenotype families within known matching relationships. Export
  coordinate-aware maps, graphs, effects and uncertainty with reproducible scales;
  retain source annotations for downstream pathway/genetic/network follow-up.

## Microbiome analyses — beyond zero-inflated models

Extend the shared source/design and phenotype-family contracts to microbial
communities. Existing count, mixed and zero-inflated estimators are components,
not a complete microbiome workflow. Initial inputs are quantified feature tables,
taxonomy, functional annotations and optional phylogenies/load measurements;
raw-read denoising, assembly and taxonomic profiling remain external operations.

- [ ] **Assay and source contracts.** Distinguish 16S/ITS amplicons, shotgun
  metagenomics and metatranscriptomics; ASVs/OTUs, taxa, strains and gene/pathway
  profiles; bacterial, archaeal, fungal and viral measurements. Require body site/
  environmental source, donor/subject, visit, collection/storage, extraction and
  sequencing batch, assay/database version and relevant exposure covariates such
  as antibiotics or diet. Keep read counts, proportions and calibrated loads distinct.
- [ ] **Interchange and taxonomy/phylogeny alignment.** Import BIOM or explicit
  count tables, taxonomy mappings, sample metadata and Newick trees, with adapters
  to phyloseq/TreeSummarizedExperiment or equivalent objects. Validate tree tips,
  branch lengths and feature identity; preserve unclassified taxa and prevent
  double counting parents together with their descendants. Harmonize databases
  and taxonomic ranks before replication or cross-study aggregation.
- [ ] **QC and contamination.** Report library size, prevalence, abundance,
  sample coverage and sequencing/extraction effects. Add negative-control and
  concentration-aware contamination assessment (for example decontam), mock
  community checks and low-biomass sensitivity. Separate removal rules from
  biological absence; record exclusions and avoid outcome-selected filtering.
- [ ] **Compositional transformations and estimands.** Add explicit reference
  log-ratios, CLR/ILR and interpretable balances with stored denominators/bases.
  Diagnose singular CLR covariance; declare zero replacement/pseudocount policies
  and sensitivity. A sequencing zero need not mean biological absence, and changing
  relative abundance need not mean changing absolute abundance.
- [ ] **Absolute microbial load.** Integrate compatible qPCR/dPCR, spike-in or
  cell-count measurements with units, calibration and uncertainty. Distinguish
  measured load from model-inferred scale, accounting for relevant copy-number/
  assay biases. Do not claim absolute changes from uncalibrated proportions alone
  without declaring the identifying assumptions of the chosen model.
- [ ] **Alpha diversity.** Implement observed richness, Shannon/Simpson and Hill
  diversity, plus phylogenetic diversity when a valid tree exists. Report sampling
  depth/coverage, uncertainty and supported standardization/rarefaction options.
  Compare diversity using covariate- and donor-aware models; rarefaction must not
  silently become the input to every differential-abundance analysis.
- [ ] **Beta diversity and ordination.** Add Bray-Curtis, presence/absence Jaccard,
  Aitchison and weighted/unweighted UniFrac where their input assumptions hold.
  Support PCoA/NMDS and constrained ordination with explicit scaling, zero handling,
  negative-eigenvalue treatment and NMDS diagnostics. Separate visualization from
  inferential evidence; benchmark large distance matrices and memory limits.
- [ ] **Community-level association.** Add PERMANOVA/distance-based models with
  declared sequential versus marginal hypotheses, effect sizes, covariates and
  dispersion diagnostics. Support design-appropriate restricted permutations for
  paired/clustered/longitudinal studies, including between-subject exposures that
  cannot be tested by permuting only within subjects. Report permutation resolution
  and multiplicity; dispersion differences can complicate location interpretation.
- [ ] **Differential abundance and prevalence.** Add composition-aware per-feature
  comparisons, continuous exposures, multigroup contrasts and supported repeated
  measures. Evaluate ANCOM-BC2 and ALDEx2 adapters against independent fixtures;
  declare each method's estimand, assumptions and stochastic settings rather than
  treating them as interchangeable. Report effects, uncertainty, adjusted tests,
  reference/scale sensitivity and prevalence/detection separately from abundance.
- [ ] **Taxonomic/functional multiplicity.** Define feature universes and testing
  families across ranks, pathways, phenotypes and contrasts. Add validated
  hierarchical/group testing where useful; a taxon and its parent or overlapping
  pathways do not provide independent corroboration. Preserve nonsignificant and
  filtered features in audit outputs.
- [ ] **Longitudinal and intervention studies.** Model within-subject changes,
  treatment-by-time effects, irregular visits, temporal dependence and paired body
  sites. Add community turnover/stability and perturbation/recovery summaries with
  uncertainty. Separate person-level and within-person effects, and handle missing
  visits and time-varying exposures under stated assumptions.
- [ ] **Functional profiling.** Import measured metagenomic gene-family/pathway
  abundances (for example HUMAnN outputs), pathway coverage and taxon-stratified
  contributions. Extend to resistome/virulence annotations with reference versions.
  Keep amplicon-predicted function distinct from measured functional profiles and
  DNA functional potential distinct from RNA activity; account for abundance when
  interpreting matched metatranscriptomic changes.
- [ ] **Host-microbiome and phenotype-family integration.** Relate community
  summaries, balances, taxa and functions to lipid/lung/kidney/other phenotype
  families and matched host RNA/protein/metabolite data. Support separate and joint
  tests, donor matching, covariates and held-out validation; preserve direction and
  distinguish association from mechanistic mediation or causal effects.
- [ ] **Compositional networks.** Evaluate specialized sparse association methods
  such as SPIEC-EASI, stability selection, donor-aware resampling and differential
  network assessment. Do not apply ordinary correlations or the generic network
  command directly to proportions and label edges ecological interactions. Report
  sensitivity to filtering, zero handling and reference choices.
- [ ] **Prediction and reproducible signatures.** Add taxon/function/balance-based
  prediction with donor/site/cohort-grouped validation, training-only preprocessing
  and feature selection, calibration and external-cohort evaluation. Quantify
  stability and transportability; avoid presenting a universal dysbiosis score.
- [ ] **Community types and source/strain analyses.** Later assess community-state
  clustering with stability and continuous-gradient alternatives, microbial source
  tracking with unknown-source uncertainty, and supplied strain-resolution profiles
  for colonization/sharing analyses. Similarity alone does not establish transmission
  direction; source panels, longitudinal design and resolution constrain conclusions.
- [ ] **Causal and genetic follow-up.** Connect suitable microbial QTL/host genetic
  findings to existing MR, colocalization and mediation only when harmonization,
  instrument validity, composition, sample overlap and confounding assumptions are
  supported. Cross-sectional taxon-metabolite associations do not establish a pathway.
- [ ] **Replication and reporting.** Produce diversity/community summaries,
  differential-feature tables, taxonomic/functional maps and evidence-linked
  reports. Compare harmonized effects across cohorts and sources with explicit
  estimands, overlap and heterogeneity; retain database versions, transformations,
  controls, seeds and tested method boundaries in CLI manifests and vignettes.
- [ ] **Delivery and validation.** First deliver table/taxonomy/QC contracts,
  compositional transforms, diversity/ordination, community tests and one validated
  differential-abundance path; then longitudinal, functional and host-integration
  workflows, followed by network/source/strain extensions. Use independent vegan,
  ANCOM-BC2/ALDEx2 or other appropriate fixtures, public replicated cohorts and
  synthetic nulls varying depth, load, zeros, contamination, dispersion, taxonomy
  aggregation and subject dependence. Check effect recovery, null calibration,
  FDR/coverage, seed reproducibility and memory/runtime; method agreement alone
  does not establish correctness.

Method references and candidate adapter documentation:

- [Compositional reference frames](https://www.nature.com/articles/s41467-019-10656-5).
- [ANCOM-BC2: multigroup and repeated-measure analysis](https://www.nature.com/articles/s41592-023-02092-7).
- [ALDEx2 method and supported interfaces](https://bioconductor.org/packages/release/bioc/manuals/ALDEx2/man/ALDEx2.pdf).
- [vegan PERMANOVA and interpretation](https://vegandevs.github.io/vegan/reference/adonis.html).
- [decontam controls and contamination assessment](https://benjjneb.github.io/decontam/).
- [HUMAnN functional profiling](https://huttenhower.sph.harvard.edu/humann/).
- [SPIEC-EASI compositional network inference](https://pubmed.ncbi.nlm.nih.gov/25950956/).

## Delivery order and validation for the new omics workflows

- [ ] **Milestone 1 — shared contracts and phenotype families.** Deliver source/design
  schemas, separate-model execution, unified reports and multiplicity, then complete
  continuous multivariate inference. Include CLI, configuration examples and vignettes.
- [ ] **Milestone 2 — replicated scRNA-seq.** Deliver annotated-input and QC-to-analysis
  paths through donor-aware differential state/abundance, pathways and replication.
  Extend to CITE-seq; implement mass-spectrometry proteomics as a separate assay branch.
- [ ] **Milestone 3 — replicated spatial assays.** From supplied measurements,
  coordinates and segmentation/annotations, deliver graphs, neighborhood/domain
  analyses and sample-aware comparisons, with a mixed-spot deconvolution branch.
- [ ] **Milestone 4 — advanced designs.** Add mixed-outcome joint tests, multivariate
  mixed models and joint replication; then validated trajectories, perturbation,
  scATAC, regulatory/signaling and multi-modal/spatial extensions.
- [ ] **Microbiome delivery track.** Reuse Milestone 1's source/design and phenotype
  contracts, then follow the ordered microbiome delivery plan above. Shared count
  models do not substitute for compositional and community-level method validation.
- [ ] **Acceptance fixtures.** Validate coefficients/covariances and joint tests
  against independent R/Python implementations; test null calibration, interval
  coverage and multiplicity under correlated traits, redundant outcomes, mixed
  outcome types, unequal missingness, cohort overlap and repeated observations.
- [ ] **Single-cell/spatial acceptance studies.** Include multi-donor public examples
  and synthetic nulls with rare populations, composition changes, batch confounding,
  unequal cells per donor, missing protein measurements, tissue boundaries, spatial
  autocorrelation, multiple sections and imperfect references. Check donor-level
  split/replication behavior and fail unsupported population inference explicitly.
- [ ] **Publication-quality delivery.** Add focused API/CLI tests, source-versioned
  adapters, numerical and workload benchmarks, runnable tutorials/vignettes, primary
  references, diagnostics and reports for each supported workflow. Document tested
  limits and incomplete branches; validate scientific behavior as well as execution.

## xWAS analysis additions

Configuration, variant follow-up and all six network method categories are now
available in the [follow-up tutorials](docs/vignettes/network-followup.md).
Further scope remains explicit:

- [ ] Native indexed whole-genome annotation stores and provider-specific catalog
  discovery; current snapshots use streaming TSV lookup with explicit schema mapping.
- [ ] Additional regulatory adapters (SCENIC/SCENIC+, GRNBoost2), covariate/cluster-aware
  network inference, paired differential networks and larger sparse/out-of-core matrices.
- [ ] Full real-data validation of local VEP/ANNOVAR installations, CBC/CPLEX solver
  adapters and biological target prioritization; do not generalize synthetic checks.
- [ ] Integrated workflow scheduling, a cross-provider lockfile manager, and drug-target
  dossiers. LLM integration is deliberately deferred pending design decisions.

Implemented differential, regional, multiplicity and imputation workflows are
documented in [modern inference validation](docs/inference-workflows-validation.md).

The first four follow-up workflows (unpartitioned LDSC, genetically predicted
TWAS/PWAS, a single genetic factor with conditional SNP tests, and prediction
score training/application/evaluation) are implemented in the source build.
See [xWAS workflow validation and scope](docs/xwas-followup-validation.md).

- [ ] **Frozen confounder projection and count-aware adjustment.** Add
  fold-owned fit/freeze/apply artifacts for prediction, plus independently
  validated RUV and ComBat-Seq contracts. Do not refit preprocessing on held-
  out folds or apply Gaussian ComBat directly to raw counts.

Further extensions of the new workflows remain explicit scope boundaries:
partitioned/two-step/liability-scale LDSC and native summary munging; native
PredictDB/FUSION model adapters and rank-truncated multi-tissue inference;
multiple genomic factors, robust DWLS and propagation of measurement-model
uncertainty into SNP effects; LD-aware Bayesian PRS, logistic score training,
grouped CV and absolute-risk calibration. Current methods do not claim full
parity with LDSC, MetaXcan, GenomicSEM, or PRS software suites.

## Rare-variant meta-analysis — further scope

Implemented engineering, calibration, conditional/VT/heterogeneous tests, and
leave-out diagnostics are documented in the
[rare-variant tutorial](docs/rare-variant-meta-analysis.md) and
[validation report](docs/rare-meta-validation.md).

The Raremetal2 format/method assessment and quantitative model-metadata contract
are completed in the [trait-model assessment](docs/raremetal2-trait-models.md).
The remaining implementations are separate:

- [ ] **Allele-aware multiallelic summaries.** Implement producer-specific
  compressed covariance decoding and allele-aware union/group/conditioning
  keys, preserving both covariance axes and global row ordering. Validate
  same-site cross-ALT covariance, missing coverage and indexed/parallel equality.
  Old position-only covariance cannot resolve multiallelic groups.
- [ ] **Unbalanced-study exact method.** Specify the pooled estimand and required
  phenotype moments, residualization and covariance scaling before implementing
  Raremetal2 `--useExact`. Validate against independent participant-level fits
  and cross-products; this is separate from binary rare-case calibration.
- [ ] **Additional rare-meta trait models and calibration.** Implement a richer
  versioned model protocol and validated binary rare-case/joint group tails,
  followed by model-specific Poisson/survival support. Reuse the existing
  conditional-GWAS score infrastructure where appropriate; normal-only exports
  and cohort p-values do not supply calibrated rare-event group inference.

## Validation for new methods

Use checked-in fixtures from the underlying R estimators, independent
likelihood/derivative/covariance checks, and reproducible workload-specific
benchmarks. Zelig often wraps other packages and must not be the sole numerical
reference. Document estimator assumptions, convergence failures, and supported
scope; scenario contrasts require additional assumptions for causal interpretation.

## Conditional GWAS — further scope

The opt-in Gaussian/logistic/probit/Poisson/Cox aggregate export, cohort-side
conditioning refits, and summary-only importer are implemented; see the
[summary schema and validation](docs/conditional-gwas-summary.md).

- [ ] **Model-specific tail calibration.** Add validated binary rare-case and
  survival calibration, with the required score cumulants/risk-set information.
  Current exports explicitly provide normal tails and no calibrated p-value.
- [ ] **Additional cohort score models.** Add dedicated mixed, cluster-robust,
  stratified/repeated-subject Cox, grouped-binomial and other outcome exporters.
  Define cross-cohort overlap handling before pooling dependent study scores.

## Remaining estimator extensions

Implemented methods and their tested scope are documented in the
[estimator extension report](docs/estimator-extensions.md). Exact diffuse
historical smoothing now scales to long series, and ARIMA regression exposes
both conditional forecasts and forecasts with joint regression/dynamic
parameter uncertainty.

- [ ] **SEM scalability and additional test variants.** Extend joint ordinal
  integration beyond four observed ordinal responses per row; add multigroup
  WLSMV, DWLS-specific modification indices, robust nested-model comparisons,
  and robust CFI/TLI. Current multigroup and ordinal modification tests use the
  joint likelihood, with explicit equality constraints and MAR marginalization.
- [ ] **Discrete conditional quantile inference.** Extend beyond the iid
  marginal order-statistic interval to regression coefficients for responses
  with mass points. Cluster/HAC density-based regression inference still assumes
  continuous responses. Automatic bandwidths do not promise finite-sample or
  simultaneous coverage.

## Interpretation and numerical limits

These boundaries are explicit; this release does not claim complete parity
with every feature of lavaan, metafor, lme4, or the other R packages.

- The legacy MR conditional-GWAS interface assumes its documented Gaussian
  score model. Cohort model-specific scores/covariance support local inference;
  they do not reconstruct a complete nonlinear likelihood away from the null.
  New nonlinear conditioning sets require cohort-side refits or a richer
  likelihood evaluation protocol.
- Multidimensional integration supports binomial-logit and Poisson-log models.
  Tensor Hermite remains node bounded. The optional defensive randomized-QMC
  method uses empirical error estimates and bounded dense precision operations;
  it is not certified integration or a replacement for large sparse methods.
- Zero-inflated deterministic multi-start certification is relative to the
  configured starts, not a mathematical global optimum over the parameter
  continuum. The random-effect integral remains a single-mode Laplace
  approximation rather than a sum over separated conditional modes.
- Nonlinear models: floating-point quantization under very large response
  offsets can prevent a score certificate even when fitted SSE is near optimal.
  Such fits return nonconvergence. Mixed effects are additive Gaussian effects,
  not a general nonlinear random-parameter likelihood.
- Selective inference: unknown-noise or response-selected penalties, dependent
  sample splits, and non-Gaussian selection require different inferential methods.
- Exact quantile regression has iid residual-density, supplied-density,
  cluster, and HAC inference paths with distinct regularity assumptions.
  Response mass points use a separate marginal interval API; conditional
  quantile-regression inference at mass points remains open.
- Kernel set tests: the saddlepoint positive-mixture fallback is approximate,
  as is opt-in analytic SKAT-O. Default parametric SKAT-O simulation is exact
  for its Gaussian score-null model up to Monte Carlo error, not an exact
  finite-sample phenotype calibration.
- Numerical identification: extreme original-coordinate deficient SVD fits
  and inconclusive sparse covariance-rank checks reject rather than returning
  a numerically unsupported estimator or variance decomposition.
