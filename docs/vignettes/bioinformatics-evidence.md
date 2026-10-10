# Bioinformatics evidence workflow

JLinAlg can extend result-first omics studies beyond replication and meta-analysis into public-source discovery, molecular genetic inference, protein and reference-sequence investigation, literature retrieval, and transparent candidate prioritization. The deterministic planner retains control. Each operation saves its inputs, source citations, queries, responses, decisions and output hashes for the project index and research notebook.

This source-build feature does not make every biological question answerable automatically. In particular, choosing a valid molecular QTL, interpreting an unusual result, establishing instrument validity, or obtaining restricted data may require the scientist. Missing evidence is reported as a prerequisite, not fabricated and not interpreted as a null finding.

## Quick start

The smallest real-study entry point remains a study contract plus the input result files and a current executable JAR:

```sh
java -jar jlinalg-0.3.6.jar study auto study.yaml
java -jar jlinalg-0.3.6.jar study verify --project .
```

`analysis.mode` defaults to `auto`. Existing discovery/replication, meta-analysis, enrichment and diagnostic operations remain available. Auto mode now also plans `source_discovery` when genetic support is relevant, `bioinformatics_investigation`, and `literature_search`. They are independent branches: a retired QTL service does not prevent an available literature or protein investigation. Optional LLM literature interpretation is a separate `literature_support` operation.

Use `--offline` to prevent built-in network access. Repeat the same command after supplying an approved amendment or restoring connectivity. Completed operations are reused only when their inputs, implementation and outputs still match. Waiting operations are retried at operation granularity; a partial public search may repeat its earlier queries. `--rerun` requests new attempts while preserving history.

For a numerical demonstration with no internet or LLM, copy [the synthetic bioinformatics example](https://github.com/robbyjo/JLinAlg/tree/main/examples/study-pipeline/bioinformatics) to a scratch directory. Run `run.ps1 -Jar /path/to/jlinalg.jar` or `sh run.sh /absolute/path/to/jlinalg.jar`. Its invented regional statistics exercise the entire numerical follow-up branch. They are not evidence for a biological hypothesis or a realism benchmark.

## Keep the study contract small

Public investigation uses dataset annotations and the hypothesis outcome. Add optional search controls only when needed:

```yaml
datasets:
  - id: blood-results
    type: transcriptomics
    stage: association_results
    species: human
    tissue: whole_blood
    genome_build: GRCh38
    annotation_release: GENCODE_v44
    files:
      results: results.tsv
      annotation: probe-annotation.tsv

hypotheses:
  - id: H1
    description: Blood molecular associations with sleep duration
    outcome: sleep duration
    inputs: {discovery: [blood-results]}
    analysis:
      mode: auto
      bioinformatics:
        phenotype_terms: [sleep duration]
        relevant_tissues: [whole blood, brain]
        max_features: 10
        fdr: 0.05
        window_bp: 500000
```

This is a fragment to merge into a complete contract with `schema_version` and `project`. For an already-resolved project, make the change through `amendment-<ID>.yaml` and scientist approval rather than editing its base contract. The [study schema](../schemas/study-v1.schema.json) describes the optional settings.

Selection uses BH across the full valid tested family in the available meta-analysis, or separately within each primary table when no meta-analysis is available. It does not recalculate FDR on only the top hits. Primary-table selection is exploratory, not pooled evidence or an independent replication claim. Conflicting or ambiguous annotations are excluded. The feature limit defaults to 10 and is capped at 100; search budgets are visible in the result records.

## Feature annotation

Local annotation tables and GTFs take priority. Declare the coordinate convention for tables. Internal annotated intervals are 1-based and closed. Unknown probe IDs are never treated as gene symbols simply because they look plausible.

For human Ensembl features with an explicit modern GENCODE release and GRCh38 build, JLinAlg can retrieve the release-specific GTF automatically, cache it, and record the source URL and checksum. This download requests a public reference release, not the project's feature list. It does not guess a genome build, lift coordinates, or download an arbitrary platform manifest. Other platforms and unsupported releases need a local `files.annotation` table/GTF.

Annotation is not a claim that the feature is causal, that an sQTL exists, or that a transcript maps uniquely to a protein isoform. The generated scientist questions retain these distinctions.

## Sources and scientific relevance

| Source | Implemented use | Scientific boundary |
| --- | --- | --- |
| GWAS Catalog REST v2 and harmonized FTP | Primary or related-trait study discovery by public phenotype label, followed by local regional extraction from selected complete harmonized GRCh38 files | REST associations are curated lead associations, not full regional results. Automatic use requires a downloadable, ancestry-matched, at-least-100,000-variant study with sorted harmonized data and a tabix index. The complete selected file and index are installed in the shared cache and verified before local extraction. |
| Pan-UK Biobank | Automatic exposure-instrument lookup for the public BMI field 21001 EUR summary statistics; other traits remain discoverable candidates until a provider adapter or reviewed custom source is declared | The complete indexed source is installed in the shared cache. Its native build is GRCh37; JLinAlg projects biallelic exposure alleles locally with the installed hg38-to-hg19 chain, verifies exact Pan-UKB alleles, and emits an audited GRCh38-coordinate adapter table. It does not lift positions by offset or treat a UKB-plus-external meta-analysis as independent of every cohort. |
| eQTL Catalogue | Audited dataset discovery request, with a documented bulk-access fallback | An unavailable/retired API is reported as such. The measured feature and molecular assay must match; a gene eQTL is not automatically a splicing QTL. |
| GTEx | Versioned tissue-availability metadata, filtered using declared study or follow-up tissues | Tissue availability is not evidence of a QTL for a particular feature. Full regional QTL input still needs selection and provenance. |
| UniProt | Reviewed gene/species protein candidates | Multiple matches remain ambiguous; absence of reviewed entries does not mean no protein exists. |
| AlphaFold DB | Predicted structure identifiers, versions, confidence metrics and model links | Prediction is not isoform equivalence, binding, druggability, safety or efficacy evidence. |
| UCSC | Build-specific reference sequence and browser links for annotated intervals | Forward genomic sequence, not patient sequence or a reconstructed spliced transcript. |
| Europe PMC | Bounded field/phenotype and gene/phenotype searches, article identifiers, bibliographic metadata and a deduplicated reference catalog | Candidate papers require reading and appraisal. This is not a systematic review or automatic confirmation of their conclusions. |

The same framework applies to different phenotypes, cohorts and tissues. For example, a blood/sleep study may declare brain as a follow-up tissue with a documented rationale. JLinAlg does not infer that all available brain QTLs are appropriate. Related phenotype suggestions from an LLM remain proposals until reviewed and incorporated into the study.

The built-in GWAS Catalog/eQTL Catalogue/GTEx discovery branch requires explicitly human datasets. It does not substitute human genetic resources for another organism. Nonhuman studies need species-matched genetic-source adapters; the reviewed protein lookup currently supports human and mouse.

Searches distinguish `no_matches`, `search_limited`, `not_found`, `endpoint_retired`, `access_denied`, `rate_limited`, `unavailable`, `offline_policy`, and `policy_denied`. HTTP 404/410, malformed responses and connectivity failures never become a scientific “no associations exist” conclusion. Search defaults are three pages per strategy, 20 regional catalogue SNPs per feature, 20 literature records per page, and a hard 200-request budget per operation. A search hitting its budget is not exhaustive. UCSC requests are limited to one per second and sequences to 100 kb by default.

## From a source candidate to genetic inference

`genetic_pairs` is the legacy/custom-source name for an analysis recipe, not a biological object. One entry binds an exact molecular exposure table, phenotype outcome table, LD reference, region, ancestry, build, sample-size and overlap declaration. It remains necessary for arbitrary sources and nonstandard designs, but supported public MR sources no longer require a hand-written entry per transcript.

For an FHS irQTL exposure and BMI outcome, declare the concise source policy once:

```yaml
analysis:
  goals: [mendelian_randomization]
  bioinformatics:
    mr:
      exposure_dataset: fhs-bmi
      exposure_ancestry: European
      outcome_provider: auto
      outcome_trait: 21001
      outcome_ancestry: EUR
      ld_provider: 1000g-30x-grch38
      ld_population: EUR
```

The actual irQTL analysis ancestry is deliberately required; JLinAlg does not infer it from the FHS cohort name. For BMI/EUR, `auto` selects the public Pan-UK Biobank field 21001 EUR analysis and its recorded sample size. A non-UKB exposure cohort plus the UKB-only outcome establishes the named-cohort basis for `sample_overlap: none`; an exposure that may include UK Biobank is rejected. The policy expands the reviewed irQTL candidate file into per-feature recipes, records every local build projection and Pan-UKB lookup, and uses the GRCh38 LD provider below. The operation still stops if fewer than two strong, matched and LD-independent instruments remain.

For colocalization, JLinAlg can use GWAS Catalog REST v2 to find studies by public phenotype label and the Catalog's harmonized FTP collection for complete summary-statistic files. The REST API itself exposes metadata and curated lead associations, not arbitrary full-summary-statistic rows. JLinAlg verifies the selected accession's harmonized metadata, GRCh38 build, sorted status and tabix index, installs the complete file and index in the shared cache, verifies the Catalog MD5 when available, and performs regional extraction locally.

### Related-phenotype colocalization policy

Related phenotypes may be scientist-declared, or an optional LLM advisor can propose up to 20 search labels. The advisor prompt requests direct consequences, biomarkers, intermediate traits, common comorbidities and physiological domains, with one to five literature identifiers per phenotype. Every suggestion is stored as `advisor_proposal_unreviewed`; its literature still needs independent verification. The suggestion authorizes only bounded public-source discovery. It does not amend the hypothesis, approve a GWAS, or establish causality.

```yaml
assistance:
  enabled: true
  provider: openai
  model: APPROVED_MODEL
  related_phenotypes: true

analysis:
  goals: [colocalization]
  bioinformatics:
    colocalization:
      exposure_dataset: fhs-bmi
      exposure_ancestry: European
      outcome_provider: gwas-catalog
      outcome_ancestry: European
      related_phenotypes: [type 2 diabetes, glycated hemoglobin, systolic blood pressure, triglycerides]
      max_related_phenotypes: 20
      max_colocalization_traits: 5
      allow_partial_regions: true
      minimum_exposure_overlap: 0.5
      ld_provider: 1000g-30x-grch38
      ld_population: EUR
      sample_overlap: unknown
```

The primary phenotype and every declared/advisor-proposed label are searched separately. Default discovery requests contain public phenotype labels only; selected genes, variants and regions are not submitted. `remote_feature_queries: true` is an explicit opt-in for feature-context searches. For each label, the bounded automatic selector prefers the largest ancestry-matched study among the retrieved exact-trait GWAS Catalog pages when it has downloadable full summary statistics, at least 100,000 variants and no obvious CNV-only label. This is a reproducible availability/sample-size heuristic, not a claim that the chosen phenotype is the most biologically important. By default at most five traits proceed to colocalization; up to 20 can be searched or explicitly analyzed.

The irQTL export may be significance-filtered or otherwise incomplete. With `allow_partial_regions: true`, JLinAlg still performs harmonization, ancestry-matched GRCh38 LD extraction, fine-mapping mechanics and prior-sensitivity colocalization when at least two usable variants remain and the declared fraction of exported QTL variants matches. The result is labeled `exploratory_partial_region_major_limitation`, not conventional complete-region colocalization. The operation records exposure and outcome match fractions independently; it never calls missing database rows negative evidence.

Every generated report/manuscript must retain this caveat: partial-data posterior probabilities are conditional on the observed selected variants; missing variants, significance filtering, unequal coverage, winner's curse and unknown sample overlap can materially bias them, and a null result is not evidence against a shared causal variant. Set `allow_partial_regions: false` when the scientific plan requires conventional full-region colocalization.

The regional adapter runs after a scientist approves an exact exposure/outcome policy or pair. Its declarations can be short, but the required facts cannot be guessed: molecular feature, phenotype, assay, tissue, build, ancestry, sample size, source citation, regional coverage and sample overlap. See the complete checked-in example for custom-source YAML.

Each `analysis.bioinformatics.genetic_pairs` entry has `id`, `feature_id`, `rationale`, `approved`, a 1-based `region`, `exposure`, `outcome`, `ld`, `sample_overlap`, and `sample_overlap_basis`. MR additionally requires `instrument_validity_rationale`. These are assertions to audit, not facts proved by JLinAlg. Approval does not replace the source paper or scientific appraisal.

Local summary files must be declared under an input dataset's `files` mapping and referenced by `dataset_id` plus `file_role`. This includes QTL, GWAS, LD and LD-allele files, so changing any of them invalidates reuse. Source records require `source_id` and `citation`.

Canonical regional summary columns are:

```text
variant_id chromosome position effect_allele other_allele beta se p_value
```

Use tab-separated or comma-separated local files, with 1-based positions and biallelic SNPs. A QTL `feature_id` column, if present, is filtered to the approved measured feature. Without it, the file must already contain only that feature's regional statistics. Do not put patient genotypes or individual phenotype records in these files. Declare `complete_region: true` only for full regional results, never a significant-only table.

A source may instead specify a public HTTPS `url`, pinned `sha256`, `license`, `source_id` and `citation`, together with its scientific metadata. JLinAlg downloads and verifies it into a project cache and retains a source snapshot. This path currently expects canonical, uncompressed regional TSV data and caps downloads at 512 MiB; it does not guess arbitrary columns or decompress multi-gigabyte genome-wide releases. Supported remote hosts are restricted to the built-in public evidence providers; redirects are refused. Other sources, including restricted or differently formatted UK Biobank/GTEx releases, must be acquired and converted locally with the conversion script and citations retained. A downloaded source is not automatically approved merely because its checksum matches.

The adapter checks the finding's identity and location, same build and declared ancestry, identical variant coordinates, allele direction and regional overlap. Conventional complete-region analysis defaults to at least 80% of the larger source's variant set. Exploratory partial-region analysis instead checks the matched fraction of the exported QTL variants and reports the much smaller outcome-region fraction separately. Palindromic SNPs are conservatively excluded rather than guessed using absent allele frequencies. Every exclusion and sign flip is recorded. MR requires explicitly nonoverlapping samples. Partial colocalization can preserve `sample_overlap: unknown`, but that uncertainty is appended to its mandatory major caveat.

### LD and automatic reference downloads

For a supplied LD matrix, rows and columns must have matching variant IDs. An `alleles_file_role` table supplies `variant_id`, `effect_allele`, and `other_allele` defining the matrix's signed orientation. Matrix symmetry, bounds, diagonal and positive semidefiniteness are checked; invalid matrices are not silently replaced or regularized.

Alternatively set `ld: {population: EUR}` (or AFR/AMR/EAS/SAS) and allow reference downloads. The existing pinned 1000 Genomes Phase 3 installer fetches a missing panel. The pipeline extracts signed Pearson A1-dosage correlation among founders, checks at least 95% observed calls per SNP and mean-imputes remaining missing reference calls. This is different from the haplotype r-squared used by the standalone clumping command. Alleles, reference coordinates and panel hashes are retained.

The downloadable bundled panel is GRCh37. For GRCh38 MR, `ld_provider: 1000g-30x-grch38` installs the complete official 30x high-coverage phased 1000 Genomes VCF and tabix index for each required chromosome in the shared cache. It then restricts genotypes locally to the requested AFR/AMR/EAS/EUR/SAS superpopulation in the official 2,504-sample Phase 3 panel, thereby excluding the additional relatives in the 3,202-sample callset; applies call-rate and MAF checks; and records source and derived-file hashes, exclusions, alleles and the resulting matrix hash. Reuse works under `--offline`; the tradeoff is explicit chromosome-scale storage.

JLinAlg will not apply GRCh37 LD to mismatching coordinates. Missing, monomorphic or incompatible reference variants remain actionable blockers. Population is never guessed; the source ancestry must agree with the selected panel. The public source definitions are documented by the [1000 Genomes 30x GRCh38 collection](https://www.internationalgenome.org/data-portal/data-collection/1000genomes_30x) and [Pan-UK Biobank per-phenotype files](https://pan.ukbb.broadinstitute.org/docs/per-phenotype-files/index.html).

### Statistical outputs

SuSiE fine-mapping runs on the harmonized overlapping SNP set for both traits and must converge. Colocalization uses those SuSiE effects, with a second shared-prior setting for sensitivity. Curated nearby lead SNPs or omics association coefficients alone are never presented as conventional colocalization. When a reviewed policy explicitly permits incomplete QTL data, the same computation may run as exploratory partial-data colocalization with the evidence-grade and mandatory caveat above.

MR selects exposure instruments using p ≤ 5e-8 and F ≥ 10 by default, then LD-prunes at r-squared ≤ 0.001. It requires at least two remaining instruments for this automatic IVW branch. Outputs include fixed/random IVW, and, with at least three instruments, Egger, weighted median and leave-one-instrument-out estimates. Weak instruments, pleiotropy, selection bias, tissue specificity and exclusion-restriction assumptions still need scientific appraisal. Numerical success is not proof of a causal drug target.

The current in-memory regional limit is 5,000 variants and the window cannot exceed 10 Mb. Multi-trait, multi-tissue and multiple-region analyses are separate declared pairs with separate provenance. No additional across-pair multiplicity adjustment is silently implied by their individual estimates.

## Questions and optional AI

Built-in missing-input questions are recorded in `QUESTIONS.md`, the question inbox and `SESSION_HISTORY.md`. Use `study questions --out answers.yaml`, then `study answer answers.yaml --actor-id PERSON --check` and `--record`. Answers record statements and may propose amendments; they do not grant execution or transmission permission.

Opt-in OpenAI/Gemini advice can produce a dedicated, typed list of up to 20 related phenotypes with literature identifiers, or separately propose cited QTL sources, up to three scientist questions, and bounded code text. Each proposal records the provider, returned model, prompt/response, retrieval citations and decision. Suggested phenotype labels may feed public GWAS discovery but remain explicitly unreviewed; they do not change the hypothesis or certify their references. Suggested questions are explicitly labelled as LLM suggestions. New scientific directions require review and amendments. Only registered analysis policies execute automatically.

Generated code is saved as a hashed `.txt` proposal and never executed by the proposal handler. Arbitrary generated-code execution, container/OS isolation, general pipeline local-model adapters, MCP and provider-specific OAuth login/refresh remain future integration work. Browser handoff is the default assistance transport when assistance is explicitly enabled. An optional OpenAI-compatible adapter accepts an API key or provider-issued OAuth access token through an environment variable; it does not borrow website sessions. These other capabilities are not made safe merely by adding a prompt or asking the scientist to approve a string.

## Privacy and offline execution

Public evidence searches transmit only constructed queries containing public identifiers, phenotype labels or reference intervals. No dataset file is uploaded by these adapters. Nevertheless, a gene/phenotype query can disclose unpublished research interests. The same deny-all policy is checked before sending these project-derived queries as before LLM requests. Additional per-request disclosure conditions fail closed when the adapter cannot satisfy them.

Use `project.data_policy.remote_transmission: deny_all` to prevent such project-derived requests; use `--offline` or `resources.network: offline` to prevent built-in network requests altogether. Fixed public-reference downloads keyed only by a declared reference release may still be permitted by deny-all, since they do not contain study findings; offline prevents those as well. An already verified local cache can still support analysis offline.

These are application-level controls, not a guarantee that unrelated software or an arbitrary user-launched script cannot transmit data. Strong non-egress assurance still requires an externally enforced network boundary, restricted filesystem access, validated disclosure bundles and institutional review. Review each result/index before giving it to an LLM; aggregate does not mean nonsensitive or automatically anonymous.

## Results and evidence packaging

`index/index.html` and its JSON counterpart show available findings, source statuses, citations and required inputs. Candidate tables keep association ranking, per-feature contextual evidence and operation availability separate. A structure or paper match is not counted as causal validation. All statistical claims should cite the result ID, source snapshot, command and hashed output—not just the dashboard.

To produce a local audit snapshot:

```sh
java -jar jlinalg-0.3.6.jar study export --project . --out evidence.zip --acknowledge-sensitive
```

The archive includes the contract, ledger, notebook, questions/answers, plans, amendments, operation artifacts and index. `EVIDENCE_BUNDLE.json` lists checksums and omissions. Original dataset files outside those evidence directories, reference caches and software binaries are not copied automatically. Intermediate outputs, answer text and LLM interactions may still be sensitive. External preprocessing and unrecorded software remain explicit audit gaps; retain their scripts, versions and input/output hashes.

The command does not upload the ZIP or assert journal/FDA readiness, anonymity, identity authentication or regulatory compliance. It refuses to overwrite an existing archive. The project records the export and its hash after the snapshot has been written.

## Validation and remaining boundaries

Automated tests use synthetic regional data and injected provider replies: successful inference, allele reversal, build mismatch, invalid LD, weak instruments, source unavailability, zero matches, schema changes, disclosure denial, proposal quarantine and archive integrity. Public service response shapes were also inspected using generic public identifiers; no private research result was transmitted and no paid LLM inference was used for that inspection.

Automatic selection and conversion of arbitrary genome-wide GWAS/QTL releases is not implemented. Catalogue discovery and exact, reviewed regional acquisition are distinct steps. A retired service or missing exact splicing-QTL mapping cannot be repaired by treating a gene eQTL as equivalent. The source index and scientist questions explain the next required input. These boundaries apply across projects, not just Splicing-BMI.

The generated references below cover the methods and public resources. Individual retrieved study/structure/article records retain their own citations and query snapshots; cite those sources in a manuscript alongside JLinAlg.

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [NHGRI-EBI GWAS Catalog (2026) — GWAS Catalog REST API v2 documentation](../CITATIONS.md#gwas-catalog-api-v2)
- [NHGRI-EBI GWAS Catalog (2026) — GWAS Catalog summary statistics documentation](../CITATIONS.md#gwas-catalog-summary-statistics)
- [eQTL Catalogue (2026) — eQTL Catalogue data access](../CITATIONS.md#eqtl-catalogue-data-access)
- [GTEx Consortium (2020) — The GTEx Consortium atlas of genetic regulatory effects across human tissues](../CITATIONS.md#gtex-atlas-2020)
- [UniProt Consortium (2026) — UniProt API queries](../CITATIONS.md#uniprot-api-documentation)
- [M. Varadi et al. (2022) — AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models](../CITATIONS.md#varadi-alphafold-2022)
- [UCSC Genome Browser (2026) — UCSC Genome Browser REST API](../CITATIONS.md#ucsc-api-documentation)
- [Europe PMC (2026) — Europe PMC RESTful Web Service](../CITATIONS.md#europe-pmc-rest)
- [Adam Frankish et al. (2021) — GENCODE 2021](../CITATIONS.md#frankish-gencode-2021) — [PMID: 33270111](https://pubmed.ncbi.nlm.nih.gov/33270111/) · [PMCID: PMC7778937](https://pmc.ncbi.nlm.nih.gov/articles/PMC7778937/)
- [1000 Genomes Project Consortium (2015) — A global reference for human genetic variation](../CITATIONS.md#thousand-genomes-2015) — [PMID: 26432245](https://pubmed.ncbi.nlm.nih.gov/26432245/) · [PMCID: PMC4750478](https://pmc.ncbi.nlm.nih.gov/articles/PMC4750478/)
- [Gao Wang, Abhishek Sarkar, Peter Carbonetto, and Matthew Stephens (2020) — A simple new approach to variable selection in regression, with application to genetic fine mapping](../CITATIONS.md#wang-susie-2020) — [PMID: 37220626](https://pubmed.ncbi.nlm.nih.gov/37220626/) · [PMCID: PMC10201948](https://pmc.ncbi.nlm.nih.gov/articles/PMC10201948/)
- [Chris Wallace (2021) — A more accurate method for colocalisation analysis allowing for multiple causal variants](../CITATIONS.md#wallace-coloc-2021) — [PMID: 34587156](https://pubmed.ncbi.nlm.nih.gov/34587156/) · [PMCID: PMC8504726](https://pmc.ncbi.nlm.nih.gov/articles/PMC8504726/)
- [Jack Bowden, George Davey Smith, and Stephen Burgess (2015) — Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression](../CITATIONS.md#bowden-mr-egger-2015) — [PMID: 26050253](https://pubmed.ncbi.nlm.nih.gov/26050253/) · [PMCID: PMC4469799](https://pmc.ncbi.nlm.nih.gov/articles/PMC4469799/)
- [Jack Bowden et al. (2016) — Consistent estimation in Mendelian randomization with some invalid instruments using a weighted median estimator](../CITATIONS.md#bowden-weighted-median-2016) — [PMID: 27061298](https://pubmed.ncbi.nlm.nih.gov/27061298/) · [PMCID: PMC4849733](https://pmc.ncbi.nlm.nih.gov/articles/PMC4849733/)
- [Geir Kjetil Sandve et al. (2013) — Ten Simple Rules for Reproducible Computational Research](../CITATIONS.md#sandve-reproducibility-2013)
- [W3C Provenance Working Group (2013) — PROV-O: The PROV Ontology](../CITATIONS.md#w3c-prov-o-2013)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
## Publication reporting and local participant summaries

Evidence acquisition must explain **why the source addresses the hypothesis**,
not merely log a URL. Public requests now retain a scientific rationale,
hypothesis, operation, citation, request/response audit and retrieval status.
Cached uses retain their own rationale. Regional genetic sources have a separate
acquisition-intent record, including local sources. Ontology and reference
downloads explain their role in annotation, enrichment or genetic follow-up.

In automatic mode, a dataset with `files.phenotype` schedules the local
`participant_characteristics` operation. It does not transmit participant data.
For example, add the following to an existing dataset:

```yaml
files:
  results: results.csv
  phenotype: phenotype.csv
  analysis_members: analysis-members.csv
columns:
  characteristics:
    id: participant_id
    membership_role: analysis_members
    membership_id_column: participant_id
    population_rationale: Exact sample identifiers exported from the fitted model design matrix.
    age: visit_age
    main_phenotype: BMI
    sex: sex
    male_value: 1
    female_value: 2
    continuous: {White blood cell count: WBC}
    categorical: {Sequencing batch: batch}
    units: {Age: years, Main phenotype: kg/m2}
```

These sex codes are an example, not a default. Without a membership file, an
explicit `analysis_population_confirmed: true` and `population_rationale` may
declare that the phenotype file already represents the analysis sample. The
declaration is recorded, not independently proven. Otherwise JLinAlg computes
**provisional phenotype-file summaries** and requests verification; it does not
equate the number of phenotype records with association-analysis N. Repeated
participant IDs require explicit visit selection. Optional `filters` map columns
to lists of allowed values; they never execute code. Age at the assay visit must
not silently become baseline age. Single-sex cohorts may use `sex_constant:
female` or `male` with a documented `sex_basis`.

Continuous summaries use available-case means and sample standard deviations
(n−1), reporting nonmissing n and missing counts. Sex percentages retain their
denominator and unmapped/missing counts. No imputation is performed. With
feature-specific missingness, provide an appropriately defined primary analysis
population and explain the varying feature-specific N separately.

Every report operation now writes a local `publication-draft.html`,
`publication-supplement.html` and `publication-claims.yaml`, without requiring an
LLM. Results begin with Table 1. Methods use `datasets[].analysis_description`
fields `model`, `outcome`, `predictor`, `covariates`, `technical_factors`,
`relatedness`, `transformations`, `sample_selection` and `source`. Unknown fields
remain explicit requests, never reconstructed from a filename or script default.
Hypotheses are stated using their descriptions, not internal IDs. Common
abbreviations are expanded at first text use; additional definitions can be
declared in `hypotheses[].metadata.abbreviations`.

The deterministic draft presents Results in the same order as the applicable
Methods, after participant characteristics. Impact ranking belongs in Discussion,
not Results. A scientist or advisory LLM selects up to three distinct principal
findings and explains their importance using effect size, uncertainty,
reproducibility, biological relevance and contrary evidence. JLinAlg does not
assign scientific importance from operation type or the smallest P-value. Missing
interpretations remain explicit writing needs, not invented findings. Null
results and major limitations remain visible; the supplement retains all evidence
layers. Publication priority is not drug-target validation. Each result links to
its source manifest and checksums.

## Publication precision and combined participant characteristics

Publication HTML displays estimates and percentages to two decimal places, with
small nonzero values in scientific notation (for example, `8.01e-21`, never
`0.00` for that value). Integer counts remain integers. This is display rounding,
not truncation of the underlying evidence: result CSV files, YAML evidence and
numerical claim bindings retain their original precision. The same formatter is
available to native study figure-label generators; existing external figures
are not rewritten by manuscript revision.

Generated study YAML can be larger than SnakeYAML's general-purpose default when
many regional evidence pairs and exclusion audits are retained. Downstream study
operations use the same safe bounded reader as the project engine: duplicate
keys and aliases are rejected, nesting is capped, and the document limit is 16
million code points. The larger bound prevents valid generated evidence from
breaking prioritization and reporting without accepting unbounded input.

Table 1 places characteristics in rows and cohorts in columns. Multiple
hypotheses can share a cohort column only when phenotype-file hashes, participant
selection fingerprints, membership hashes, N, sex coding, units and overlapping
summaries agree. Different populations remain separate columns; JLinAlg does not
sum repeated cohort counts. Reduced available-case n is shown only where it
differs from N. Supplemental Table 1 retains the per-dataset details and missing
counts, alongside a full-precision `participant-characteristics.csv` export.

Explicit binary definitions are useful when several outcomes share a descriptive
table. For example, the following uses the investigator-specified definition
**BMI > 30**, not BMI ≥ 30:

```yaml
columns:
  characteristics:
    id: participant_id
    age: Age
    sex: Sex
    male_value: 1
    female_value: 2
    sex_basis: Investigator-confirmed data dictionary convention.
    main_phenotype: BMI
    continuous: {BMI: BMI}
    units: {Age: years, BMI: kg/m2}
    derived_binary:
      Obesity:
        column: BMI
        operator: gt
        threshold: 30
        definition: Obesity is defined as BMI > 30 kg/m2.
    analysis_population_confirmed: false
```

The definition is configurable; it is not a global obesity default. Missing BMI
is excluded from the binary denominator, not counted as non-obese. Operators
`gt`, `ge`, `lt` and `le` require an explicit finite threshold and definition.
Sex codes are never inferred from the presence of 1 and 2. Declare their meanings
from a dictionary or scientist statement, and document the basis. A descriptive
definition does not establish the definition used in previously fitted models.

To refresh demographics without refitting unrelated association analyses, run a
separate native study with the same hypothesis IDs, phenotype files and updated
`characteristics` settings. Set `analysis: {mode: manual, goals: []}` to request
only the baseline validation, participant summaries, diagnostics and report;
omit unrelated irQTL inputs and disable assistance/network. Retain its ledger.
Preserve the original `columns.validation` and `analysis.harmonization` settings:
reducing the requested analysis goals does not authorize changing the existing
invalid-row or common-feature policy.
Then add this explicit, hash-pinned reference to the manuscript revision:

```yaml
participant_project:
  path: /absolute/path/to/participant-summary-project
  ledger_sha256: EXACT_SHA256_OF_ledger/events.jsonl
  rationale: Update descriptive coding and outcome summaries as requested by the scientist.
```

JLinAlg verifies both projects and imports only current participant-characteristic
operations matching existing hypothesis IDs. Association, enrichment, genetic
and literature evidence still comes from the original project. The revision
preserves the refresh contract and ledger, source manifests, hashes and rationale.
Keep both projects with the manuscript. An unverified analysis population remains
provisional after this refresh; formatting cannot establish missing membership.

## Local isoform-ratio QTL review

A dataset may declare a local, optionally bzip2/gzip-compressed irQTL result file:

```yaml
files:
  irqtl: irqtl-cis-results.csv.bz2
columns:
  irqtl:
    genome_build: GRCh38
    effect_allele: ALT
    complete_region: false
    export_p_threshold: 1.0e-7
    rationale: Match the measured transcript isoform ratio to its own cis genetic associations for instrument assessment.
```

The current streaming adapter expects columns `SNP`, `SNP_Chr`, `SNP_Pos`, `EAF`,
`N`, `Tx`, `Is_Cis`, `RSq`, `Fx`, `T`, and `log10P`, with SNP IDs such as
`chr1:123:G:A:rs123`. It checks exact transcript versions, parses declared ALT
effects, retains log10 probabilities to avoid underflow, and derives candidate
standard errors as `abs(Fx/T)`. That derivation requires the statistic definition
to be correct. Flags highlight rare alleles, near-perfect fits, extreme statistics,
weak instruments and the usual MR instrument threshold. Flags are review aids,
not proof of artifacts; their absence is not instrument approval. At most 50,000
matching records are written, in source order, with truncation disclosed.

The output is **candidate review, not by itself MR or colocalization inference**. The
candidate-review operation is complete when its audited candidate files exist;
it no longer reports itself as waiting for downstream MR inputs. MR still needs
suitable phenotype GWAS effects, LD, ancestry, overlap assessment and instrument
validity review, supplied either by `analysis.bioinformatics.mr` or custom
`genetic_pairs`. Conventional colocalization needs unfiltered regional statistics
on both sides, including nonsignificant variants. A P < 10⁻⁷ irQTL export cannot
stand in for complete regional evidence. `analysis.bioinformatics.colocalization`
can nevertheless run a deliberately labeled exploratory partial-data analysis and
propagate its mandatory caveat into reporting. Gene-expression QTLs are not silently
substituted for isoform-ratio QTLs. Sharing a consortium or genome build does not
prove independence or matching ancestry.

For a separately reviewed `genetic_pairs` MR analysis, a filtered exposure is
permitted when its declared `export_p_threshold` includes the prespecified
instrument threshold. Outcome statistics must not be selected by their outcome
P-value: declare complete regional data or `selection_basis:
exposure_instrument_lookup`. Harmonization coverage is measured against the
exposure-variant set for MR; regional inference retains its full-region coverage
check. These declarations do not themselves acquire missing data or approve
instruments. Automatic source preparation rejects unsupported non-SNP alleles,
invalid standard errors and duplicate rsIDs row by row, and records exclusion
counts instead of rejecting an otherwise usable feature. If LD clumping leaves
one independent strong instrument, the pipeline reports a Wald ratio with an
explicit limitation: heterogeneity, pleiotropy, Egger, weighted-median and
leave-one-instrument-out diagnostics are not estimable.

## Literature depth and excerpts

The deterministic literature step starts with four field-level strategies:
field background, hypothesis relevance, similar studies, and potentially
conflicting findings. These searches run even when no gene passes the analysis
threshold. They combine declared phenotype labels with molecular-field terms,
inferred from dataset type or supplied explicitly in each hypothesis:

```yaml
analysis:
  bioinformatics:
    phenotype_terms: [body mass index, obesity]
    literature_terms: [RNA splicing, isoform ratio, transcript isoform]
    max_pages: 3
```

`literature_terms` replaces the inferred molecular terms; `phenotype_terms`
replaces the default outcome label. Use synonyms appropriate to the question,
not patient information or a confidential project narrative. Unknown assay
types need explicit field terms. These are bounded title/abstract keyword
searches, not an automated understanding of every aspect of the hypothesis.
Review tissue, species, design and phenotype relevance when screening results.

For eligible genes, three additional strategies address direct gene–phenotype
matches, molecular mechanisms, and replication/contrary findings. The adapter
retrieves Europe PMC core records with available abstracts and records cursor
pagination and search limits. It preserves bibliographic provenance and literal
excerpts of at most twelve words, with offsets in explicitly normalized abstract
text and links to hashed source responses. Excerpts are lexical candidates for
review, not a claim of hypothesis support; unavailable abstracts remain unavailable.
This is not an exhaustive or systematic review. Queries intentionally do not
select only positive results, but vocabulary, indexing and bounded budgets can
still miss important work.

Each literature operation writes `reference-catalog.yaml` and a browsable
`reference-catalog.html`. Publication drafts and manuscript revisions also
consolidate their native literature evidence into these catalogs. Entries retain
available title, authors, year, journal, DOI, PMID/PMCID, source link, literal
excerpt, metadata conflicts and every retrieval context. Context includes the
strategy, query, scientific rationale, page, timestamps and hashed response
or cached-response audit link. The manuscript snapshot additionally identifies
the source operation and manifest; keep the source analysis bundle for its raw
response files.

Records are deduplicated by normalized identifiers, not similar titles. Matching
identifiers connect occurrences without discarding their query histories;
conflicting bibliographic metadata stays visible for review. An article without
a usable identifier remains a candidate but cannot be cited by the revision
command. Catalog IDs are deterministic for the available identifier set, not
permanent identifiers across later metadata corrections. References remain
`unreviewed_candidate_not_hypothesis_support` after retrieval.

Coverage is reported for each strategy and operation/phenotype scope. A completed
empty search differs from missing inputs, denied transmission, offline access,
failed requests and truncated results. Even complete bounded coverage requires
scientific review: assess synonyms and database/date coverage, read the relevant
source text, examine contrary findings and compare the closest prior studies.
Neither a paper count nor an empty search establishes sufficient scholarship or
that a study is the first of its kind. Extend the search through a study amendment
when needed. The reporting-only revision command does not acquire missing papers.

Opt-in LLM review is instructed to distinguish quotations, interpretation and
proposed experiments; give design, species, tissue, trait and contrary evidence;
and justify publication priorities separately from therapeutic opportunities.
Generated advice still requires source verification and scientist review. The
new participant summaries are **not automatically added to remote LLM payloads**.
Literature and draft advisors receive at most 20 catalog candidate previews,
ordered by reference ID rather than relevance. Editorial review receives all
cited catalog entries plus at most 20 additional field/cohort/method-search candidates. These
previews omit full abstracts, local paths and query-audit records; they cannot
replace review of the complete local catalog. New sources suggested or retrieved
by an advisor remain separately audited candidates, not automatically verified
manuscript references. No remote provider is enabled by this writing policy.

The [eQTL Catalogue Data Access page](https://www.ebi.ac.uk/eqtl/Data_access/)
reports that its REST API has been retired (checked 2026-10-02). JLinAlg records
that service state rather than retrying the withdrawn endpoint or treating it as
no biological matches. FTP/tabix acquisition remains a separately reviewed source
route, not an implemented automatic replacement in this adapter.

## Revise a manuscript without rerunning the analysis

The current source build provides a reporting-only, offline revision command:

```powershell
java -jar jlinalg.jar study manuscript revision.yaml --project SOURCE_PROJECT --out NEW_REVISION_BUNDLE --actor-type human --actor-id scientist-1
java -jar jlinalg.jar study manuscript-verify --project NEW_REVISION_BUNDLE
```

Use the rebuilt executable JAR; the version label alone does not establish that
an older release contains this command. The source project remains unchanged.
The command verifies its ledger and outputs, reconstructs current evidence from
authoritative artifacts rather than trusting an editable index, and creates a
new bundle. It does not refit models, amend hypotheses, query databases or call
an LLM. A supplied LLM-authored editorial request should use `--actor-type llm`
and `--model` with the known model identity, or an explicit unknown declaration.

### Develop the Introduction from the field to the hypotheses

Native first drafts, revision checks and LLM drafting/editorial instructions share
the same four- or five-paragraph structure:

1. One or two paragraphs summarizing the relevant field and established findings.
2. One paragraph explaining why the study hypothesis is relevant to that field.
3. One paragraph explaining why the question matters, what closely related studies
   have already done, and precisely how this study differs.
4. One final paragraph stating the overall hypotheses and their scientific motivation.

Use citations for background, rationale and prior-work comparisons. A different
assay, tissue, population or validation design is a specific contribution to
explain, not proof of novelty. Do not invent a prespecified motivation from the
significant results; changed hypotheses belong in the study amendment history.
The finished Introduction uses connected prose without paragraph mini-headings
or editorial/audit commentary.

Generated bundles provide `introduction-template.yaml`, a revision request with
the four roles and explicit placeholders, and `introduction-review.yaml`, which
lists missing roles, citation needs and literature coverage. Fill the template
with reviewed prose and cataloged references. For example, a paragraph under
`sections` can have:

```yaml
- id: intro-prior-work
  section: introduction
  introduction_role: importance_and_gap
  text: >-
    [SCIENTIST INPUT REQUIRED: Compare the closest studies and explain the
    specific contribution, using sources that actually support the comparison.]
  references: []
  sources: []
  rationale: Develop the importance and prior-work comparison after source review.
```

The ordered roles are `field_background`, `hypothesis_relevance`,
`importance_and_gap`, and `hypotheses_and_motivation`. The first permits two
paragraphs; every other role permits one. Blank lines separate paragraphs.
Unknown roles, roles outside the Introduction, excess paragraphs and reversed
role order are rejected before a revision bundle is created. Missing roles,
placeholders and missing citations instead remain explicit writing needs.
Existing untagged prose is preserved and flagged for classification; it is not
silently rewritten. A structurally complete Introduction still needs scientific
review: JLinAlg cannot verify source support merely by checking citation IDs.

### Start Methods with the cohorts and then the performed analyses

The shared Methods model begins with a description and design/recruitment
citation for **each participating cohort**, followed by association analysis,
discovery/replication where applicable, meta-analysis, enrichment/ontology,
Mendelian randomization, network analysis, fine mapping where applicable,
colocalization, sensitivity analyses and reproducibility. Include the methods
actually performed, not every possible pipeline branch. A blocked MR or
colocalization step must not read as a completed analysis.

For established cohorts, select the design publication matching the particular
cohort or subcohort. Examples include the
[Framingham Offspring design paper](https://pubmed.ncbi.nlm.nih.gov/474565/),
[CARDIA design and recruitment paper](https://pubmed.ncbi.nlm.nih.gov/3204420/),
[WHI design paper](https://pubmed.ncbi.nlm.nih.gov/9492970/), and
[JHS recruitment and cohort-description paper](https://pubmed.ncbi.nlm.nih.gov/16317982/).
These are starting points for source selection, not interchangeable references
for every generation, ancillary sample or examination. In particular, a cohort's
original enrollment count must never become the current study's analysis N.
Describe the analyzed visits, recruitment subset, inclusion/exclusion criteria,
tissue and sample selection from this study's actual provenance.

The cohort contract now accepts public names and selected references:

```yaml
cohorts:
  - id: fhs-offspring
    name: Framingham Offspring Study
    literature_name: Framingham Offspring Study
    references: [pmid:474565]
    description: >-
      [SCIENTIST INPUT REQUIRED: describe this study's included examination,
      eligibility criteria and sample selection; do not copy historical N.]
```

`literature_name` is an explicit public search label. Otherwise JLinAlg uses an
explicit `name`, or expands the recognized identifiers FHS, CARDIA, WHI, JHS,
`fhs-offspring` and `fhs-gen3` to public names. Unknown private cohort identifiers
and descriptions are not used as search queries; the index requests a public
name instead. Literature searches now include cohort design/recruitment queries
and exact lookups of declared DOI, `pmid:` or `pmcid:` references. A search result
is still a candidate: JLinAlg does not select the first hit as the definitive
cohort paper. A declared reference is inserted into the writing template only
after its identifier is present in the retrieved catalog.

Method references may also be supplied in `datasets[].analysis_description.references`
for externally supplied associations, or in the relevant analysis configuration,
for example `analysis.meta_analysis.references`. These fields select reporting
references; they do not change the statistical estimator. New native association,
meta-analysis, enrichment and regional genetic runs record `details.method_report`
with method identity, settings or links to their authoritative logs, software and
available references. The literature branch searches these recorded method terms
and reference identifiers as well. Cohort/method query rationale, status and
provenance are retained separately from Introduction search coverage.

For every performed analysis, cite the method and report the relevant settings:

| Method | Reporting requirements |
| --- | --- |
| Association | Model family/link, formula orientation, covariates, technical factors, transformations, missing-data handling, sample selection, pedigree/relatedness adjustment and multiple-testing rule. Distinguish externally supplied results from a model fitted by JLinAlg. |
| Meta-analysis | Effect scale, weighting, fixed/random-effects estimator, heterogeneity, minimum cohorts, inference/confidence settings and overlap handling. |
| Enrichment and ontology | Selected features/genes, threshold, mapped tested-gene universe, mapping rules, database release, test and multiplicity correction; cite the statistical method and resource. |
| MR | Molecular exposure/outcome and tissue, GWAS/QTL sources, instrument P/F thresholds, LD reference and clumping scope/r2, allele harmonization, sample overlap, estimators and sensitivities. |
| Network | Network source/version, node/edge definitions, weights, filtering thresholds, algorithm and null model. The writing model does not add an executable network adapter. |
| Colocalization | Region/build, full regional inputs, sample sizes, LD, single/multiple-signal model, priors, posterior interpretation threshold and sensitivity to priors. |

Explain parameter choices where a justification was recorded. Merely using a
default is not a scientific justification. New meta-analysis records capture
backend defaults at execution; old runs with missing settings remain flagged
for recovery from original logs, rather than acquiring today's defaults. The
regional adapter records actual instrument thresholds and estimators, fine-mapping
settings/logs, and colocalization priors. These reporting additions do not refit
historical analyses or retroactively justify decisions.

Each new manuscript bundle contains `methods-template.yaml`,
`methods-inventory.yaml` and `methods-review.yaml`. Use the template's ordered
`methods_role` values and preserve source IDs. Cohort paragraphs identify
`cohort_ids`; analytical paragraphs carry `parameters` and
`parameter_justifications` mappings. For example, an illustrative paragraph
under `sections` may contain:

```yaml
- id: methods-meta
  section: methods
  methods_role: meta_analysis
  text: >-
    [SCIENTIST INPUT REQUIRED: describe the verified estimator, parameters
    and supporting method citation in publication prose.]
  sources: [H1.meta_analysis]
  references: []
  parameters:
    source: Verified method_report and original command log for H1.meta_analysis
  parameter_justifications:
    estimator: Not recorded; request the scientist's rationale.
  rationale: Document the actual pooled analysis without inventing its rationale.
```

Fill this from the actual run and catalog; an empty citation list remains a
review need. The role order is recorded in the template policy, with `cohorts`
first and `reproducibility` last. Invalid or reversed roles are rejected before
bundle creation. Missing cohort descriptions, citations, parameters or
justifications are review needs, not fabricated text. Existing untagged Methods
paragraphs and their attribution remain intact; classify or explicitly replace
them in a new revision. Structural checks do not establish source support or
certify that prose accurately describes an executed analysis.

LLM drafting and editorial prompts use the same model. Automatic assistance does
not receive raw cohort descriptions, patient files or the local methods inventory.
Explicitly submitted paragraph parameters and justifications are part of the
hash-approved editorial payload; inspect them for sensitive text and paths before
authorizing disclosure. Offline and deny-all policies still block remote queries.

## Results and Discussion writing model

Results begin with native Table 1, followed by association, discovery and
replication, meta-analysis, enrichment and ontology, Mendelian randomization,
network analysis, fine mapping, colocalization, and sensitivity analyses, as
applicable. These are the same analytical roles and order used by Methods; an
unperformed analysis does not become a result. Important null and discordant
findings belong with their analysis, not hidden in the supplement. There is no
five-paragraph Results limit.

Tag each Results paragraph with `results_role` and its primary `sources`.
`supporting_evidence` places bioinformatics context beside the finding it informs.
Each entry names a current evidence `source` and explains the `relationship`,
including the molecular feature, tissue and phenotype match, concordance or
conflict, and interpretive limits. Paragraph `references` cite the corresponding
retrieved literature. These links record an interpretation; they do not certify
independence, mechanistic support or isoform-specific causality. A partly blocked
operation may retain useful retrieved records, but its incomplete status remains
visible and its usable source-level result must be identified and reviewed.
Native drafting inventories candidate context without attaching every database
match to every finding.

Discussion begins with one connected paragraph about the most impactful finding
and why it matters, followed by one paragraph each for the second and third.
Use `discussion_role: principal_finding`, consecutive `impact_rank` values 1–3,
`result_ids` linking the relevant Results paragraphs, and `impact_reason`.
Express that reason in the scientific prose as well as the audit field. Consider
effect magnitude and uncertainty, replication, biological relevance, sensitivity
and conflicting evidence; a small P-value or completed operation is not by itself
a reason for high impact. Different findings can originate in the same analysis,
but repeated Results links trigger a check that the priorities are genuinely
distinct. If fewer than three findings merit discussion, remove unused template
slots and explain why in the top-level `discussion_selection_note`.

After those paragraphs, use `future_potential`, `next_steps`, and `limitations`
in that order. Distinguish possible applications and proposed follow-up from work
actually performed; explain how limitations constrain interpretation. Finish
with `overall_conclusion`, either as the final Discussion paragraph or in the
separate `conclusion` section, not both. It should synthesize the supported
contribution, not introduce a new finding or promise clinical utility.

Every new native report and manuscript revision includes:

- `results-discussion-template.yaml`: editable ordered paragraphs. The three
  impact slots are explicitly unfilled, not automatically ranked analyses.
- `results-discussion-inventory.yaml`: eligible analysis sources, available
  bioinformatics and literature context, and unresolved evidence.
- `results-discussion-review.yaml`: missing structure, reasons, evidence links
  and unclassified older prose requiring review.

For example, these illustrative paragraphs can be included under `sections`
after replacing source IDs and prose with the actual verified findings:

```yaml
- id: results-meta
  section: results
  results_role: meta_analysis
  title: Pooled molecular associations
  text: >-
    [SCIENTIST INPUT REQUIRED: report estimates, uncertainty, heterogeneity,
    and important null or discordant results.]
  sources: [H1.meta_analysis]
  references: []
  supporting_evidence:
    - source: H1.bioinformatics_investigation
      relationship: >-
        [Replace with the exact feature-level relationship and limitations;
        a protein annotation alone is not isoform-specific validation.]
  rationale: Present pooled results in Methods order with relevant context.
- id: discussion-primary
  section: discussion
  discussion_role: principal_finding
  impact_rank: 1
  result_ids: [results-meta]
  impact_reason: >-
    [SCIENTIST INPUT REQUIRED: explain why this is the most impactful finding.]
  text: >-
    [SCIENTIST INPUT REQUIRED: interpret this finding, explain its importance,
    compare relevant literature, and retain uncertainty and contrary evidence.]
  sources: [H1.meta_analysis]
  references: []
  rationale: Distinguish scientific importance from statistical significance.
```

Unknown sources, mismatched Results roles, reversed typed ordering, invalid impact
ranks, unknown Results links and duplicated closing roles are rejected before a
revision bundle is created. Missing prose, rationale, classification or evidence
coverage remains a review need. Untagged older paragraphs retain their order,
text and attribution; JLinAlg does not silently rewrite an earlier manuscript.
Structural completeness never certifies scientific validity. LLM drafting and
editorial prompts follow this same model; suggestions still require acceptance
in a new audited revision. Explicit relationship text and impact reasons enter
the hash-approved editorial payload, but the full local evidence inventory and
patient files are not automatically transmitted.

## Creating a manuscript revision

A small revision request can contain:

```yaml
schema_version: 1
id: manuscript-01
title: Molecular associations with the study phenotype
rationale: Integrate the latest verified findings and retain uncertainty.
sections:
  - id: main-association
    section: results
    results_role: meta_analysis
    title: Pooled molecular associations
    text: >-
      The analysis tested {{H1.meta_analysis#/tested_features}} features;
      {{H1.meta_analysis#/BH_FDR_0_05_features}} met the displayed false
      discovery rate threshold. These associations do not establish causality.
    sources: [H1.meta_analysis]
    rationale: Put the completed primary analysis before exploratory context.
comments:
  - id: check-population
    comment: Confirm that Table 1 represents the analysis population.
    status: pending
```

Replace `H1` with the actual hypothesis ID. Bindings use a hypothesis/operation
key plus a JSON pointer into its native `key_findings`. Missing, ambiguous,
malformed or nonscalar bindings are errors, not guessed values. The claim map
records the exact substituted values and their source manifests. Ordinary prose
is attributed editorial input and is **not** automatically fact-checked.

Allowed sections are `introduction`, `methods`, `results`, `discussion` and
`conclusion`. Results start with native participant characteristics, retaining
unverified population and missing-coding warnings. Typed Results paragraphs must
follow Methods order, each requiring matching completed analysis evidence and an
editorial rationale. Untagged older paragraphs remain intact and flagged for
classification. Important caveats remain in the manuscript, while every
current evidence layer is retained in the supplement. A completed operation is
not necessarily positive evidence. Optional `references` lists DOIs (also accepted
as DOI URLs), `pmid:123`, `pmcid:PMC123`, or `ref-...` catalog IDs already present
in the current native literature results; a missing reference is rejected.
Retrieval establishes provenance, not scientific support.

For another revision, add a `previous_revision` mapping with `manifest` pointing
to the earlier bundle's `metadata/manifest.yaml` and its `sha256`. Supply only
changed paragraph IDs: omitted paragraphs and their author attribution are
retained, while replaced paragraphs have before/after records. To remove a
paragraph, use `remove_sections: [{id: paragraph-id, rationale: reason}]`.
To move existing prose or place a new paragraph before an older one, optionally
supply `paragraph_order`, listing **every retained paragraph ID exactly once**.
The change log records that explicit ordering decision without changing the
authorship of untouched paragraphs. Omit it on later revisions to retain the
current order; an old ordering request is not replayed on new additions.
Comments merge by ID and have `pending`, `addressed` or `deferred` status;
addressed comments require a `response`. Numerical bindings refresh against the
current evidence, so reread all retained interpretations after scientific changes.

The bundle includes `artifacts/manuscript.html`, `supplement.html`, a deduplicated
`literature-review.html`, `reference-catalog.html`, `reference-catalog.yaml`,
`introduction-template.yaml`, `introduction-review.yaml`, `evidence.yaml`,
`methods-template.yaml`, `methods-inventory.yaml`, `methods-review.yaml`,
`results-discussion-template.yaml`, `results-discussion-inventory.yaml`,
`results-discussion-review.yaml`,
`publication-claims.yaml`, and a local HTML index. Metadata preserves the request, resolved revision, previous drafts
as inert text, source manifests, original actor attribution and editorial change
log. The running packaged JAR is copied to `scripts/`; a class-directory test run
explicitly reports that no packaged JAR is available. A notebook, timestamps,
software identity and checksum manifest document generation. Keep the original
analysis bundle: large result tables and public response caches are referenced,
not duplicated into this reporting snapshot.

Never edit a checksummed generated draft in place. Free-form edits in an external
editor are not imported automatically; put revised text into the next request's
section IDs. Scientific changes still require study amendments and analysis
execution. This reporting command cannot establish missing sample membership,
sex coding, model covariates, valid instruments or regional colocalization input.
No patient-level rows are read by the revision command, but editorial text and
aggregate results may remain sensitive; local generation does not authorize
sharing. Integrity verification is not publication or regulatory certification.

### Request an optional LLM editorial pass

Each revision prepares `artifacts/editorial-payload.json` locally, without sending
it. This contains the manuscript sections, editorial comments, bounded aggregate
association/replication/enrichment/irQTL previews, and cited literature metadata.
It excludes input files, patient rows, automatic Table 1 summaries, local paths,
raw response caches and free-text scientist answers. Manuscript prose or comments
can themselves contain sensitive information: inspect the **whole exact payload**.
An allowlist is not a guarantee of anonymization or protection against disclosure.

After that inspection, authorize the payload by its exact SHA-256:

```powershell
$payload = 'REVISION_BUNDLE/artifacts/editorial-payload.json'
$approved = (Get-FileHash -Algorithm SHA256 -LiteralPath $payload).Hash.ToLowerInvariant()
# Review $payload before using its hash as disclosure approval.
java -jar jlinalg.jar study manuscript-editorial --project REVISION_BUNDLE --out NEW_REVIEW_BUNDLE --approved-payload-sha256 $approved --actor-id scientist-1
java -jar jlinalg.jar study manuscript-verify --project NEW_REVIEW_BUNDLE
```

Omitting `--provider` selects **browser handoff**, not a remote request. No model,
key, account, browser installation or internet connection is needed to prepare
the packet. Open `NEW_REVIEW_BUNDLE/artifacts/handoff/index.html` and review the
exact `prompt.txt` before manually submitting anything. `--offline` still allows
local preparation and import; it does not authorize browser transmission.

### Complete the browser handoff

The local packet is plain HTML and UTF-8 text, with no JavaScript, external assets,
forms, automatic clipboard access or remote requests. Keep it unchanged. To open
it explicitly:

```sh
java -jar jlinalg.jar study handoff-open --handoff NEW_REVIEW_BUNDLE/artifacts/handoff
# Alternatively select --browser chrome, chromium, firefox, edge, links or none.
```

`auto` first requests the system-default browser; if unavailable, it looks for
installed Chrome, Chromium, Firefox or Edge in a graphical environment, then
Links in an interactive terminal. An explicit selection does not silently choose
another family. It opens only the verified local packet, reports the selection
or manual-opening path, and never installs a browser. Launch success does not
prove that a page rendered. A headless server can prepare the packet for an
authorized workstation instead.

Browser choices were evaluated against their documented capabilities:

| Option | Role in this workflow |
| --- | --- |
| Chrome, Chromium, Firefox or Edge | Recommended for the actual modern AI website; use the installed browser that works with the scientist's account. No website compatibility or login is assumed. |
| Linux Links | Lightweight local packet viewer and terminal fallback. Its [feature page](https://links.twibright.com/features.php) describes HTML 4 without CSS, and its [manual](https://links.twibright.com/user_en.html) says JavaScript is deactivated. Do not depend on it for AI Mode/Gemini/ChatGPT web applications. |
| Browsh | A possible terminal browser for modern pages, but its [architecture requires Firefox and an extension](https://www.brow.sh/docs/introduction/). Not bundled or required; AI-site compatibility has not been tested. |
| Automated browser drivers | Not part of this handoff implementation. Website changes, authentication, disclosure controls and provider terms need separate evaluation. No CAPTCHA bypass, private endpoint replay or cookie extraction. |

If policy permits, manually copy **only the exact reviewed prompt** into the
chosen AI website. Do not upload the study folder, phenotype files, full evidence
index or other attachments. The approved payload excludes patient files, but
free-text manuscript content and unpublished aggregate findings can still be
confidential. Provider terms matter even for free services.

Copy `handoff/response-template.yaml` **outside the frozen packet**. Fill in its
provider, displayed model (`not_disclosed` if unknown), clean HTTPS website URL,
quoted ISO-8601 submission/receipt timestamps, complete response, references and
review rationale. Reference entries have `title` and `url`. Record only URLs
without query strings, fragments or credentials; never include login tokens.
Confirm `exact_prompt_submitted`, `disclosure_approved` and `no_patient_data`
only when true. If the prompt was edited, prepare a new packet; do not attest an
exact match. `decision` is `pending_review`, `accepted_for_revision` or `rejected`.

```sh
java -jar jlinalg.jar study handoff-import response.yaml --handoff NEW_REVIEW_BUNDLE/artifacts/handoff --out NEW_RECEIPT --actor-id scientist-1
java -jar jlinalg.jar study handoff-verify NEW_RECEIPT/receipt.yaml
```

The new receipt preserves the prompt, response, request identity, checksums,
provider/model declarations, timestamps, named human actor, disposition and
`SESSION_HISTORY.md`. Response HTML and code are displayed as escaped text,
never executed. Existing packets/receipts cannot be overwritten. Importing an
accepted answer still does not edit a manuscript or execute a new analysis.

**Audit limit:** JLinAlg verifies local integrity and request matching. It did
not observe the browser interaction, authenticate the provider/model, establish
that all conversation turns were supplied, or independently verify citations.
External submission times, no-patient-data assertions and response attribution
are human attestations, not certified telemetry. Store only one exact-prompt
exchange per receipt; create another packet for follow-up context rather than
silently omitting the conversation. The archive explains this gap explicitly.

### Resume pipeline assistance

For a pipeline, enable the desired purposes and optionally predeclare receipt
paths before the first run. Omitted provider means `browser`; enabling assistance
is still explicit and does not itself authorize disclosure:

```yaml
assistance:
  enabled: true
  next_steps: true
  literature: true
  response_files:
    - responses/H1-advisor/receipt.yaml
    - responses/H1-literature/receipt.yaml
```

The project index links each waiting operation's `handoff/index.html` and states
the required input. Import each response into its declared **new** receipt
directory, then repeat the same `study auto` command. Missing receipt paths are
allowed initially. Matching uses the exact prompt, purpose, hypothesis and policy
identity, not filenames alone. Completed unchanged upstream operations are reused;
changed response artifacts invalidate assistance caching. Receipts are copied
into the operation's immutable outputs. Declare more paths by an approved
amendment if they were not planned initially. Do not edit the base study in place.

An unrelated receipt is not applied, duplicate matches are rejected, a rejected
response stays waiting, and malformed suggestions cannot authorize execution.
For next-step advice the imported answer must satisfy the requested JSON
`proposals` contract. Literature references are labelled human-transcribed and
unverified, not provider-retrieved evidence. Empty citations are not proof of no
relevant literature. Editorial receipt suggestions are incorporated by creating
a new `study manuscript` revision and citing the receipt in its rationale.

### Optional API key and OAuth access token connections

Use an API connection only when explicitly selected. `--provider openai` uses
the OpenAI Responses API; `--provider gemini` uses the Gemini API. Each needs its
model and environment-variable credential. The portable alternative is
`--provider openai_compatible`, which uses non-streaming Chat Completions:

```sh
java -jar jlinalg.jar study manuscript-editorial --project REVISION_BUNDLE --out NEW_API_REVIEW --approved-payload-sha256 HASH --provider openai_compatible --model PROVIDER_MODEL --endpoint https://gateway.example.org/v1/chat/completions --auth api_key --api-key-env RESEARCH_AI_KEY --actor-id scientist-1

# Only when the endpoint officially accepts a token issued by its OAuth service:
java -jar jlinalg.jar study manuscript-editorial --project REVISION_BUNDLE --out NEW_OAUTH_REVIEW --approved-payload-sha256 HASH --provider openai_compatible --model PROVIDER_MODEL --endpoint https://gateway.example.org/v1/chat/completions --auth oauth --oauth-token-env RESEARCH_AI_ACCESS_TOKEN --actor-id scientist-1
```

The corresponding study keys are `provider: openai_compatible`, `model`,
`endpoint`, `auth: api_key` plus `api_key_env`, or `auth: oauth` plus
`oauth_token_env`. Use one authentication method. Endpoint URLs must use HTTPS
and end in `/chat/completions`, without credentials, query parameters or
fragments. Redirects are disabled. Tokens are read from the environment and
sent only as bearer credentials to that explicit endpoint; echoed tokens are
redacted from stored responses. Never place credentials in YAML, CLI arguments,
manuscripts or browser transcripts. The [OpenAI authentication reference](https://developers.openai.com/api/reference/overview)
documents bearer authentication; individual compatible gateways may differ.

OAuth here means **using an already issued API access token**. JLinAlg does not
implement provider registration, browser/device login, consent or token refresh,
and does not claim that a ChatGPT/Google consumer subscription grants API access.
Use the provider's supported authentication tooling. An expired/denied token is
recorded, with no silent API-key/provider fallback. Generic Chat Completions does
not guarantee web search or retrieval citations; the literature branch remains
waiting if the provider returns no usable retrieval evidence. Model capabilities,
account access, terms, retention and costs must be checked for that endpoint.

The command refuses a mismatched or unanchored payload. Source `deny_all` policy
still wins over payload approval. Offline restrictions block API calls, not local
browser packet preparation/import. Additional project
disclosure constraints unsupported by the adapter remain blocking; approval does
not silently relax the study policy. The reviewed payload and policy are preserved
with the source-revision hash and approving actor. Exact prompts, credential-redacted
responses, requested/returned model, timestamps, request disposition and editorial
decision are retained. Missing credentials, denials, timeouts and refusals appear
in `artifacts/index.html`; an archived attempt is not necessarily successful review.

The model proposes section-ID replacements, reasons, source checks and scientist
questions. It cannot change the draft, numerical evidence or analytical plan.
Accept selected suggestions by entering those sections in a **new** manuscript
revision and recording the editorial-review bundle in the rationale. Previous text
and attribution remain auditable. No generated code or HTML is executed, no new
web search is requested by this editorial pass, and no fallback or automatic retry
occurs. A new attempt can incur another provider charge.

The OpenAI adapter uses `store: false` as described in the
[Responses API documentation](https://developers.openai.com/api/docs/guides/migrate-to-responses).
That setting does not by itself establish zero retention; review the provider's
[data controls](https://developers.openai.com/api/docs/guides/your-data) and your
institution's disclosure requirements. Leaving the editorial pass unused keeps
manuscript revision entirely local.

For an already installed **local** model, the explicit editorial-only adapter
also supports a keyless loopback chat-completions endpoint. This is not the
default and does not require or install any particular model application:

```powershell
java -jar jlinalg.jar study manuscript-editorial --project REVISION_BUNDLE --out NEW_LOCAL_REVIEW --provider local --model YOUR_LOADED_MODEL --endpoint http://127.0.0.1:1234/v1/chat/completions --approved-payload-sha256 $approved --actor-id scientist-1
```

Only literal IPv4 loopback, an explicit port and the exact chat-completions path
are accepted. DNS names, remote addresses, credentials in URLs, query strings,
proxies and redirects are rejected or disabled. No API key is read or sent. A
remote-transmission `deny_all` policy does not prohibit this loopback-only local
inference; `--offline` blocks even local HTTP and unsupported extra disclosure
constraints still block the request. JLinAlg does not install, start or download
models, and a loopback address cannot prove that a user-managed server does not
itself forward data elsewhere. Use a trusted local-only model configuration.
The recorded model name is the server's self-report, not a verified weights hash.
Model identity/weights and server version should be included in deployment
provenance. A missing or stopped server is recorded as unavailable, not successful
editorial review. This local adapter is currently editorial-only.

### Local editorial review with llama.cpp

`--provider llama_cpp` explicitly selects the same keyless, loopback-only transport
and records the declared runtime as llama.cpp. No LM Studio dependency is needed.
Use an already installed, trusted `llama-server` and an instruction-tuned local
GGUF model. The [upstream server documentation](https://github.com/ggml-org/llama.cpp/tree/master/tools/server)
describes its Chat Completions endpoint and model aliases. For example:

```sh
llama-server -m /path/to/review-model.gguf --host 127.0.0.1 --port 8080 --alias editorial-local --ctx-size 32768
java -jar jlinalg.jar study manuscript-editorial --project REVISION_BUNDLE --out NEW_LLAMA_REVIEW --provider llama_cpp --model editorial-local --endpoint http://127.0.0.1:8080/v1/chat/completions --approved-payload-sha256 HASH --actor-id scientist-1
```

Choose a context size the model and hardware support, large enough for the full
prompt plus output; 32768 is an example, not a memory-capacity guarantee. JLinAlg
requests at most 4000 output tokens and rejects truncated completions rather than
presenting them as successful reviews. The alias avoids leaking model filesystem
paths in request metadata. Record the GGUF checksum, server build/version and
launch arguments separately; the alias is not proof of model identity.

Use local weights, not a cloud-forwarding server or remote RPC backend. JLinAlg
does not download or launch llama.cpp or its models. `deny_all` permits this local
route; `--offline` blocks local HTTP too. Do not expose the unauthenticated server
on a public/LAN interface. An actual local inference run requires a working server;
adapter regression tests alone do not establish model quality or hardware fitness.

### ChatGPT sign-in and an external Codex editorial client

A Codex login is not an API key to copy into JLinAlg. The supported
[Codex non-interactive client](https://learn.chatgpt.com/docs/non-interactive-mode)
can use its own saved sign-in for a reviewed manuscript-only request. Keep
credentials inside that client; record its version, exact prompt, response,
model identity when available, timestamps and decisions. Disable unrelated tools
and do not run the editor in a patient-data directory. Feed only the approved
aggregate payload. This external-client workflow is distinct from JLinAlg's
manual browser receipts: never label an automated exchange as human-transcribed.

OpenAI also documents a separate
[Sign in with ChatGPT flow for eligible open-source apps](https://developers.openai.com/siwc/token-sharing-open-source).
That requires app registration, user consent, authorized scopes and token lifecycle
handling. It is not implemented by JLinAlg's generic `--auth oauth` option, and
existing Codex session credentials are not automatically a JLinAlg authorization.
The documented ChatGPT-plan route uses streaming Responses requests; do not point
the non-streaming Chat Completions adapter at it or reuse an unrelated session token.
