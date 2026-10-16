# Auditable study pipeline tutorial

For follow-up beyond replication and meta-analysis, see the [bioinformatics evidence tutorial](bioinformatics-evidence.md): source discovery, reviewed regional MR/fine-mapping/colocalization, proteins, genomic sequence, literature and local evidence archives. It includes a separate runnable offline synthetic example.

> **Source-build tutorial.** This walkthrough runs the checked-in, result-first
> two-cohort example through validation, resolution, plan approval, execution,
> verification, and resumable reuse. It also shows where an optional LLM may assist
> and how to keep patient-level data outside a remote LLM boundary.

## What you will build

The example begins with two association-result tables rather than participant
matrices. JLinAlg will standardize the inputs, compare discovery and replication
effects, perform fixed-effect meta-analysis, run deterministic diagnostics, and
write a traceable report. Every committed step is represented in the event
ledger and the readable research notebook.

The tutorial is intentionally small enough to inspect by hand. The same contract
and audit structure apply when a project begins with processed molecular
matrices, has more cohorts or tissues, or stops at a later evidence stage.

The checked-in fixture uses `analysis.mode: manual` to stay small and completely
offline. In a production study, omit that field or set `mode: auto`. Explicit
goals then become the minimum requested work, and JLinAlg adds eligible
association, replication, meta-analysis, enrichment, genetic-evidence,
literature, diagnostic, prioritization, and report branches from the objective
and declared inputs.

The default pipeline does not contact an LLM. An explicitly enabled `assistance`
block defaults to an audited manual browser handoff for bounded advice,
literature review and a first manuscript draft. Explicit alternatives include
OpenAI Responses, Gemini GenerateContent and OpenAI-compatible Chat Completions
with an API key or provider-issued OAuth bearer token. Browser transcripts and
citations are human-attested, not authenticated provider output. JLinAlg retains
control. MCP, provider-specific OAuth login/refresh, general pipeline local-model
adapters and autonomous LLM orchestration remain future work. The `llm` actor
option alone never enables a call.

## Effect-size quality control before downstream analysis

Numerically finite coefficients can still be extreme: a nearly constant predictor,
ill-conditioned design, separation, unit mismatch or unstable fit can produce a
very large coefficient and SE. A large raw coefficient is not itself proof of a
bad association, and summary statistics cannot establish the cause. Do not hide
these points by clipping a plot or choose exclusions to improve correlation or
significance. Inspect predictor variation, units and upstream fit diagnostics.

JLinAlg now applies an auditable, conservative effect-size QC screen to supplied
association results during harmonization and to newly fitted association results
before annotation. The defaults are:

~~~yaml
analysis:
  harmonization:
    effect_outliers:
      enabled: true
      robust_z: 100
      min_features: 20
      # max_abs_effect: 1000  # optional; choose only with known effect units
~~~

For each dataset separately, using all numerically eligible association rows,
the center is the median beta and the robust scale is
`1.482602218505602 * median(abs(beta - median(beta)))`. If this scale is zero
or nonfinite, use `(Q75 - Q25) / 1.348979500392164` with linearly interpolated
quantiles. Exclude estimates with `abs(beta - center) / scale > robust_z`.
This is a configurable QC heuristic, **not** a calibrated hypothesis test or a
universal biological plausibility cutoff. The conservative default of 100 is
intended to catch gross extremes, not trim an ordinary heavy tail. Uniform
changes of effect units do not change the relative selection. It uses neither
P values nor agreement with other cohorts. Relative screening is explicitly
skipped below `min_features` or without a positive finite robust scale; no
epsilon denominator is invented. An optional `abs(beta) > max_abs_effect`
rule still applies to these families. Set `enabled: false` for an explicitly
unfiltered sensitivity run.

The operation writes `effect-outlier-summary.tsv`,
`effect-outlier-exclusions.tsv` (dataset, feature, beta, SE, P, robust score and
reason) and `effect-outlier-qc.yaml` (effective settings, counts, fallbacks,
skips and caveats), all included in its output manifest and Methods inventory.
The source files are untouched. Excluded estimates are removed from the generated
standardized tables, so comparisons/correlations, discovery/replication, FDR,
meta-analysis, sensitivity analyses, enrichment, MR and colocalization selection
consume the same QC-eligible results. With `common_universe: intersection`,
the intersection is computed **after** QC: a feature excluded in one supplied
cohort is absent from every supplied cohort's downstream table. Otherwise only
the flagged dataset estimate is removed; valid estimates in other datasets
remain eligible under the usual minimum-cohort rules. The testing family is
the retained eligible family, not the old unfiltered set. If separate matrix
association operations generate new tables, QC is per operation/dataset;
harmonization does not retroactively impose a shared intersection on those
later outputs.

Because selection uses estimated effects, it is not an independent-filtering
procedure that guarantees unchanged FDR calibration merely by avoiding P values.
Prespecify the policy where possible and retain an unfiltered sensitivity analysis.

Record a post hoc QC change as an amendment, review exclusion counts and
unit/model diagnostics, and compare a separately approved unfiltered sensitivity
run before changing scientific conclusions. This software update does not
invalidate caches, overwrite sealed results, redraw historical figures, or repair
an already completed manuscript: a new approved analysis is necessary. Never
claim old results used this screen based on the current software defaults.

## Probe annotation and genomic locations

The source build assesses annotation during result harmonization (or after an
association scan). It writes `feature-annotation.tsv` and
`annotation-summary.yaml`, both covered by the operation manifest. A feature can
be `located`, `missing_location`, `unmapped`, `ambiguous`, or
`incomplete_across_datasets`. Unresolved annotation is reported while association
and replication calculations continue.

Use existing dataset fields plus two optional declarations:

~~~yaml
datasets:
  - id: blood-results
    type: transcriptomics
    stage: association_results
    feature_id_type: probe
    genome_build: GRCh38
    coordinate_system: 1-based-closed
    annotation_release: laboratory-probes-2026-09
    files:
      results: blood-associations.csv
      annotation: probes.tsv
    columns:
      feature_id: Isoform_ID
      effect: Fx_BMI
      standard_error: SE_BMI
      p_value: P_BMI
      annotation_feature_id: Probe
      annotation_gene_id: EnsemblGene
      annotation_gene_symbol: Symbol
      annotation_chromosome: Chromosome
      annotation_start: Start
      annotation_end: End
~~~

Column mappings are unnecessary when the annotation table already has
`feature_id`, `gene_id`, `gene_symbol`, `chromosome`, `start`, and `end`.
`position` may replace start/end for a point annotation. Gene ID/symbol columns
are optional for coordinate-only annotations. `genome_build` may instead be a
table column, but conflicting table and dataset declarations fail. Table
coordinates require `coordinate_system`; BED intervals must declare
`0-based-half-open` and have both start and end. Output intervals are normalized
to `1-based-closed`. An omitted convention cannot establish a usable location.

`files.annotation` also accepts `.gtf` and `.gtf.gz`. GTF gene and transcript rows
provide their defined 1-based coordinates; declare `genome_build` separately.
An explicitly versioned transcript ID must match exactly. Unversioned IDs can
match the declared annotation, with multiple mappings retained as ambiguous.
The [GENCODE GTF format](https://www.gencodegenes.org/pages/data_format.html)
defines these fields. A transcript span is useful for locus discovery; it is
not necessarily the exact splicing junction or intron quantified in the study.
An sQTL follow-up must match the actual measured molecular feature.

For human Ensembl gene/transcript IDs, a declared modern GRCh38 GENCODE release
can trigger the existing public metadata download to resolve gene names. This
metadata provider does **not** supply genomic intervals: provide a local GTF or
annotation table for those. Lifted GRCh37 and nonhuman resources require explicit
local annotations. Resource/source hashes and declarations are recorded.
JLinAlg does not infer genome build from a gene symbol, CpG probe, or rsID.
GWAS result tables with explicit chromosome/position columns can supply their
own coordinates once build and coordinate convention are declared.

Set `feature_id_type: gene_symbol` only when the study IDs really are gene
symbols. Arbitrary probe IDs are never silently interpreted as gene symbols.
The small tutorial uses fictional gene symbols and therefore still needs real
annotation before any genomic follow-up.

## Reading findings and choosing the next analysis

`index/index.html` contains bounded previews of the primary result tables,
discovery/replication comparisons, meta-analysis, and available enrichment.
Feature rows include mapped gene IDs/symbols and mapping status. Full tables and
their SHA-256 hashes are linked. The displayed BH-FDR 0.05 counts are descriptive;
the approved plan retains the actual hypothesis selection/testing rules.
Each retained run has its own result ID; an older run is not silently replaced.

Waiting operations display their required inputs. For non-variant omics results,
variant scoring is `not_applicable` and contributes to a separate skipped count.
Auto-mode source discovery records phenotype-matched GWAS studies and nearby
catalogue SNPs, QTL service availability and tissue candidates. Reviewed canonical
regional inputs can drive automatic genetic inference; arbitrary full-genome
source selection/conversion remains separate. See the [bioinformatics workflow](bioinformatics-evidence.md).
Manual plans without source discovery still report `catalog_search_status: not_run`,
not a negative search. Explicit `analysis.operations` adapters remain available.

Diagnostics write `next-steps.yaml` with each proposal's trigger, required
inputs, and source artifact. Current rules propose scale/model audits for
effect-scale mismatches, uncertainty and power review for near-zero replication,
and fixed/random-effects review for widespread heterogeneity. With at least
three explicitly identified cohorts and a completed primary meta-analysis,
leave-one-cohort-out replication and meta-analysis run automatically as
registered influence diagnostics. Confirm comparable estimands before pooling; fitting
a between-cohort slope is not a justified way to correct measurement scales.
These are proposals, not completed analyses. Record the selected change through
an `analysis_requests` amendment, specify the operation or execute externally
with provenance, and approve the resulting plan. Subgroup analyses additionally
need enough cohorts and relevant metadata; result-only inputs cannot refit
participant-level covariate or batch-adjustment models.

## Advisory AI while JLinAlg controls execution

An optional advisor can suggest related phenotypes or sensitivity analyses
without receiving execution authority. JLinAlg constructs a bounded request,
asks an official API for structured advice, and checks proposed operations against registered methods,
available inputs, and the approved plan. Gemini supports JSON-schema-constrained
[structured outputs](https://ai.google.dev/gemini-api/docs/structured-output);
schema conformance alone does not validate scientific claims.

A general-topic request might contain only a public phenotype term, modality,
and tissue, with no cohort names, gene hit list, statistics, file paths, study
description, or patient records. More detailed advice may require a separately
permitted aggregate-results request. Even a topic-only request discloses the
research topic; keeping execution in JLinAlg does not remove that transmission.
`deny_all` must prohibit all remote advisor requests. A local model can be used
under a no-network policy. Provider terms depend on the selected service and
configuration; assess the [Gemini API terms](https://ai.google.dev/gemini-api/terms)
when implementing the connector.

Advice should include proposed action, rationale, evidence sources, required
inputs, expected benefit, and limitations. Exact prompts/responses, model ID,
timestamps, disclosure manifests, and accepted/rejected decisions enter the
audit trail. JLinAlg may execute already-approved, eligible operations; new
phenotypes, hypotheses, or changes to the estimand require an amendment. An
empty valid response, an unavailable provider, denied transmission, and a search
that was never run have different recorded dispositions. When public catalogs
are downloaded as whole reference snapshots and queried locally, private locus
and phenotype selections do not need to be sent as query parameters.

The implemented adapter sends allowlisted aggregate association/meta-analysis
previews, phenotype labels and up to ten unambiguously annotated top genes.
It does not yet implement a separate topic-only mode or validate phenotype
ontology IDs automatically. Retrieval citations are retained, but their scientific
claims still require review. See the configuration and privacy boundaries below.

## Data dictionaries, optional assistance, and first drafts

An optional local YAML or JSON dictionary supplies the same `columns` mappings
as the study contract. For example, `dictionary-bmi.yaml` may contain:

~~~yaml
columns:
  feature_id: Isoform_ID
  effect: Fx_BMI
  standard_error: SE_BMI
  p_value: P_BMI
descriptions:
  Fx_BMI: Association coefficient for BMI, on the declared analysis scale.
  SE_BMI: Standard error of that coefficient, not its variance.
~~~

Reference it using `datasets[].files.data_dictionary: dictionary-bmi.yaml`.
Explicit `datasets[].columns` values take precedence and each override is
recorded. Resolution snapshots the dictionary, mappings and SHA-256. Editing a
dictionary after plan approval stops execution: save a versioned file and amend
the dataset reference instead. Free-text descriptions are local context, never
instructions to choose a model, change an estimand, or transmit patient data.
CSV/TSV dictionaries are not currently supported. For a file containing BMI and
obesity results, use one dictionary/view per phenotype and explicit mappings.

Optional assistance is disabled unless enabled in the approved study:

~~~yaml
project:
  data_policy:
    llm:
      remote_transmission: aggregate_results  # or deny_all
      remote_patient_level_data: deny
assistance:
  enabled: true
  provider: browser          # default if omitted; no key, account or network to prepare
  next_steps: true
  literature: true
  manuscript: false           # explicit true permits an unreviewed HTML draft
  response_files:             # predeclare one new receipt directory per exchange
    - responses/H1-advisor/receipt.yaml
    - responses/H1-literature/receipt.yaml
  allow_analyses:
    - fixed_random_comparison
    - leave_one_cohort_out
~~~

Browser handoff prepares a local HTML/text packet. The scientist reviews the exact
prompt, submits it manually in an installed modern browser, then imports the
response and repeats the run command. Linux Links can view the packet, but is not
a dependable client for modern AI websites. See the [browser handoff and optional
API authentication tutorial](bioinformatics-evidence.html#complete-the-browser-handoff)
for commands, required attestations, immutable receipts and audit limitations.
No content is submitted by the handoff adapter. Assistance remains opt-in.

For an explicitly selected API provider, supply credentials securely through the
environment. Each requested advisory operation makes at most one provider request
per attempt; calls/search may incur
provider charges. Repeating a timed-out request can incur another charge if the
provider processed it before the local timeout. There is no browser scraping,
login automation, provider fallback, or execution of model-generated code.

The registered fixed/random comparison runs locally when explicitly listed in
`allow_analyses`, even with `enabled: false` or an unavailable LLM. Leave-one-
cohort-out influence analysis is automatic for a hypothesis with at least three
explicitly identified cohorts and a completed primary meta-analysis. It removes
each cohort in turn, reruns the approved meta-analysis and, where discovery and
replication roles remain, reruns replication and its plots. The resulting
`sensitivity-comparison.tsv` records the full set of omissions, replication-rate
and mean-I-squared deltas. A cohort is an `odd_one_candidate` only when omission
improves both available criteria. That candidate is featured in a labeled
sensitivity subsection, never silently substituted for the all-cohort primary
analysis. Promotion to primary requires a prespecified rationale or a reviewed
study amendment. Missing prerequisites appear in the index. These remain
exploratory influence diagnostics; comparable effects and valid cohort-dependence
assumptions are still required.

Other AI suggestions receive `requires_scientist_amendment` or
`unsupported_requires_scientist_review`. New phenotypes, hypotheses and changes
of estimand are never silently approved. Configure an existing project through
an amendment, not by rewriting its base `study.yaml`:

~~~yaml
schema_version: 1
amendment:
  id: advisor-01
  base_study_hash: CURRENT_RESOLVED_STUDY_HASH
  rationale: Permit a bounded review of aggregate findings.
changes:
  configure:
    assistance:
      enabled: true
      provider: browser
      next_steps: true
~~~

Save as `amendment-advisor-01.yaml`, use `study amend ... --check`, then approve
with `--approved-by PERSON`. `changes.configure` supports `assistance`,
`resources`, and project metadata/privacy; normal dataset/hypothesis changes
continue to use `changes.update`.

Every attempted advisory operation records `llm/interaction-001/metadata.yaml`
and its disposition. When applicable, it also stores the exact `query.txt`,
credential-free `request.json`, redacted `response.json`, `answer.txt`, and
`decision.yaml`. Metadata includes UTC start/end times, requested and returned
model, provider, response ID, endpoint, HTTP status and request/response hashes.
API keys are headers, never URL parameters or audit fields. The operation
records an attempted request separately from confirmed HTTP receipt: DNS/TLS
failures and timeouts leave `request_sent: unconfirmed`, not proof of delivery
to a model. A proxy HTTP response also does not establish downstream delivery.
The operation
manifest hashes these files; `SESSION_HISTORY.md` records attempted remote calls
with the LLM actor/provider/model. Reused answers are attributed to their original
attempt, not misrepresented as a new call.

Literature review uses [OpenAI web search](https://developers.openai.com/api/docs/guides/tools-web-search)
or [Gemini retrieval metadata](https://ai.google.dev/api/generate-content).
Only citations returned in provider retrieval metadata enter the reference list;
references invented in model prose are not accepted as retrieved evidence.
No citations, a failed request, and a request not sent are distinct states.
This is a bounded exploratory search, not a systematic review or independent
verification that a cited paper supports each claim.

With `manuscript: true`, `manuscript.html` contains a clearly labeled unreviewed
first draft, supplied aggregate findings and links to `source-map.yaml`, the
operation manifest and model record. The prompt requests Background, Methods,
Results, Discussion and Conclusion, and missing facts must remain placeholders.
The draft may be incomplete or wrong: every number, method, citation and causal
claim needs scientist review. No submission, authorship or ethics approval is
inferred. Instructions and returned HTML-like text are escaped, not executed.

## Publication citations and parameter reporting

### Explaining observed results from the study contract

JLinAlg records `study-context-review.yaml` for advisor requests and manuscript
revisions. It compares declared tissues, compartments, assays, preparation,
population metadata, upstream models, covariates, technical factors, scales,
annotation versions and cohort relationships. Missing values are unknown, not
evidence that cohorts used identical settings. The inventory lists omitted
fields for local review rather than silently assuming that all context is covered.

With `assistance.enabled: true` and `assistance.next_steps: true`, the existing
advisor step asks the configured LLM to consider this context against **all
observed patterns**, including heterogeneity, discordant replication, null or
unexpected findings. The advisor must identify the result and metadata fields,
separate facts from explanations, identify competing explanations/confounding,
and propose a sensitivity check or scientist question. Editorial review receives
the same context under its exact-payload approval. For example, a PBMC/whole-blood
contrast can motivate tissue-sensitive follow-up; it cannot establish a tissue
effect independently of cohort if only one cohort used PBMCs.

The structured projection, exact prompt, response, model, timestamps and decisions
are audited. Free-text notes, arbitrary dictionaries, file contents/paths,
donor/specimen identifiers and credentials are not automatically sent. Relevant
omitted context requires a scientist-reviewed summary; a field allowlist is not
a guarantee that user-entered metadata is non-sensitive. Existing deny-all,
offline and manual-browser-handoff rules remain in force. No LLM is invoked if
assistance is disabled. Explanations never authorize cohort exclusion, model
changes or analysis execution by themselves.

### Claim-adjacent references

Native `study manuscript` revisions use numbered citations in first-appearance
order, reuse the same number for repeated references, and link superscript
brackets to the **References** section at the end of the paper. The reference
list follows Nature-style numerical ordering and bibliographic fields; the
brackets are an explicit display choice, not a claim of exact journal typesetting.
The separate reference catalog retains retrieval provenance and uncited candidates.

Every statement drawn from a journal needs a citation, including Methods,
bioinformatics context beside Results, and Discussion—not only Introduction.
Place the token beside the supported claim, before punctuation:

~~~yaml
citation_policy: strict
sections:
  - id: methods-fdr
    section: methods
    methods_role: discovery_replication
    citation_basis: literature
    text: >-
      We used Benjamini–Hochberg adjustment[[cite:10.1111/j.2517-6161.1995.tb02031.x]],
      with FDR < 0.05.
    sources: [H1.discovery_replication]
    rationale: Cite the performed method without changing the analysis.
    parameters: {method: BH, alpha: 0.05}
~~~

`citation_basis` may be `literature`, `study_results`, `study_design`,
`interpretation`, or `not_applicable`. Mixed prose containing a journal-derived
claim is `literature`. `supporting_evidence` entries have the same citation fields
for their `relationship` text. Multiple references in one token are separated
by semicolons. DOI, `pmid:ID`, and retrieved reference IDs are supported.

Strict mode rejects missing classifications, uncited declared literature prose,
unplaced reference assignments and unresolved identifiers **before creating the
revision output**. Compatibility mode (`citation_policy: review`, the default
for older requests) preserves older prose, displays assigned references and
records gaps in `citation-review.yaml`. It does not certify the draft as ready.
No deterministic checker can establish that a citation supports every sentence
or detect all undeclared borrowing; scientist review is still required.

Acquire missing references through a separate audited native command:

~~~yaml
# references.yaml
schema_version: 1
references:
  - id: 10.1111/j.2517-6161.1995.tb02031.x
    rationale: Cite the multiple-testing method used in discovery and replication.
~~~

~~~powershell
java -jar $jar study references references.yaml --project . --out references-001
Get-FileHash references-001/metadata/manifest.yaml -Algorithm SHA256
~~~

This retrieves exact public DOI/PMID records from Europe PMC and, when a DOI has
no matching Europe PMC record, publisher-deposited Crossref metadata. The source
project's resolved transmission/download policy and `--offline` apply. Only
public identifiers are transmitted; manuscript text, local rationales, study
files and patient rows are not sent. Every request, reason, response, timestamp
and checksum is retained, including blocked or unresolved requests. Retrieval
is not evidence that the article supports a claim. Import the resulting manifest
and its exact checksum in the manuscript request:

~~~yaml
reference_bundles:
  - manifest: references-001/metadata/manifest.yaml
    sha256: REPLACE_WITH_THE_RECORDED_SHA256
~~~

Reference bundles are verified and copied into the new immutable revision.
The manuscript command remains offline; it neither fetches silently nor changes
the source study. Include the saved JAR and acquisition request when archiving
the full research bundle.

Report all effective parameters and their operators explicitly—for example,
`BMI > 30 kg/m2`, not an ambiguous “obesity threshold.” Customary values such as
`FDR < 0.05` need no special justification. List uncommon or deviating choices in
`uncommon_parameters` and supply their recorded reasons in
`parameter_justifications`; missing reasons remain review needs. These fields do
not change the analysis. Dataset tissues are preserved in the cohort inventory:
PBMC and whole blood are different sample compartments. Such a contrast may
contribute to heterogeneity but cannot be isolated from other cohort differences
without an appropriate sensitivity analysis.

## 1. Build and copy the example

From the repository root, build the executable JAR:

~~~powershell
.\gradlew.bat executableJar
$jar = (Resolve-Path "build/cli/jlinalg-0.3.6.jar").Path
~~~

Copy the example to a disposable working directory so the checked-in files stay
clean:

~~~powershell
$project = Join-Path $env:TEMP "jlinalg-study-tutorial"
Copy-Item "examples/study-pipeline" $project -Recurse
Set-Location $project
~~~

If that destination already exists, choose another empty directory. JLinAlg
protects approved contracts and completed artifacts against silent overwrite.

The full lifecycle can be run in one command:

~~~powershell
java -jar $jar study auto study.yaml --project . --actor-id scientist-1
~~~

That command validates input files, resolves defaults, plans, records the
command-line approval, executes or resumes every eligible branch, verifies the
audit chain, and rebuilds the indexes. The remaining sections unpack the same
lifecycle so each boundary is visible.

The input contract is `study.yaml`; its two tabular inputs are
`data/discovery.tsv` and `data/replication.tsv`. The contract explicitly declares
that the cohorts are independent, so this example has no pending relationship
assumption.

## 2. Confirm beginning state and validate inputs

A missing `state/project-status.yaml` means that the project is at its beginning,
not that an analysis completed without a status file:

~~~powershell
java -jar $jar study status --project .
java -jar $jar study validate study.yaml --check-files
~~~

`validate` checks the versioned structure, IDs, cross-references, analysis goals,
and—because `--check-files` is present—the declared paths. It reads but does not
resolve or execute the project.

Before continuing, inspect the input tables. Each row has a feature ID, effect,
standard error, and p-value. A result-first meta-analysis requires compatible
effect definitions and units; syntactically compatible columns are not evidence
that two scientific estimands are interchangeable.

If one physical table contains several phenotypes, declare one logical dataset
per phenotype and reuse the file path with different column mappings. The
[result-first contract vignette](auditable-omics-pipeline.md#multi-phenotype-result-files-and-strict-feature-universes)
shows how to handle R row names, invalid-row exclusion, a strict cross-cohort
feature intersection, discovery-subset FDR, and direction-independent validation.

## 3. Resolve defaults and inspect the contract

Resolve the human contract into an immutable, fully explicit contract:

~~~powershell
java -jar $jar study resolve study.yaml --actor-id scientist-1
~~~

This creates:

- `resolved-study.yaml`, including expanded defaults and assumptions;
- `ledger/events.jsonl`, the append-only machine-readable event history;
- `SESSION_HISTORY.md`, the readable research notebook;
- `state/project-status.yaml`, a materialized progress view; and
- `AGENTS.md`, a derived working brief for a person or future agent.

Review all five. `study.yaml` remains the human-authored source, while
`resolved-study.yaml` is the exact contract from which the plan is compiled.
Directly editing the base study after resolution cannot replace the approved
contract; project changes use amendments.

## 4. Compile and approve the plan

Compile a deterministic operation graph and inspect it before approval:

~~~powershell
java -jar $jar study plan --project . --actor-id scientist-1
java -jar $jar study status --project .
Get-Content plans/*/analysis-plan.yaml
~~~

The plan identifies each operation, its datasets, dependencies, executor,
configuration, expected outputs, and registered diagnostics. Planning does not
run an analysis.

Approval is a separate human event:

~~~powershell
java -jar $jar study approve --project . --approved-by scientist-1
~~~

If another project contains default independence assumptions, the approval
command stops until they are reviewed. `--accept-assumptions` is an explicit
acknowledgment; it should not be added merely to make a command pass.

## 5. Execute the approved operations

Run the plan:

~~~powershell
java -jar $jar study run --project . --actor-id scientist-1
~~~

The example completes five operations. Under
`runs/run-<plan-id>/operations/`, inspect:

- standardized association tables and their source hashes;
- `replication.tsv`, including adjusted discovery and replication p-values,
  sign concordance, and the final `replicated` indicator;
- `meta-analysis.tsv`, its log, and captured console output;
- `diagnostics.yaml`, including cross-cohort correlation, OLS scale/slope,
  sign concordance, discovery-selected and replicated counts, meta-analysis
  FDR/heterogeneity counts, the unique count replicated in any validation
  cohort, and the current unusual-observation list; and
- `report.md`, whose result IDs link the report to operation manifests.

Every operation directory also has a `manifest.yaml` with its configuration,
inputs, output hashes, timestamps, executor, and state. The run-level
`run-manifest.yaml` summarizes the operation states. The exact directory names
come from the plan and should be discovered from the manifest rather than
hard-coded in a downstream script.

For `mode: auto`, functional and genetic follow-up branches are independent.
JLinAlg runs local GO and Reactome enrichment when feature mapping is available;
`annotation_release: GENCODE_v48`, for example, permits an audited public
GENCODE transcript-ranking and HGNC-symbol mapping download. Public references
are fetched by the HTTP client embedded in the executable JAR; curl, Python,
and a separate downloader are not required. If eligible variant-level work requires LD, the
default is the checksum-verified `1000g-phase3` package with five ancestry
panels. It is downloaded only when that branch is eligible, and
`resources.ld.population: auto` waits for an ancestry choice rather than
silently assuming EUR. A missing regional molecular-QTL or outcome-GWAS table
becomes a precise `request.yaml`; it does not block enrichment or reporting.

## 6. Verify and safely resume

Verify the evidence chain:

~~~powershell
java -jar $jar study verify --project .
java -jar $jar study run --project . --resume --actor-id scientist-1
java -jar $jar study verify --project .
~~~

Verification recalculates hashes for the base and resolved contracts, plans,
amendments, observations, derived notebook files, manifests, and operation
outputs. Resume reuses only verified complete artifacts. A changed output is an
error, not a cache hit.

`--resume` is now optional. Both `study run` and `study auto` resume by default.
Every invocation retains a separate immutable attempt, reusing only completed
operations whose input hashes and implementation identity match. Waiting and
failed operations are retried; changed upstream outputs invalidate dependants.
The index shows the latest attempt by default and offers a historical-attempt
toggle; historical blockers do not count as current blockers.

Use `study auto study.yaml --rerun` for a complete recomputation, bypassing result
reuse without deleting earlier evidence. `--overwrite` is deliberately not used.
`--offline` (or `resources.network: offline`) prevents built-in network requests;
available local work continues. Missing resources and blocked providers are
shown in the index. Repeating without `--offline` retries unfinished work.
An HTTP 403 is recorded as an endpoint/authorization denial, not asserted to be
a firewall diagnosis. A timeout or failed search never becomes “no association”
or “no literature.” No concurrent execution is allowed in the same project;
do not amend it while a run is active.

The rebuilt index has four complementary entry points:

- `index/index.html` for interactive search;
- `index/index.json` for artifact-level retrieval;
- `index/evidence-summary.json` for a portable machine handoff; and
- `index/LLM_CONTEXT.md` for paper drafting or drug-prioritization retrieval.

The LLM context embeds the compact deterministic diagnostics and the leading
GO, Reactome, or candidate-prioritization rows, while retaining the result ID,
artifact path, and SHA-256 needed to verify each value. It deliberately does not
embed entire result tables. A drafter should retrieve additional rows from the
indexed artifact rather than treating the context file as a new source.

The last file is an index, not a source of truth. Any LLM-generated claim must
follow its result ID to the hashed result table and operation manifest.

The command exit statuses are:

| Status | Meaning |
| --- | --- |
| `0` | Requested work completed or verified. |
| `2` | Validation, integrity, or execution failed. |
| `3` | The plan is waiting for an explicitly configured operation or resource. |

`study run --allow-partial` converts the last case to process success for batch
schedulers, but the project remains visibly `waiting_for_input`.
`study auto` instead treats documented scientific prerequisite gaps as a
successful completion of all currently available work and returns nonzero only
for validation, execution, or integrity failures.

## 7. Make a study change through an amendment

Copy the current `resolved_study_hash` from `study status` into an amendment.
An amendment contains only the change and its rationale:

~~~yaml
schema_version: 1
amendment:
  id: sensitivity-2026-10-01
  title: Add a prespecified heterogeneity sensitivity analysis
  base_study_hash: "<resolved_study_hash from study status>"
  proposed_at: 2026-10-01T16:30:00Z
  proposed_by:
    type: human
    id: scientist-1
  source:
    type: meeting
    id: methods-meeting-2026-10-01
  rationale: Test whether the pooled finding is sensitive to cohort heterogeneity.
changes:
  add:
    analysis_requests:
      - id: SA1
        type: sensitivity_analysis
        for_hypothesis: H1
        description: Compare fixed and random-effects pooling.
        reason: Cohort effects may not share one common effect.
~~~

Preview before approval:

~~~powershell
java -jar $jar study amend amendment-sensitivity-2026-10-01.yaml --project . --check
java -jar $jar study amend amendment-sensitivity-2026-10-01.yaml --project . --approved-by scientist-1
~~~

Applying an amendment supersedes the active plan. Compile, review, and approve a
new plan before more execution. Prior contracts, plans, runs, and rejected or
retired directions remain in the evidence history.

## 8. Record an unusual result

Suppose a generated effect-comparison plot contains an off-diagonal streak.
Keep the plot inside the project and record the observation without deciding in
advance whether it is biology or artifact:

~~~powershell
java -jar $jar study observe --project . `
  --id off-diagonal-1 --type off_diagonal_streak `
  --artifact plots/effects.png --hypothesis H1 `
  --description "Off-diagonal effect streak" `
  --explanation "processing batch" `
  --explanation "biological subgroup" `
  --proposed-direction sensitivity_analysis `
  --amendment-id check-batch `
  --actor-type human --actor-id scientist-1
~~~

The observation stores the plot hash, competing explanations, proposed
direction, hypothesis, and actor. If a local or remote LLM drafted the review,
use `--actor-type llm` and identify the provider/model or local model hash with
`--model`. That is provenance only: the proposed direction still needs a human-
approved amendment.

## Documented scientist questions

The source build now records deterministic questions in `QUESTIONS.md` and
`SESSION_HISTORY.md`, with their exact wording, reason, response options,
hypothesis/operation scope, input hashes, UTC timestamp, rule, and actor.
`state/questions.yaml` is the machine-readable view; `index/index.html` makes
questions and answers searchable alongside the analysis evidence. Documented
means available for review, not proof that someone read the question.

The small result-first example may need no questions. A matrix study without a
phenotype file/outcome, a shared-parent cohort comparison, or a genetic follow-up
without feature annotation or GWAS/QTL inputs will exercise the first rules.
Missing core required study fields still cause validation errors; general
incomplete-contract intake is not implemented. Unsupported older blockers remain
visible in operation requests even if they have no question rule yet.

In the project directory, inspect questions and export a blank answer batch:

~~~powershell
java -jar $jar study questions --project .
java -jar $jar study questions --project . --out answers-meeting-01.yaml
~~~

No provider call or analysis runs, and an existing export is never overwritten.
Open `QUESTIONS.md` beside the answer file. Retain the exported `project_id`,
question ID, revision, and hash. Choose your own batch `id`. Fill only responses
you are ready to record; remove the other response entries from this batch.
For example, a factual phenotype-file response has this shape (the IDs/hash
below are placeholders and must come from your export):

~~~yaml
schema_version: 1
id: meeting-01-answers
project_id: my-project
amendment_id: phenotype-metadata-01
meeting_id: lab-meeting-01
responses:
  - question_id: Q-EXPORTED-ID
    question_revision: 1
    question_hash: EXPORTED-HASH
    disposition: answer
    value: metadata/phenotypes.tsv
    rationale: The data manager confirmed this is the released phenotype file.
    evidence:
      - file: metadata/README.txt
        selector: Release and phenotype definitions
        note: Source description reviewed by the scientist.
~~~

The phenotype/annotation resource must exist inside the project. Its hash is
recorded without copying patient-level contents into the answer notebook.
Evidence `file` entries, in contrast, are copied into immutable local snapshots:
use short UTF-8 metadata documents, not patient records (maximum 1 MiB each).
A `citation` or explicit source `note` can be used without a file; citations
are recorded as scientist-supplied, not automatically verified.

Every response needs a rationale. `disposition` is `answer`, `unknown`, `defer`,
or `decline`; `value` must be absent/null for the last three. Mandatory
prerequisites cannot be declined. Unknown/deferred questions remain unresolved
and do not unblock an analysis. Declining optional review records work not
pursued, not a negative finding. Recommendations are never selected by default.

Check the complete batch, then record it explicitly:

~~~powershell
java -jar $jar study answer answers-meeting-01.yaml --project . --actor-id scientist-1 --check
java -jar $jar study answer answers-meeting-01.yaml --project . --actor-id scientist-1 --record
~~~

`--check` writes nothing. `--record` validates the entire batch before recording
it, preserves the submitted file, and appends a ledger event. Both require an
explicit local human actor. A local name is attribution, not authenticated
identity. Neither command runs an analysis or contacts an LLM.

Supported factual changes create a proposed amendment, not an immediate change
to `resolved-study.yaml`. Include a user-chosen `amendment_id` for those batches;
omit it when recording only unknown/defer/decline or narrative review responses.
The CLI receipt names the proposal, when one exists:

~~~powershell
java -jar $jar study amend amendments/proposed/amendment-phenotype-metadata-01.yaml --project . --check
java -jar $jar study amend amendments/proposed/amendment-phenotype-metadata-01.yaml --project . --approved-by scientist-1
java -jar $jar study auto study.yaml --project . --actor-id scientist-1 --offline
~~~

Inspect the amendment impact before approval. A narrative response about GWAS/QTL
source selection or diagnostic follow-up is documented but does not generate an
arbitrary analytical plan: prepare a supported amendment and any required local
inputs separately. Related/overlapping cohorts still require an appropriate
analysis; acknowledging dependence or passing `--accept-assumptions` does not
make independent-cohort meta-analysis valid.

To correct an unresolved answer, export with `--include-answered`, use a new
batch ID, and retain its generated `supersedes` reference to the previous answer.
Old answers remain auditable. Repeating an identical recorded batch is safe;
reusing its ID for changed content is rejected. Changed study/input hashes or
question revisions reject stale answers rather than silently applying them.
After an amendment the existing cache conservatively re-evaluates all operations
under the changed study hash.

Questions, answers, and small evidence snapshots are checked by `study verify`.
If only a derived notebook/status file was lost, this can rebuild the views:

~~~powershell
java -jar $jar study verify --project . --repair-views
~~~

Repair refuses missing or modified authoritative records. It does not alter the
ledger or recover corrupt evidence. The local HTML/JSON index contains answer
text and must be treated as sensitive. Raw answers and evidence snapshots are
omitted from `LLM_CONTEXT.md` and existing assistance previews; no file becomes
approved for upload simply because it is called an index.

See the [implemented scope and remaining design](study-questions-and-decisions.md#implemented-first-tranche).
LLM question generation, interactive prompting, authenticated remote approval,
and generated-code execution remain future work.

## Patient-level data and LLM non-egress

### What can be claimed today

LLM calls are opt-in. The adapter never uploads an input file or gives the model
filesystem, shell, sample-matrix or phenotype-row access. It constructs a small
allowlisted result preview; raw feature identifiers, project descriptions,
dictionary prose and local paths are omitted. Gene symbols, numeric aggregate
results, phenotype labels and retrieved reference metadata can be transmitted.
These are still sensitive unpublished findings. The approved input metadata
must correctly classify results; this is not general-purpose PHI detection.
`deny_all` blocks every remote advisory request, including topic-only requests,
and `--offline` also blocks built-in reference downloads.

JLinAlg cannot truthfully prove that unrelated software, a separately launched
agent, a shell command, or a user copied no data elsewhere. A defensible project
claim therefore needs both application evidence and enforcement outside the
application.

### A controlled aggregate-results disclosure design

A future LLM integration should never give a model process unrestricted access
to the study directory. Instead, JLinAlg should implement a release boundary:

~~~text
restricted project data
        |
        v
artifact classification + identifier/DLP checks + small-cell rules
        |
        v
separate disclosure bundle (aggregate, redacted artifacts only)
        |
        v
human approval + immutable disclosure manifest
        |
        v
sandboxed connector with read-only bundle and allowlisted network destination
        |
        v
LLM response, response hash, and provider/model attribution
~~~

The implemented adapter uses `aggregate_results` as the default eligibility
policy when explicitly enabled: bounded result previews may be disclosed, but
patient/sample-level data may not. This is permission for a controlled artifact
class, not blanket permission to upload every file in `results/` and not an
instruction to transmit results when no LLM operation needs them.

For this policy, an aggregate result is identified by its registered artifact
type, schema, content checks, and lineage—not by its filename. Eligible examples
include association and meta-analysis summary tables, enrichment and fine-mapping
results, summary diagnostics, and redacted figures. Participant matrices,
phenotype/covariate rows, sample identifiers, individual residuals or predictions,
and per-person scores remain prohibited even if a tool wrote them to a result
directory. Prompts and study metadata are also manifested because unpublished
hypotheses and findings can be sensitive despite containing no patient rows.

A project may set `deny_all` so no project content is transmitted to a remote
LLM. This is appropriate for embargoed, commercially sensitive, regulated, or
otherwise restricted unpublished work. The additional controls below are a
future high-assurance contract, not all implemented today:

~~~yaml
project:
  data_policy:
    classification: restricted_patient_level
    llm:
      remote_transmission: aggregate_results  # default; or deny_all
      remote_patient_level_data: deny
      eligible_remote_artifacts:
        - association_results
        - meta_analysis_results
        - enrichment_results
        - fine_mapping_results
        - summary_diagnostics
        - redacted_figures
      require_human_approval: true
      minimum_cell_count: 10
~~~

`remote_transmission` denial and patient-level exclusion are implemented for the
bounded adapter. Per-request human approval, artifact-specific disclosure lists,
small-cell checks, and OS isolation remain future work. If
`require_human_approval: true`, `eligible_remote_artifacts`, or
`minimum_cell_count` is present, the current adapter refuses to transmit and
records `additional_disclosure_review_required`; it does not silently ignore
these stricter policies. Do not remove a required control merely to get a call
to run. Application logs are not an independent non-egress guarantee.

For `aggregate_results`, the system should create an immutable
`llm-disclosure-manifest.yaml` before each request. It should record:

- hashes and media types of the exact bytes released;
- the source result IDs and transformations used for redaction or aggregation;
- direct-identifier, quasi-identifier, free-text, and minimum-cell checks;
- destination, provider, reported model, purpose, timestamp, and approver;
- the policy version and pass/fail decision; and
- request and response hashes, without duplicating restricted content into an
  unrestricted log.

Calling a table “aggregate” or “de-identified” is not by itself proof of
negligible re-identification risk. For regulated health data, the project should
tie release criteria to the applicable legal and institutional process—for
example, HIPAA Safe Harbor or Expert Determination in the United States—and
retain the determination as evidence. Small cells, rare combinations, genomic
variants, dates, geography, and free text require particular care.

### Independent enforcement and verification

For high-assurance studies, run the connector in a separate OS account,
container, or sandbox that can read only the approved disclosure directory.
Use default-deny outbound networking with an allowlist for the selected provider,
and retain proxy, DNS, firewall, or equivalent egress logs. Here, default-deny
describes connector filesystem and network capabilities; it does not change the
project-policy default that verified aggregate results are eligible. A local
model should run with networking disabled and record hashes for its weights,
runtime, prompt, and supplied artifacts.

Assurance should be tested, not inferred from configuration alone:

1. Review the connector and its transitive dependencies; publish an SBOM and
   signed/reproducible release artifacts where practical.
2. Run integration tests with canary identifiers in forbidden files and verify
   that requests, logs, crash reports, and telemetry never contain them.
3. Verify that the connector cannot read outside the disclosure directory and
   cannot reach an unapproved endpoint.
4. Reconcile disclosure manifests with network logs and the research ledger.
5. Periodically obtain an independent security/privacy review for regulated or
   high-consequence deployments.

Provider agreements, retention controls, training opt-outs, and a BAA or DPA may
be necessary, but they complement rather than replace technical data
minimization and egress enforcement. If an external tool cannot supply adequate
evidence, the project should record a provenance and privacy-assurance gap rather
than claim end-to-end non-egress.

## Current boundaries

- Direct deterministic adapters currently cover result standardization,
  bulk-matrix association, discovery/replication comparison, configurable
  fixed/random meta-analysis, GENCODE feature mapping, local GO/Reactome
  enrichment, diagnostics, transparent candidate evidence tables, and
  traceable reporting.
- Configured downstream operations can call declared tools; absent tools or
  inputs produce a hashed request instead of an invented result. Eligible
  LD-dependent work automatically resolves the default 1000 Genomes reference.
- Optional provider adapters support bounded advice, retrieval-cited literature
  review and explicitly permitted unreviewed HTML drafts. No provider was called
  while testing this feature: provider tests use simulated transport responses.
- Raw single-cell, spatial, microbiome, and multimodal data wait for a registered
  modality-specific operation rather than being analyzed as a generic bulk
  matrix.
- Full disclosure-bundle isolation, DLP/small-cell certification, general pipeline
  local models, MCP and provider-specific OAuth login/refresh remain unimplemented.
  Browser handoff is the default transport when assistance is enabled; optional
  API adapters accept API keys or provider-issued OAuth access tokens. The separate
  manuscript editorial command supports exact-payload approval and explicit loopback models;
  see the [bioinformatics vignette](bioinformatics-evidence.md). No regulatory
  compliance or publication-readiness certification is claimed.

For the full contract, amendment, provenance, and operation-registry design, see
the [auditable omics research pipeline vignette](auditable-omics-pipeline.md).
The runnable files are in the
[study-pipeline example](https://github.com/robbyjo/JLinAlg/tree/main/examples/study-pipeline).

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [Geir Kjetil Sandve et al. (2013) — Ten Simple Rules for Reproducible Computational Research](../CITATIONS.md#sandve-reproducibility-2013)
- [W3C Provenance Working Group (2013) — PROV-O: The PROV Ontology](../CITATIONS.md#w3c-prov-o-2013)
- [Adam Frankish et al. (2021) — GENCODE 2021](../CITATIONS.md#frankish-gencode-2021) — [PMID: 33270111](https://pubmed.ncbi.nlm.nih.gov/33270111/) · [PMCID: PMC7778937](https://pmc.ncbi.nlm.nih.gov/articles/PMC7778937/)
- [Gene Ontology Consortium (2021) — The Gene Ontology resource: enriching a GOld mine](../CITATIONS.md#go-consortium-2021) — [PMID: 33290552](https://pubmed.ncbi.nlm.nih.gov/33290552/)
- [Marija Milacic et al. (2024) — The Reactome Pathway Knowledgebase 2024](../CITATIONS.md#milacic-reactome-2024) — [PMID: 37941124](https://pubmed.ncbi.nlm.nih.gov/37941124/) · [PMCID: PMC10767911](https://pmc.ncbi.nlm.nih.gov/articles/PMC10767911/)
- [1000 Genomes Project Consortium (2015) — A global reference for human genetic variation](../CITATIONS.md#thousand-genomes-2015) — [PMID: 26432245](https://pubmed.ncbi.nlm.nih.gov/26432245/) · [PMCID: PMC4750478](https://pmc.ncbi.nlm.nih.gov/articles/PMC4750478/)
- [U.S. Food and Drug Administration (2024) — Electronic Systems, Electronic Records, and Electronic Signatures in Clinical Investigations: Questions and Answers](../CITATIONS.md#fda-electronic-records-2024)
- [U.S. Department of Health and Human Services, Office for Civil Rights (2012) — Guidance Regarding Methods for De-identification of Protected Health Information in Accordance with the HIPAA Privacy Rule](../CITATIONS.md#hhs-hipaa-deidentification-2012)
- [National Institute of Standards and Technology (2024) — Artificial Intelligence Risk Management Framework: Generative Artificial Intelligence Profile](../CITATIONS.md#nist-genai-profile-2024)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
