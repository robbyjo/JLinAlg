# Gene-set, disease and phenotype enrichment

The source build provides Java APIs and a command-line enrichment workflow with
an analysis-specific background. It supports ordinary over-representation
analysis (ORA) and the equivalent-CpG/fractional-count `gsameth` method for EWAS.
Database downloads are explicit; analyses use local files and do not contact
annotation services. Build with `./gradlew check executableJar` (Windows:
`.\gradlew.bat check executableJar`). These features may postdate the latest
published release JAR.

## A complete example

The repository includes a small synthetic results table, probe annotation and
GMT collection in `examples/enrichment`. Run from the repository root:

```shell
java -jar build/cli/jlinalg-0.3.5.jar --enrichment Custom --enrichment-db examples/enrichment/sets.gmt --input examples/enrichment/results.csv --input-id probe --universe-selection "status == 'OK' && n_studies >= 3" --selection "p_value < bonferroni(0.05)" --annot examples/enrichment/annotation.csv --annot-id IlmnID --gene-col genes --annot-cols chromosome --gene-id-type symbol --min-set-size 1 --out build/enrichment-example/ora.tsv
```

The six eligible probes define the feature-level testing family. Bonferroni uses
`0.05 / 6`, before mapping to genes. One eligible probe has no gene mapping; it
still counts toward that threshold. Two probes are selected, mapping to three
unique genes. Failed or insufficient-cohort probes are excluded from both
selection and background. The two eligible terms include a zero-hit term.
The small size filter is for this demonstration; the default is 10–500 genes.

Four files are written: the results table, `.log`, `.features.tsv`, and
`.manifest.properties` appended to its filename. Existing outputs are refused.
The feature audit records each eligible probe, selection status, gene mappings
and requested annotation columns. The manifest records source/input hashes,
settings, counts and the correction family. The screen and log announce defaults.

## Define the universe before looking at enrichment

With `--selection`, `--input` must be the complete results table for the chosen
analysis, not just its significant rows. Its eligible IDs define the universe.
For a meta-analysis, use the actual success/minimum-cohort requirements:

```shell
java -jar build/cli/jlinalg-0.3.5.jar --enrichment GO:BP --enrichment-db databases/go-human --input meta.csv --input-id probe_id --universe-selection "status == 'OK' && n_studies >= 3" --selection "FDR < 0.05" --annot annotation.csv --annot-id probe_id --gene-col gene_symbol --gene-id-type symbol --out go-bp.tsv
```

Adjust column names to the table you actually produced. Eligibility must reflect
which hypotheses could have been selected in that analysis. Select one phenotype
and contrast before running enrichment; duplicate eligible feature IDs are an
error. Rows omitted from a supplied table cannot be recovered automatically.

`--input-cols` is an alias for the single `--input-id` column. Annotation joins
are explicit. `--gene-col` names a column in `--annot`, or in the universe table
when no annotation file is supplied. Multiple genes in a cell are separated by
literal `;` by default (`--gene-separator`). Repeated edges and genes are
deduplicated. `--annot-cols all` retains all annotation columns in the audit.

Map input and database IDs to the **same namespace and species**. The tool does
not guess aliases or convert symbols, Ensembl IDs, Entrez IDs and CURIEs.
`--strip-gene-version true` explicitly strips a terminal `.digits` from input
gene IDs. It does not convert namespaces. The database manifest declares its
namespace, and a conflicting `--gene-id-type` is rejected. A zero total database
match is an error; partial mapping losses are reported.

The gene universe is the union of genes mapped from all eligible probes. Genes
without membership in this particular collection remain in the universe; they
are not silently discarded. This is an explicit background convention, which
can differ from tools that intersect the background with annotated genes. Keep
the convention fixed across comparisons. Features without any gene mapping are
reported and omitted from gene-level testing. They still belong to the original
feature-level selection family.

If input is already a selected-only list, supply a complete, analysis-eligible
background explicitly. Headerless files contain one ID per line:

```shell
java -jar build/cli/jlinalg-0.3.5.jar --enrichment Custom --enrichment-db sets.gmt --input selected-genes.txt --background tested-genes.txt --out enrichment.tsv
```

Headered tables require `--input-id` and `--background-id`, respectively. Selected
IDs must belong to the eligible background. `--selection` and `--background`
are mutually exclusive. There is no whole-genome default background.

## Selection expressions and missing results

Quote the whole expression for your shell. Supported expressions include:

| Expression | Meaning |
| --- | --- |
| `FDR < 0.05` | Use an existing input column; no recomputation of that FDR |
| `p_value < bonferroni(0.05)` | Alpha divided by evaluable eligible feature count, before gene mapping |
| `p_value < bonferroni` | Same, with alpha 0.05 |
| `p_value < bonferroni(0.05, 850000)` | Use an explicitly supplied original testing-family size |
| `p_value < 0.01 && abs(beta) > 0.2` | Joint significance and effect-size selection |
| `is_finite(p_value) && p_value < 0.05` | Explicitly guard numeric evaluation |

The language supports parentheses, `&&`, `||`, `!`, comparisons, arithmetic,
`abs()`, `is_missing()` and `is_finite()`. Column names containing spaces can be
enclosed in backticks; string values use single or double quotes. It is a small
parser, not an R, Java or shell interpreter. Names are case-sensitive. Logical
operators short-circuit; missing comparisons otherwise yield an unknown result.
Blank, `NA`, `NaN`, `NULL` and `.` are missing. Invalid numeric text and arithmetic
such as division by zero fail explicitly.

By default, an unevaluable eligibility/selection expression is an error.
`--missing-selection exclude` excludes unevaluable rows before resolving the
default Bonferroni count. Prefer an explicit eligibility rule for failed analyses.
For example, put `is_finite(p_value)` in `--universe-selection` if nonfinite tests
could not have been selected. Putting that guard only in `--selection` instead
keeps those rows as eligible, unselected background. This distinction is deliberate.

The existing input `FDR` and output enrichment `FDR` belong to different hypothesis
families. If the upstream correction family included hypotheses absent from this
file, supply its original count to `bonferroni(alpha, m)` or use its adjusted
column. A permissive selection such as `FDR < 0.5` is accepted, but should be
scientifically justified before inspecting the enrichment results.

## Statistical methods

### Ordinary ORA: the default

`--enrichment-method ora` tests over-representation with the inclusive upper
hypergeometric tail, equivalent to a one-sided Fisher exact test. For universe
size N, selected size n, background term size K and overlap k, the table is:

| | In term | Outside term |
| --- | ---: | ---: |
| Selected | k | n − k |
| Unselected | K − k | N − K − n + k |

The p-value is `P[X >= k]` for `X ~ Hypergeometric(N, K, n)`. This is a competitive
null: conditional on n, genes in the chosen universe are equally likely to be
selected. A negative-binomial count model is not the default for this question.
It models a different sampling process and does not automatically address gene
selection bias or overlapping terms.

ORA reports expected overlap `n*K/N`, fold enrichment `k/(n*K/N)`, and the
cross-product odds ratio. Boundary odds ratios can be infinite or undefined;
undefined values are written as `NA`, with no pseudocount added. Empty selections
produce p-values of one. Hypergeometric and FDR calculations preserve log
probabilities, so a displayed double p-value of zero may still have a finite
`log_p_value`.

### EWAS: `gsameth`

CpGs do not give every gene the same chance of selection: genes differ in probe
coverage, and individual probes can map to several genes. For an EWAS, explicitly
choose `--enrichment-method gsameth` and provide the normalized, platform-appropriate
probe-to-gene annotation used for **all eligible CpGs**:

```shell
java -jar build/cli/jlinalg-0.3.5.jar --enrichment Custom --enrichment-db examples/enrichment/sets.gmt --enrichment-method gsameth --input examples/enrichment/results.csv --input-id probe --universe-selection "status == 'OK' && n_studies >= 3" --selection "p_value < bonferroni(0.05)" --annot examples/enrichment/annotation.csv --annot-id IlmnID --gene-col genes --min-set-size 1 --out build/enrichment-example/ewas.tsv
```

The implementation follows missMethyl's `equiv.cpg=TRUE`, `fract.counts=TRUE`,
`prior.prob=TRUE` path:

1. A probe mapped to r unique genes contributes `1/r` equivalent CpGs to each.
2. Selected probes contribute the same fractions, capped at one per gene.
3. Genes are ordered by equivalent coverage. A rank-based tricube moving average
   with span 0.5 estimates selection probability from binary gene selection.
   Coverage ties use lexicographic gene ID order (R comparison uses `LC_COLLATE=C`).
4. A term's bias odds are its mean fitted probability divided by the outside-term
   mean. The two-colour Wallenius noncentral hypergeometric distribution tests
   its weighted overlap. The draw count is the number of selected unique genes.
   Fractional overlap is truncated to an integer for the inclusive tail, matching
   the BiasedUrn interface used by the R reference.

Wallenius probabilities are computed by deterministic successive-sampling dynamic
programming. This is distinct from Fisher's noncentral hypergeometric distribution.
There is no sampling, normal approximation or fallback to an ordinary central
hypergeometric tail when a biased tail fails numerically. Interior underflow fails
explicitly. The central odds=1 case uses the central distribution directly.
The default per-term work cap is 100 million state updates; increase it explicitly
with `--wallenius-max-updates` if needed. The method can be costly for very large
selected lists and terms. The smoother is quadratic in the number of genes.

The `*_unweighted` effect summaries remain descriptive ordinary gene-count
summaries; `weighted_overlap` and `selection_odds` identify the inputs to the
biased test. This method adjusts probe-coverage and multi-gene mapping bias.
It does not eliminate arbitrary CpG/gene correlation, annotation uncertainty or
all selection-model misspecification. Its estimated probability-weighting null
does not turn downstream FDR guarantees into unconditional guarantees for every
EWAS design. Platform mapping and upstream QC remain the analyst's responsibility.

See the [GOmeth/GOregion methods paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8186068/)
and [missMethyl source](https://github.com/bioc/missMethyl/tree/c5d0a518aaacd98c9f31c4b6ca48239246363ad4).

### GWAS and TWAS boundaries

GWAS Catalog enrichment here means **overlap of a supplied gene list with Catalog
trait-associated gene sets**. The default Catalog membership is the deduplicated
**union of reported and mapped genes**. It is announced on screen and in the log;
`--catalog-gene-source reported` or `mapped` selects one source explicitly.
Catalog associations pass `p <= 5e-8` by default, configurable with
`--catalog-p-threshold`. Invalid or missing p-values are rejected and counted.
Trait identifiers default to mapped ontology IDs (`--trait-id-type efo`);
`label` instead uses reported trait labels. EFO mode retains other mapped
ontology identifiers too. Raw source files preserve association provenance.

Union is a transparent inclusion policy, not a causal-gene assignment method.
SNP-to-gene mapping, LD, gene size and correlated molecular predictors can bias
ordinary overlap tests. This CLI does not turn GWAS summary statistics into an
LD-aware competitive pathway test, and it does not claim MAGMA, stratified LDSC,
CAMERA or phenotype-permutation inference. Likewise, ordinary TWAS gene overlap
does not adjust correlated prediction models. Use an appropriate upstream
gene-level analysis and interpret the specific overlap null being tested.

## Overlap, hierarchy and FDR

The default `--fdr-type BH` applies Benjamini–Hochberg to **every term passing the
background-only size filters**. Terms with zero selected hits stay in the family
with p=1. Sizes are measured after intersecting each term with the universe.
Filtering by observed hit count before correction would change this family and
is intentionally not provided. `GO:BP`, `GO:MF`, `GO:CC` each define a separate
requested family; `GO` tests all three aspects together. Separate invocations
do not automatically control one combined family across all databases or traits.

BH has its standard guarantee under independence or suitable positive dependence.
Overlap alone does not prove that these assumptions hold. `--fdr-type BY` uses
the Benjamini–Yekutieli harmonic factor for arbitrary dependence, provided the
individual null p-values are valid; it is usually more conservative. Both are
monotone adjustments performed in log space. `NONE` returns raw p-values in the
adjusted-value columns and is labeled in the run metadata.

GO propagation follows `is_a` and `part_of` in `go-basic.obo`, excluding `NOT`
annotations and obsolete terms; alternate term IDs are resolved. It does not
propagate `regulates`. HPO/Mondo propagation stays within the requested ontology's
term namespace. Cycles in the applicable hierarchy fail explicitly. GO/HPO/Mondo
uninformative root IDs are excluded independently of selection. `--propagate false`
adds no ancestors, but cannot undo ancestry already present in a source such as
Reactome's all-level mappings. `--evidence EXP,IDA,IPI` optionally filters GAF
evidence; the default retains all positive source evidence codes.

These are flat term-level FDR procedures. Hierarchy-aware testing, Focused BH,
DAGGER, `elim`/`weight`, and redundancy-aware discovery families are not implemented
in this version. Showing a few representative significant terms does not inherit
the same FDR claim for that filtered subset. Retain the full corrected table,
and describe any later term clustering as interpretation. See the
[R BH/BY documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)
and [Focused BH paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10281705/).

## Download and reuse annotation databases

Choose a **new directory for each snapshot**. Files are staged, checked, hashed
with SHA-256 and parsed before installation. A failed download does not publish
an installed package. Loading a package verifies its source hashes. Retain the
whole directory and manifest for reproducibility; `latest` is resolved at download
time where the provider exposes a release identifier, and otherwise the retrieval
timestamp, URLs, HTTP metadata and hashes identify the snapshot.

| Selector | Built-in source and default gene IDs | Access and scope |
| --- | --- | --- |
| `GO`, `GO:BP`, `GO:MF`, `GO:CC` | GO release, go-basic OBO + organism GAF; symbols, optionally UniProt | Public downloads; automatic human/mouse/rat; other taxa through explicit source import |
| `Reactome` | All-level pathway mappings and hierarchy; Ensembl, optionally Entrez/UniProt | Public downloads; filter species; current or explicit release |
| `HPO` | HPO release, genes-to-phenotype + hp.obo; Entrez | Public human phenotype associations, not a disease-gene ontology by themselves |
| `GWASCatalog` (alias `GWAS`) | Catalog ontology-annotated associations; symbols | Public gene-trait overlap source; reported/mapped union by default |
| `Monarch:Mondo` | Monarch causal gene-to-disease associations + Mondo OBO; CURIEs such as `NCBIGene:7157` | Public downloads; `--monarch-evidence correlated` chooses a different source package |
| `OpenTargets` | Direct overall target-disease associations and disease labels; Ensembl | Public release Parquet; default overall score >=0.5; may download over 1 GB |
| `MSigDB:<collection>` | User-provided GMT, e.g. `MSigDB:H` | Free registration/terms apply; obtain the file first, then import; an already accessible source URL can also be used |
| `KEGG` | Pathway links and names; Entrez for human/mouse/rat | Public REST is restricted to academic users; explicit `--kegg-academic true` or authorized local import |
| `Custom` | GMT or headered term-to-gene table | User-provided file or accessible HTTPS URL |

Other than KEGG, several sources offer public automatic downloads, but **MSigDB
still requires its registered download workflow**. An ontology file alone (for
example Mondo or Disease Ontology) does not supply gene memberships. Use an
association source with that ontology, or import a compatible custom gene set.
Do not treat disease association evidence as proof of a causal gene. Open Targets
overall scores combine evidence and the threshold is an explicit analyst choice,
not a significance test. `--association-min-score` changes that threshold.

```shell
java -jar build/cli/jlinalg-0.3.5.jar --enrichment GO --download databases/go-human --species human --gene-id-type symbol
java -jar build/cli/jlinalg-0.3.5.jar --enrichment Reactome --download databases/reactome-human --gene-id-type ensembl
java -jar build/cli/jlinalg-0.3.5.jar --enrichment HPO --download databases/hpo-human
java -jar build/cli/jlinalg-0.3.5.jar --enrichment GWASCatalog --download databases/gwas-catalog
java -jar build/cli/jlinalg-0.3.5.jar --enrichment Monarch:Mondo --download databases/monarch-mondo
java -jar build/cli/jlinalg-0.3.5.jar --enrichment OpenTargets --download databases/open-targets --association-min-score 0.5
java -jar build/cli/jlinalg-0.3.5.jar --enrichment MSigDB:H --download databases/msigdb-h --source-file h.symbols.gmt --gene-id-type symbol
```

For an authorized academic KEGG user:

```shell
java -jar build/cli/jlinalg-0.3.5.jar --enrichment KEGG --download databases/kegg-human --species human --kegg-academic true
```

KEGG uses a rate-limited REST download; clusterProfiler support does not imply a
redistributable bundled KEGG database. No annotation database is shipped in the
JAR. Follow each provider's current terms and citation instructions:
[GO](https://geneontology.org/docs/download-go-annotations/),
[Reactome](https://reactome.org/download-data),
[HPO](https://hpo.jax.org/data/annotations),
[GWAS Catalog](https://www.ebi.ac.uk/gwas/downloads),
[Monarch](https://monarchinitiative.org/downloads),
[Open Targets](https://platform-docs.opentargets.org/licence),
[MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp),
[KEGG API](https://www.kegg.jp/kegg/rest/).

Downloaded species/namespace and Monarch evidence source are fixed by the
package. Install another package to change them. Gene-list-independent parsing
filters such as Catalog p threshold, Catalog gene source, GO evidence and Open
Targets score can be changed when analyzing the retained raw source; effective
settings are recorded. `--db-release` selects a provider release where supported.
Older GO archives predating the current GAF naming convention can be imported
with explicit source/ontology URLs. User source imports default to GMT for `.gmt`
or `.gmt.gz`; otherwise specify the actual `--db-format` if it is not `term2gene`.
For compressed source URLs, use URL paths ending in `.gz` so decompression can
be identified. URL redirects are followed; failed transfers have bounded retries.

For a custom disease-gene table, the required UTF-8 TSV/CSV columns are
`term_id` and `gene_id`, with optional `term_name`; each row is one membership.
GMT uses tab-separated term ID, description, then one gene per field. gzip is
accepted. Input tables support quoted cells, but not embedded physical newlines.

```shell
java -jar build/cli/jlinalg-0.3.5.jar --enrichment Custom --download databases/my-disease-sets --source-file disease-genes.tsv --db-format term2gene --gene-id-type entrez
```

Open Targets uses DuckDB JDBC to read local Parquet; its native libraries are
included in the executable dependency. Other enrichment modes do not load that
native reader. A supported native platform and writable temporary directory are
needed for Open Targets. On recent JDKs, `--enable-native-access=ALL-UNNAMED` before
`-jar` suppresses the native-access warning. No R installation is needed at runtime.

## Java API

```java
import java.util.List;
import java.util.Set;
import org.jlinalg.enrichment.EnrichmentAnalysis;
import org.jlinalg.enrichment.GeneSet;

var results = EnrichmentAnalysis.ora(
    Set.of("A", "B", "C", "D"), Set.of("A", "B"),
    List.of(new GeneSet("T1", "Example", Set.of("A", "C"))),
    1, 500, EnrichmentAnalysis.Fdr.BH);
System.out.println(results.get(0).pValue());
```

`EnrichmentAnalysis.gsameth` accepts the complete eligible feature-to-gene mapping,
selected feature IDs, gene sets, size limits, correction and Wallenius work cap.
`Gsameth.weights` exposes equivalent coverage, fractional selected contributions
and fitted probabilities for inspection. `Wallenius.upperTail` provides the
inclusive biased-urn tail independently of database loading. `adjustLog` exposes
BH/BY/NONE adjustment in original input order.

## Validation and limits

Regression tests cover hand-enumerated Fisher/Wallenius probabilities, independent
R `BiasedUrn` tails, the actual limma rank smoother, and an end-to-end missMethyl
`gsameth` fixture with multi-gene probes. The R comparator uses missMethyl commit
`c5d0a518aaacd98c9f31c4b6ca48239246363ad4`, limma 3.68.5 and BiasedUrn 2.0.12;
only the array-manifest adapter is replaced by explicit fixture mappings.
The test compares weighted overlap, bias odds, p-values and BH results. These
fixtures establish agreement for those cases, not universal finite-sample
calibration or a certified error bound for all biased-urn parameters.
An additional 90-digit Python Decimal calculation propagates the full count
distribution and sums its final tail. This checks relative accuracy down to
`2.4468354213672173e-106`, including a case where BiasedUrn's returned tail is zero.
Regenerate it with `python tools/generate_enrichment_precision.py`.

CLI tests cover analysis-derived eligibility, missing values, pre-mapping
Bonferroni counts, annotation joins, zero-hit families, namespace checks,
source formats, ontology closure/cycles, snapshot integrity and local Parquet.
Network access is not part of the automated test suite; live provider checks are
recorded below. Reproduce the offline checks with:

```shell
./gradlew check javadoc executableJar
```

To regenerate the R fixtures, install the pinned reference packages and run
`Rscript tools/generate_enrichment_fixtures.R` from the repository root. Inspect
the generated changes before accepting new reference versions.

### Recorded validation: 2026-09-13

The full `check javadoc` gate discovered 788 tests: 785 passed, three optional
CHOLMOD tests were skipped, and there were no failures or errors. `executableJar`
was built separately. The documented ORA and EWAS commands both completed.
The ORA example returns T1 p=0.05 and BH=0.10; the EWAS example returns weighted
overlap 2 and p approximately 0.99343 for T1. These tiny examples illustrate
the behavior and are not a realistic calibration study.

| Live downloaded collection | Snapshot | Usable sets before analysis-specific size filters |
| --- | --- | ---: |
| GO:BP, human symbols | GO pipeline 2026-08-05 | 13,798 |
| Reactome, human Ensembl | Current download, retained hashes | 2,839 |
| HPO, human Entrez | v2026-09-01 | 11,933 |
| GWAS Catalog, union | ZIP last modified 2026-09-04, retained hashes | 13,923 |
| Monarch:Mondo, human CURIEs | Current causal associations + Mondo v2026-09-01 | 9,132 |
| Open Targets, Ensembl | 26.06, 14 direct-association Parquet parts + disease labels | 8,658 |

Every listed source completed download, parsing, hash validation and package
installation. Open Targets was tested through the executable JAR on Windows;
both its newer `associationScore` and older `score` schemas have regression
coverage. A separate offline Catalog analysis confirmed the default union policy
in its screen output and log, testing 1,256 terms against a synthetic ten-gene
background. Counts depend on source releases and settings and are not expected
to remain constant. Live KEGG requests were not made without academic-user
authorization; its local parser is tested. MSigDB local GMT import, namespace
checking and corruption rejection are tested without a registered session.
