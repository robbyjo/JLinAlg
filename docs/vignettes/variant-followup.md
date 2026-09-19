# Variant annotation, consequences and evidence scoring

These source-build commands implement three distinct operations: annotation
lookup, consequence prediction/import, and transparent evidence integration.
They require no LLM. All `--out` arguments name fresh directories. Examples are
synthetic; the coordinates and annotations do not assert actual biology.

## Inputs and variant identity

Canonical TSV/CSV inputs contain `genome_build`, `chrom`, `pos`, `ref`, `alt`.
Use GRCh37 or GRCh38 and positive, one-based **VCF** coordinates. Alleles must be
distinct biallelic A/C/G/T sequences. Symbolic alleles and unsplit multiallelic
records are rejected. The optional `variant_id` must equal
`GRCh38:1:101:A:G`, including build and alleles. A lowercase `chr` prefix is
removed and alleles are uppercased for identity matching. rsIDs alone are not keys.

Normalize indels against the chosen reference before importing, for example with
`bcftools norm -f reference.fa -m -any`. This workflow does not left-align, liftover,
infer strand, validate a REF allele against FASTA during local lookup, or repair
incompatible representations. Different representations do not silently match.
Keep transcript-level records: several annotations for one variant are legitimate.

## 1. Install and query a local annotation snapshot

```shell
java -jar build/cli/jlinalg-0.3.6.jar variant-db --source-file examples/followup/annotations.tsv --genome-build GRCh38 --release synthetic-1 --license CC0 --out build/followup/annotation-db
java -jar build/cli/jlinalg-0.3.6.jar variant-annotate --input examples/followup/variants.tsv --genome-build GRCh38 --database build/followup/annotation-db --out build/followup/annotated
```

Inspect `annotations.tsv`: the first example variant retains both transcripts;
the second has `not_found`. Missing annotation is not a benign classification.
Annotation columns are prefixed with `annotation.`. The join does not collapse
genes/transcripts or automatically merge observational association statistics.

`variant-db --source-url URL --source-sha256 SHA256` downloads instead of copying
a local source. Supply the provider's exact release URL and independently obtained
SHA256, `--release`, `--license`, and `--genome-build`. `--compression gzip`
decompresses gzip/BGZF sources. A canonical source is a TSV with the five key
columns and any additional annotation/score columns.

For other TSV schemas, use `--format table --chrom-column '#chr' --pos-column
'pos(1-based)' --ref-column ref --alt-column alt` with names from the actual file.
This supports precomputed-score tables such as a suitably selected dbNSFP release;
check which assembly its position column represents. Original fields are retained
as `source.*`. `--format vcf` imports a normalized, biallelic VCF and retains its
ID and INFO string. It does not interpret INFO arrays as independent scores.
Restricted databases are user-supplied; provider licensing and access conditions
still apply, including commercial-use differences among score collections.

The current backend validates and **streams the entire local TSV**, retaining
only query matches. It is not a genome-scale index: use a regional/preselected
snapshot for frequent small queries against very large resources. Installation
requires disk space for both the download and normalized table. Checksums are
verified before analysis, adding another sequential read.

## Remote POST and provider APIs

For a controlled service that implements the canonical TSV contract:

```shell
java -jar build/cli/jlinalg-0.3.6.jar variant-annotate --input examples/followup/variants.tsv --genome-build GRCh38 --backend post --endpoint https://YOUR-SERVICE/annotations --out build/followup/remote
```

Replace the example URL with an actual service. JLinAlg posts **only the five
canonical variant columns** as `text/tab-separated-values`, expects a canonical
annotation TSV, and stores request/response files. This is a documented service
contract, not a claim that arbitrary web forms accept TSV. There is no browser
automation. HTTP errors, redirects, malformed rows and build mismatches fail the
run. Use a final URL, `--timeout`, `--max-download-mb`, and, if needed,
`--token-env NAME` for bearer authentication. No implicit switch to remote occurs.

The `vep-rest` consequence backend implements the actual Ensembl JSON POST API,
in batches of at most 200 variants. It chooses the GRCh37 endpoint for GRCh37,
checks the returned assembly and input identities, and retains raw JSON. It
requires network access and sends coordinates/alleles to the selected service:

```shell
java -jar build/cli/jlinalg-0.3.6.jar variant-consequence --input REAL_VARIANTS.tsv --genome-build GRCh38 --engine vep-rest --out build/followup/vep-rest
```

The hosted service can change releases. Archived response hashes make the
performed lookup auditable; rerunning against a live endpoint need not reproduce
an older release. Use local VEP/cache versions when that matters.

## 2. Consequence prediction or import

Import the runnable synthetic VEP result:

```shell
java -jar build/cli/jlinalg-0.3.6.jar variant-consequence --input examples/followup/variants.tsv --genome-build GRCh38 --engine import --annotations examples/followup/vep.tsv --format vep --out build/followup/consequences
```

Supported imports are canonical TSV, VEP tab output beginning with
`#Uploaded_variation`, and ANNOVAR's `*_multianno.vcf` (`--format annovar-vcf`).
VEP Uploaded_variation / VCF ID must preserve the canonical variant ID supplied
by JLinAlg. Imported files must have been generated for the declared build; a
VEP text file without explicit assembly metadata cannot independently prove it.
The ANNOVAR import additionally verifies coordinates and REF/ALT. VEP transcripts
are retained; ANNOVAR aggregated gene fields/INFO remain aggregated and the
transcript is `NA` when it is not explicitly available.

Run an installed VEP with a compatible local cache and FASTA:

```shell
java -jar build/cli/jlinalg-0.3.6.jar variant-consequence --input REAL_VARIANTS.tsv --genome-build GRCh38 --engine vep --executable /path/to/vep --cache /path/to/cache --cache-version VERSION --fasta /path/to/GRCh38.fa --out build/followup/vep-local
```

JLinAlg writes VCF IDs, invokes VEP with `--offline --cache --tab --everything`,
and imports the resulting transcript records. Choose an executable matching your
platform; an external Perl environment is still required by VEP.

Run a separately installed ANNOVAR:

```shell
java -jar build/cli/jlinalg-0.3.6.jar variant-consequence --input REAL_VARIANTS.tsv --genome-build GRCh38 --engine annovar --executable /path/to/table_annovar.pl --perl /path/to/perl --database /path/to/humandb --protocol refGene --operation g --out build/followup/annovar
```

The adapter uses `-vcfinput`, hg38/hg19 as appropriate, and reads the resulting
multianno VCF. Protocols/operations are explicit comma-separated ANNOVAR options.
Install compatible databases using ANNOVAR's downloader. JLinAlg does not bundle
or automatically register for ANNOVAR. External commands use argument arrays,
not a shell. Failed processes/timeouts do not publish partial results.

## 3. Integrate evidence without disguising uncertainty

```shell
java -jar build/cli/jlinalg-0.3.6.jar variant-score --input examples/followup/evidence.tsv --components examples/followup/components.tsv --out build/followup/scores
```

`components.tsv` specifies `column`, positive `weight`, `minimum`, `maximum`, and
`direction` (`higher` or `lower`). A component is linearly scaled to [0,1] using
those **predeclared** bounds, reversed for `lower`, and multiplied by its normalized
weight. `priority_score` is the sum of contributions. Values outside the bounds
are rejected rather than silently clipped. All original fields and component
contributions remain in the output; rows retain input order, with no transcript
or gene aggregation.

`NA`, empty or `.` component values produce a missing contribution and an `NA`
total score. `weight_coverage` reports available weight; there is no automatic
reweighting that rewards poorly annotated variants. Other nonfinite values fail.

This is a user-defined prioritization index, **not** a pathogenicity/causality
probability or a published combined test. PIP and coloc evidence may overlap;
weights do not make them independent. Avoid counting multiple correlated
functional predictors as separate corroboration. Select candidates with the
underlying evidence, tissue context and sensitivity to weights visible.

## Validation and supported scope

Automated tests cover YAML precedence, checksum failure, duplicate/mismatched
variants, transcript-preserving lookup, missing evidence, score arithmetic,
local HTTP POST/VEP JSON fixtures, VEP/ANNOVAR imports and atomic failure.
The live public VEP smoke check on 2026-09-19 was unavailable (direct HTTPS POST
returned HTTP 503; Java's connection ended before a response). Live-service
success is therefore unverified on this host; local protocol fixtures passed.
Installed VEP/ANNOVAR execution requires the user's tool and database setup.
Their local full-cache execution is not established by the import/HTTP tests.

See [configuration](project-configuration.md) and [network follow-up](network-followup.md).

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [Kai Wang, Mingyao Li, and Hakon Hakonarson (2010) — ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data](../CITATIONS.md#wang-annovar-2010)
- [William McLaren et al. (2016) — The Ensembl Variant Effect Predictor](../CITATIONS.md#mclaren-vep-2016)
- [Geir Kjetil Sandve et al. (2013) — Ten Simple Rules for Reproducible Computational Research](../CITATIONS.md#sandve-reproducibility-2013)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
