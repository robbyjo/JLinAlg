# ACAT rare-variant tests

Use ACAT-V when a set may contain only a few associated variants. Use ACAT-O
when the genetic architecture is not known in advance: it Cauchy-combines
Burden, SKAT, and ACAT-V under Beta(1,25) and Beta(1,1) weighting. Both tests
operate on nuisance-adjusted scores and covariance, so the same APIs work with
one cohort, pooled independent cohorts, or an individual-level JLinAlg null
model.

## 1. Build the current executable

```powershell
./gradlew.bat executableJar --no-daemon --no-parallel
```

The examples below use the source-built
`build/cli/jlinalg-0.3.6.jar`. ACAT support is included in v0.3.6 and later.

## 2. Export one score file per cohort

The small checked-in example is quantitative and synthetic. Real analyses
should use a prespecified rare-variant mask and compatible cohort null models.

```powershell
java -jar build/cli/jlinalg-0.3.6.jar rare-score --vcf examples/rare-meta/cohort-a.vcf --pheno examples/rare-meta/cohort-a.tsv --id sample --response trait --genome-build GRCh38 --cov-window 100 --out build/rare-demo/a
java -jar build/cli/jlinalg-0.3.6.jar rare-score --vcf examples/rare-meta/cohort-b.vcf --pheno examples/rare-meta/cohort-b.tsv --id sample --response trait --genome-build GRCh38 --cov-window 100 --out build/rare-demo/b
```

The cohort manifest is `examples/rare-meta/cohorts.tsv`; paths in a manifest are
resolved relative to the manifest. Cohorts must be independent and use the same
trait scale, genome build, null-model contract, and allele convention.

## 3. Run ACAT-V and ACAT-O

```powershell
java -jar build/cli/jlinalg-0.3.6.jar rare-meta --cohorts examples/rare-meta/cohorts.tsv --genome-build GRCh38 --groups examples/rare-meta/groups.txt --test acat-v,acat-o --maf 0.5 --acat-mac-threshold 10 --out build/rare-demo/acat
```

The realistic default for `--maf` remains 0.05; this miniature example needs
0.5 to retain its deliberately common synthetic variants. ACAT-V defaults to
published Beta(1,25) weighting. `--weights equal|mb|beta|raremetal-beta`
selects a custom ACAT-V weighting, using the same values as burden coefficients
and Cauchy component weights as in the reference R custom-weight interface.
Canonical ACAT-O deliberately ignores that option and always reports
`weights=canonical`, because its six components are fixed by the paper.

Variants with pooled minor-allele count at or below
`--acat-mac-threshold` are collapsed into one burden component. JLinAlg derives
pooled MAC from the cohort allele counts used for the pooled MAF; it does not
approximate MAC by multiplying one global sample count. Variants above the
threshold contribute two-sided marginal score p-values. `acat-v.tsv` reports
the normalized Cauchy weight and p-value of every retained component.
`acat-o.tsv` reports all six set-test p-values.

ACAT-V is especially useful for sparse alternatives. ACAT-O trades a small
omnibus penalty for robustness across sparse, same-direction, and
mixed-direction alternatives. Neither result is a signed effect estimate; run
`--test burden,acat-v,acat-o` when an explicitly defined burden beta and SE are
also scientifically meaningful.

## 4. Use score-summary APIs

MAF and MAC arrays must be in exactly the same retained-variant order as the
score state, whose scores and covariance must use the corresponding minor-allele
orientation. When `ScoreMetaAnalysis.pool` filters columns, apply
`pooled.indices()` to both arrays first.

```java
import org.jlinalg.settest.*;

double[] scores = {1.8, -0.7, 2.1, 0.3};
double[] covariance = {
    4, .4, .2, 0,
    .4, 3, .3, .1,
    .2, .3, 2, .2,
    0, .1, .2, 1.5
};
double[] maf = {.0005, .002, .01, .03};
double[] mac = {2, 8, 40, 120};
var state = new SetTestScoreState(scores, covariance, scores.length);

var acatV = SummarySetTests.acatV("GENE1", state, maf, mac);
var acatVEqual = SummarySetTests.acatV(
    "GENE1", state, maf, mac, 1, 1, 10);
var acatO = SummarySetTests.acatO("GENE1", state, maf, mac);
var combined = Acat.combine(new double[]{.02, .0004, .2, .1, .8});
```

For participant-level workflows, `SetTests.acatV(...)` and
`SetTests.acatO(...)` prepare dosages, compute MAF/MAC, orient the score state to
the minor allele, and request one shared projection from any
`SetTestScoreNullModel`. The canonical ACAT-O
path then performs only the two eigendecompositions required by its two SKAT
components. `SetTests.acatVCustomWeights(...)` uses each `WeightedVariant`
weight as both the burden coefficient and Cauchy weight.

## 5. Accuracy and performance contract

The generic combiner evaluates very small p-values with the reciprocal Cauchy
tail and evaluates the final upper tail without subtracting from one. It rejects
invalid p-values, negative weights, all-zero weights, and an unresolved mixture
of exact zero and exact one endpoints.

`src/test/R/acat-reference.R` independently regenerates four score-summary
fixtures with R 4.6.1, the published ACAT/ACAT-V formulas, and CompQuadForm
Davies tails. ACAT and ACAT-V match those fixtures at tight floating-point
tolerance. ACAT-O inherits JLinAlg's existing positive-mixture tail; the tests
therefore state a separate bounded tolerance for the known Imhof-versus-Davies
SKAT component difference instead of attributing it to ACAT.

The reproducible command

```powershell
./gradlew.bat benchmarkAcat "-Djlinalg.benchmark.acat.variants=40" "-Djlinalg.benchmark.acat.batch=3" "-Djlinalg.benchmark.acat.warmups=2" "-Djlinalg.benchmark.acat.measurements=5" --no-daemon --no-parallel
```

measured 6.734 ms per fused ACAT-O versus 8.555 ms for composing the six public
tests separately on this development host: 1.270x faster, with p-value
difference `2.78e-16`. Preparation and data matching are outside both timed
regions. These host- and workload-specific medians are not a universal speed
guarantee.

The current CLI contract is quantitative Gaussian score meta-analysis. ACAT
does not repair conservative or otherwise miscalibrated component p-values;
binary rare-case tails, overlapping cohorts, and phenotype-model mismatch must
be handled before combination. The Cauchy calibration is the published
dependence-robust approximation, not an exact finite-sample permutation test.

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [Yaowu Liu et al. (2019) — ACAT: A fast and powerful p value combination method for rare-variant analysis in sequencing studies](../CITATIONS.md#liu-acat-2019) — [PMID: 30849328](https://pubmed.ncbi.nlm.nih.gov/30849328/) · [PMCID: PMC6407498](https://pmc.ncbi.nlm.nih.gov/articles/PMC6407498/)
- [Michael C. Wu et al. (2011) — Rare-variant association testing for sequencing data with the sequence kernel association test](../CITATIONS.md#wu-skat-2011) — [PMID: 21737059](https://pubmed.ncbi.nlm.nih.gov/21737059/) · [PMCID: PMC3135811](https://pmc.ncbi.nlm.nih.gov/articles/PMC3135811/)
- [Seunggeun Lee et al. (2012) — Optimal unified approach for rare-variant association testing](../CITATIONS.md#lee-skato-2012) — [PMID: 22863193](https://pubmed.ncbi.nlm.nih.gov/22863193/) · [PMCID: PMC3415556](https://pmc.ncbi.nlm.nih.gov/articles/PMC3415556/)
- [B. E. Madsen and S. R. Browning (2009) — A groupwise association test for rare mutations using a weighted sum statistic](../CITATIONS.md#madsen-browning-2009) — [PMID: 19214210](https://pubmed.ncbi.nlm.nih.gov/19214210/) · [PMCID: PMC2633048](https://pmc.ncbi.nlm.nih.gov/articles/PMC2633048/)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
