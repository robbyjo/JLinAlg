# Auditable result-first study example

This small project starts from two supplied association-result tables. It runs
input validation and standardization, a cross-cohort effect scatterplot with
four-decimal Pearson correlation, discovery/replication comparison and Venn
diagram, fixed-effect meta-analysis, deterministic diagnostics, and a traceable report.
The explicit cohort relationship avoids a pending default assumption. This tiny
offline teaching fixture uses `analysis.mode: manual` so it does not download
public annotation or enrichment resources.

From the repository root, first build the executable JAR:

```powershell
.\gradlew.bat executableJar
```

Copy this directory before running it if you want to keep the checked-in example
clean. In the copied directory, run:

```powershell
$jar = "D:\path\to\JLinAlg\build\cli\jlinalg-0.3.6.jar"
java -jar $jar study auto study.yaml --project . --actor-id scientist-1
```

The run creates `resolved-study.yaml`, `plans/`, `runs/`, the hash-chained
`ledger/events.jsonl`, derived `project-status.yaml`, `SESSION_HISTORY.md`,
`AGENTS.md`, searchable `index/index.html`, machine-readable
`index/evidence-summary.json`, and `index/LLM_CONTEXT.md`. Re-run the same command
to verify and reuse completed operation artifacts and retry unfinished work.
Each invocation keeps a separate immutable attempt; use `--rerun` to recompute
everything without deleting earlier evidence, or `--offline` to prohibit built-in
network requests. The index defaults to the latest attempt.

The run's operation directories include `effect-comparisons.tsv` and
`effect-scatterplots/*.png|pdf`, plus `replication-venn.tsv` and
`replication-venn/*.png|pdf`. In a hypothesis with several result datasets,
every pair is plotted. Reusing one cohort ID for two phenotype-specific result
datasets labels that pair `within_cohort`; distinct cohort IDs label it
`cross_cohort`.

When a supported rule needs scientist input, the run also creates `QUESTIONS.md`
and `state/questions.yaml`, with the exact questions and answers reflected in the
session notebook and searchable index. Inspect with `study questions --project .`;
export blank answers with `--out answers.yaml`. Check and record answers using
`study answer answers.yaml --project . --actor-id scientist-1 --check` or
`--record`. Supported factual answers produce proposed amendments requiring
separate review and approval. This minimal example may have no questions;
genetic follow-up or missing matrix metadata exercise the first rules. See the
[question workflow tutorial](../../docs/vignettes/study-pipeline-tutorial.md#documented-scientist-questions).

`data-dictionary.yaml` supplies the column roles and explanations for both result
files. Its contents and hash are snapshotted when the study is resolved. To
change a dictionary later, version it and approve a dataset amendment.

Optional API advice, literature review and manuscript drafts require explicit
configuration and permission; this fixture enables none of them. See the
[tutorial](../../docs/vignettes/study-pipeline-tutorial.md#data-dictionaries-optional-assistance-and-first-drafts)
for privacy settings, recorded interactions, and preapproved local sensitivities.

For a real end-to-end project, omit `mode: manual` or set `mode: auto`. Explicit
goals then become a minimum rather than a ceiling: JLinAlg adds eligible
meta-analysis, functional enrichment, genetic-evidence, literature, diagnostics,
and reporting branches. Public reference downloads default to automatic; an
unavailable scientific input is indexed as a blocker while independent branches
continue.

For `study auto`, exit 0 includes a successfully indexed partial run with input
requests; inspect its recorded state. Exit 2 indicates validation/execution
failure. `study run` additionally returns 3 for waiting operations unless
`--allow-partial` is supplied.
