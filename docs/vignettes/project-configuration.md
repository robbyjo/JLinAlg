# Project configuration and reproducible follow-up

The source build supports local YAML configuration for every CLI command. No LLM,
model account, API key or external agent is needed. Build the executable with
`./gradlew check executableJar` (`.\gradlew.bat` on Windows).

## Precedence

Settings resolve in this order, from lowest to highest priority:

1. The command's built-in defaults.
2. User/machine settings in `~/.jlinalg/config.yaml` (or `--local-config FILE`).
3. `jlinalg.yaml` in the current project directory (or `--config FILE`).
4. Explicit CLI options.

The project configuration **overrides local configuration**. `--config` replaces
automatic project-file discovery, not the user defaults. Discovery checks the
current working directory only; it does not search ancestors. `--no-config`
disables both file layers and cannot be combined with explicit configuration paths.

## File schema

```yaml
schema_version: 1
commands:
  network:
    method: sparse
    matrix: {path: discovery.tsv}
    lambda: 0.15
    rule: and
    seed: 19
  variant-annotate:
    input: {path: variants.tsv}
    genome-build: GRCh38
```

The only root keys are `schema_version`, `defaults`, and `commands`. `defaults`
contains CLI options common to the commands you intend to run; prefer command
sections because commands reject inapplicable options. Use `association` for the
original formula CLI and `enrichment` for `--enrichment` commands. Other section
names equal their CLI subcommand, including `meta-analysis` and `variant-score`.

Option names omit the leading `--`. Scalar values become one CLI option; lists
become repeated options, for example `cohort: [A=cohort-a.tsv, B=cohort-b.tsv]`.
Maps merge by key; lists replace earlier lists. A null value removes an inherited
option. Boolean `help`, `overwrite`, `no-log`, `resume`, and `version` are switches;
other booleans become explicit `true`/`false` values. Quote strings that YAML 1.1
would interpret as booleans or numbers, such as `"on"` or `"001"`.

Use `{path: relative/file.tsv}` for a path resolved relative to the configuration
file defining it. Plain strings retain CLI semantics and are relative to the
working directory; `NAME=PATH` cohort strings are also passed unchanged. CLI paths
always remain relative to the working directory. Environment-variable expansion,
shell commands, YAML object tags, duplicate keys and aliases are not supported.
Store credentials outside these files; HTTP commands accept an environment-variable
**name** through `token-env`, never a literal token.

A machine file can supply R installation details without embedding them in a
shared project:

```yaml
schema_version: 1
commands:
  network:
    rscript: 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe'
    r-library: {path: r-library}
```

Use these settings only for R-backed network methods. Native methods reject
`rscript` and `r-library` instead of silently ignoring irrelevant configuration.

## Inspect and run the supplied project

From the repository root:

```shell
java -jar build/cli/jlinalg-0.3.6.jar config --command network --config examples/followup/jlinalg.yaml
java -jar build/cli/jlinalg-0.3.6.jar network --config examples/followup/jlinalg.yaml --lambda 0.2 --out build/followup/configured-sparse
```

The first command prints source paths/hashes and the configured argument expansion;
method-specific built-in defaults remain the command's responsibility and are
recorded in follow-up manifests. The second command uses the example matrix path
relative to its YAML file and overrides its lambda on the CLI.

## Provenance and output collisions

Configured runs with `--out` save `OUT.config.yaml`, including effective arguments
and the paths/SHA256 of configuration sources. New variant/network commands save
this sidecar even with `--no-config`. Existing unconfigured commands retain their
previous behavior. No sidecar is written after a failed run.

Variant/network `--out` values name **fresh directories**, not files. Results are
staged beside the destination, published on success, and accompanied by
`manifest.yaml` with input/output hashes, operation, settings and timestamps.
The configuration sidecar is adjacent to that directory. Failed runs print their
error and do not publish a partial result directory. For R adapters, successful
directories also retain the script, normalized inputs, external log and session
versions. Keep both the directory and the configuration sidecar with your paper.

Installed variant databases have a manifest that locks release, genome build,
source checksum and normalized-table checksum. Annotation verifies this checksum
on every use. Updating a database requires a new directory; there is no automatic
remote fallback or silent refresh. A general cross-workflow scheduler/lockfile
manager and LLM integration are not part of this addition.

## Related workflows

- [Variant annotation, consequences and scoring](variant-followup.md)
- [Networks and candidate regulators](network-followup.md)

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [Geir Kjetil Sandve et al. (2013) — Ten Simple Rules for Reproducible Computational Research](../CITATIONS.md#sandve-reproducibility-2013)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
