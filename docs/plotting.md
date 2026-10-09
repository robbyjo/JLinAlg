# Command-line plotting

JLinAlg can render publication-ready plots directly from CSV or TSV tables.
The output filename selects a raster format (`png`, `jpg`, `bmp`, `gif`, or
`tiff`) or a vector format (`svg`, `eps`, or `pdf`). Existing files are never
silently replaced.

## Quick start

Given [`examples/plotting/measurements.csv`](../examples/plotting/measurements.csv):

```text
age,bmi,sex,weight,visit
41,24.2,F,58,baseline
53,31.7,M,91,baseline
47,27.5,F,68,follow-up
```

create a grouped scatter plot:

```powershell
java -jar build/cli/jlinalg-0.3.6.jar plot `
  --type scatter --input examples/plotting/measurements.csv `
  --x age --y bmi --group sex `
  --title "BMI by age" --out bmi-by-age.png
```

The same command can write a vector PDF or SVG by changing `--out`:

```powershell
java -jar build/cli/jlinalg-0.3.6.jar plot `
  --type violin --input examples/plotting/measurements.csv `
  --x sex --y bmi --out bmi-by-sex.pdf
```

## Chart types and mappings

| `--type` | Required mappings | Notes |
|---|---|---|
| `scatter` | `--x`, `--y` | Optional `--group` creates colored series and a legend. |
| `line`, `step`, `area` | `--x`, `--y` | Uses input row order within each optional group. |
| `connected-scatter` | `--x`, `--y` | Draws both points and connecting lines. |
| `bubble` | `--x`, `--y`, `--size` | Sizes must be positive; grouping is supported. |
| `histogram` | `--x` | `--bins` defaults to 30. |
| `density` | `--x` | Gaussian kernel density with an automatic bandwidth. |
| `dot` | `--x`, `--y` | Repeated category rows are represented by their mean. |
| `box`, `violin` | `--x GROUP`, `--y VALUE` | `--group` can be used instead of `--x`; the violin also marks its IQR and median. |
| `bar` | `--x`, optional `--y` | Without `--y`, rows are counted. Repeated rows are summed. Group with `--position grouped` or `stacked`. |
| `pie`, `donut` | `--x`, optional `--y` | Without `--y`, rows are counted; repeated slice labels are summed. |
| `heatmap` | `--x COLUMN`, `--y ROW`, `--fill VALUE` | The long-form grid must be complete. `--cluster rows`, `columns`, or `both` adds average-linkage Euclidean dendrograms and reorders the cells. `vertical` and `horizontal` are accepted aliases for rows and columns. |
| `spider`, `radar` | `--x AXIS`, `--y VALUE` | Values must already be scaled to `[0,1]`; optional `--group` draws multiple polygons. |
| `circular-bar` | `--x CATEGORY`, `--y VALUE` | Optional `--group` controls bar colors and the legend. |
| `parallel` | `--id LABEL`, `--columns A,B,C` | Each numeric dimension is independently scaled; optional `--group` colors observations. |
| `dendrogram` | `--id LABEL`, `--columns A,B,C` | Average-linkage Euclidean clustering; `--orientation vertical` or `horizontal`. |
| `treemap` | `--x LABEL`, `--y VALUE` | Positive areas; optional `--group` adds a first hierarchy level. |
| `network` | `--source`, `--target` | Deterministic force-directed layout; optional positive `--weight` and edge `--group`. |
| `chord` | `--source`, `--target` | Weighted radial chord links and colored node arcs. |
| `edge-bundling` | `--source`, `--target` | Radial links share group-specific bundle control points; `--group` chooses bundles. |
| `sankey` | `--source`, `--target` | Positive `--weight` controls flow width. Input must be a directed acyclic graph. |
| `venn` | `--x ITEM`, `--group SET` | Exactly two set labels are required. Repeat an item in both sets to place it in the overlap. Circle area is illustrative; the three printed region counts are exact. |

Missing values (`NA`, `N/A`, `null`, `.`, and blank fields) are omitted with a
reported row count. Use `--na error` to make any incomplete mapped row fail the
command. Non-numeric text in a numeric mapping always fails.

## Axes, legends, and appearance

Axis titles default to their mapped column names. Override them with
`--x-label` and `--y-label`. The common controls are:

```text
--title TEXT
--width 800 --height 600
--theme default|ggplot2|minimal|classic
--legend auto|none|right|bottom|inside-ne|inside-nw|inside-se|inside-sw
--no-grid
--log-x --log-y
--x-range MIN,MAX --y-range MIN,MAX
--number-pattern 0.###
--point-size 8
--point-shape circle|square|diamond|triangle-up|triangle-down|cross|plus|none
--line-width 2.0
--line-pattern solid|dashed|dotted|dash-dot
--colors '#0072B2,#D55E00,#009E73'
```

`--marker-size` and `--dot-size` are aliases for `--point-size`;
`--dot-pattern` is an alias for `--point-shape`. Palettes are applied in the
listed order to groups or series. For a heatmap, the palette is interpreted as
a continuous multi-stop gradient; its default is blue–white–red.

Use quotes around a hexadecimal palette in PowerShell because `#` otherwise
starts a comment.

## Study-pipeline effect and replication figures

`study auto` generates `effect-comparisons.tsv` plus PNG and vector PDF
scatterplots for every unordered pair of result tables in a hypothesis. Feature
IDs are matched exactly. Each plot uses common x/y limits, includes the `y = x`
identity line, and shows `Pearson r` to exactly four decimals in its legend.
Different declared cohort IDs produce a `cross_cohort` comparison; two result
tables with the same cohort ID produce a `within_cohort` comparison, which is
the intended contract for correlated phenotypes analyzed in one cohort.

Discovery/replication operations additionally write `replication-venn.tsv` and
one PNG/PDF Venn diagram per discovery/replication pair. The left set is the
prespecified discovery-selected set. The right set contains features satisfying
the approved replication P-value and direction criterion. The overlap is the
replicated set. The table records the multiple-testing scope and same-sign rule,
so a `discovery_subset` analysis is not mistaken for a symmetric genome-wide
comparison.

A standalone two-set diagram uses the same membership contract:

```powershell
java -jar build/cli/jlinalg-0.3.6.jar plot `
  --type venn --input examples/plotting/venn-memberships.csv `
  --x feature --group set --title "Discovery and replication" `
  --out discovery-replication.pdf
```

## Scope

The plotting layer covers the statistical, matrix, hierarchy, network, chord,
flow, and circular charts above, multi-series axes and legends, themes, and
consistent raster/vector export. The R Graph Gallery also catalogs geospatial,
3D, animated, interactive, word-cloud, cartogram, and data-art visualizations.
Those still require specialized data or runtime models and are not represented
by this command yet.
