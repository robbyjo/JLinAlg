# Plotting example

From the repository root, build the executable and render the example table:

```powershell
.\gradlew.bat executableJar
java -jar build/cli/jlinalg-0.3.6.jar plot `
  --type scatter --input examples/plotting/measurements.csv `
  --x age --y bmi --group sex --title "BMI by age" `
  --out build/plotting-example.pdf
```

Change the output extension to `png` or `svg` without changing the plot
definition. See [the complete plotting guide](../../docs/plotting.md).

Cluster rows and columns of a complete long-form heatmap:

```powershell
java -jar build/cli/jlinalg-0.3.6.jar plot `
  --type heatmap --input examples/plotting/heatmap.csv `
  --x column --y row --fill value --cluster both `
  --out build/clustered-heatmap.svg
```

Render the same weighted edge table as a network, chord diagram,
grouped radial edge bundle, or directed Sankey flow:

```powershell
java -jar build/cli/jlinalg-0.3.6.jar plot `
  --type sankey --input examples/plotting/flows.csv `
  --source source --target target --weight weight --group kind `
  --out build/participant-flow.pdf
```

`profiles.csv` is long-form input for a grouped `spider` plot.
`observations.csv` is wide-form input for `parallel` and `dendrogram` plots.
