/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.awt.Color;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import org.jlinalg.plot.DiagramPlots;
import org.jlinalg.plot.LegendPosition;
import org.jlinalg.plot.LinePattern;
import org.jlinalg.plot.Plot;
import org.jlinalg.plot.PlotOptions;
import org.jlinalg.plot.PlotTheme;
import org.jlinalg.plot.PointShape;
import org.jlinalg.plot.Plots;
import org.jlinalg.plot.XYPlotBuilder;

/** General CSV/TSV-to-plot command line. */
final class PlotCli {
    private PlotCli() { }

    static int run(String[] arguments, PrintStream output, PrintStream error)
            throws IOException {
        Options options = Options.parse(arguments);
        if (options.help) {
            output.print(help());
            return 0;
        }
        DelimitedData table = DelimitedData.read(options.input);
        PipelinePaths.requireFreshOutputs(options.output);
        Counter omitted = new Counter();
        Plot plot = create(table, options, omitted);
        plot.save(options.output);
        output.println("Wrote " + options.type + " plot to " + options.output);
        if (omitted.value > 0)
            error.println("jlinalg plot: omitted " + omitted.value
                + " row(s) with missing required values");
        return 0;
    }

    private static Plot create(DelimitedData table, Options options,
            Counter omitted) {
        return switch (options.type) {
            case "scatter", "line", "connected-scatter", "area", "step" ->
                xy(table, options, omitted);
            case "bubble" -> bubble(table, options, omitted);
            case "histogram" -> Plots.histogram(numbers(table, require(options.x,
                "--x"), options, omitted), options.bins,
                options.plotOptions(options.x, "Count"));
            case "density" -> Plots.density(numbers(table, require(options.x,
                "--x"), options, omitted), options.plotOptions(options.x, "Density"));
            case "dot" -> dot(table, options, omitted);
            case "bar" -> bar(table, options, omitted);
            case "pie", "donut" -> pie(table, options, omitted);
            case "box", "violin" -> distribution(table, options, omitted);
            case "heatmap" -> heatmap(table, options, omitted);
            case "spider", "radar" -> spider(table, options, omitted);
            case "circular-bar" -> circularBar(table, options, omitted);
            case "parallel" -> parallel(table, options, omitted);
            case "dendrogram" -> dendrogram(table, options, omitted);
            case "treemap" -> treemap(table, options, omitted);
            case "network" -> DiagramPlots.network(edges(table, options, omitted),
                options.plotOptions("", ""));
            case "chord" -> DiagramPlots.chord(edges(table, options, omitted),
                options.plotOptions("", ""));
            case "edge-bundling" -> DiagramPlots.edgeBundle(
                edges(table, options, omitted), options.plotOptions("", ""));
            case "sankey" -> DiagramPlots.sankey(edges(table, options, omitted),
                options.plotOptions("", ""));
            case "venn" -> venn(table, options, omitted);
            default -> throw new IllegalArgumentException("unsupported plot type: "
                + options.type);
        };
    }

    private static Plot xy(DelimitedData table, Options options, Counter omitted) {
        String xName = require(options.x, "--x");
        String yName = require(options.y, "--y");
        Map<String, List<Point>> groups = points(table, options, false, omitted);
        XYPlotBuilder builder = Plots.xy().options(options.plotOptions(xName, yName));
        groups.forEach((name, values) -> {
            double[] x = values.stream().mapToDouble(Point::x).toArray();
            double[] y = values.stream().mapToDouble(Point::y).toArray();
            switch (options.type) {
                case "scatter" -> builder.points(name, x, y);
                case "line" -> builder.line(name, x, y);
                case "connected-scatter" -> builder.lineAndPoints(name, x, y);
                case "area" -> builder.area(name, x, y);
                case "step" -> builder.step(name, x, y);
                default -> throw new IllegalStateException("unexpected XY type");
            }
        });
        return builder.build();
    }

    private static Plot bubble(DelimitedData table, Options options,
            Counter omitted) {
        String xName = require(options.x, "--x");
        String yName = require(options.y, "--y");
        require(options.size, "--size");
        Map<String, List<Point>> groups = points(table, options, true, omitted);
        Map<String, Plots.BubbleData> data = new LinkedHashMap<>();
        groups.forEach((name, values) -> data.put(name, new Plots.BubbleData(
            values.stream().mapToDouble(Point::x).toArray(),
            values.stream().mapToDouble(Point::y).toArray(),
            values.stream().mapToDouble(Point::size).toArray())));
        return Plots.bubble(data, options.plotOptions(xName, yName));
    }

    private static Plot dot(DelimitedData table, Options options, Counter omitted) {
        String category = require(options.x, "--x");
        String value = require(options.y, "--y");
        Aggregation aggregation = aggregate(table, category, value, null,
            false, options, omitted);
        double[] means = aggregation.singleValues(true);
        return Plots.dotPlot(aggregation.categories(), means,
            options.plotOptions(category, value));
    }

    private static Plot bar(DelimitedData table, Options options, Counter omitted) {
        String category = require(options.x, "--x");
        Aggregation aggregation = aggregate(table, category, options.y,
            options.group, options.y == null, options, omitted);
        PlotOptions plotOptions = options.plotOptions(category,
            options.y == null ? "Count" : options.y);
        if (aggregation.groups().size() == 1)
            return Plots.barChart(aggregation.categories(),
                aggregation.singleValues(false), plotOptions);
        Map<String, double[]> series = aggregation.series(false);
        return options.position.equals("stacked")
            ? Plots.stackedBarChart(aggregation.categories(), series, plotOptions)
            : Plots.groupedBarChart(aggregation.categories(), series, plotOptions);
    }

    private static Plot pie(DelimitedData table, Options options, Counter omitted) {
        String category = require(options.x, "--x");
        Aggregation aggregation = aggregate(table, category, options.y, null,
            options.y == null, options, omitted);
        double[] values = aggregation.singleValues(false);
        Map<String, Double> slices = new LinkedHashMap<>();
        for (int index = 0; index < values.length; index++)
            slices.put(aggregation.categories().get(index), values[index]);
        PlotOptions plotOptions = options.plotOptions("", "");
        return options.type.equals("donut") ? Plots.donutChart(slices, plotOptions)
            : Plots.pieChart(slices, plotOptions);
    }

    private static Plot distribution(DelimitedData table, Options options,
            Counter omitted) {
        String value = require(options.y, "--y");
        String group = options.group != null ? options.group : options.x;
        int valueIndex = table.column(value);
        int groupIndex = group == null ? -1 : table.column(group);
        Map<String, List<Double>> selected = new LinkedHashMap<>();
        for (String[] row : table.rows()) {
            Double number = number(row[valueIndex], value);
            String label = groupIndex < 0 ? "all" : text(row[groupIndex]);
            if (number == null || label == null) {
                missing(options, omitted, group == null ? value : value + ", " + group);
                continue;
            }
            selected.computeIfAbsent(label, ignored -> new ArrayList<>()).add(number);
        }
        if (selected.isEmpty()) throw new IllegalArgumentException("plot has no complete rows");
        Map<String, double[]> groups = new LinkedHashMap<>();
        selected.forEach((label, values) -> groups.put(label,
            values.stream().mapToDouble(Double::doubleValue).toArray()));
        PlotOptions plotOptions = options.plotOptions(group == null ? "" : group, value);
        return options.type.equals("box") ? Plots.boxPlot(groups, plotOptions)
            : Plots.violinPlot(groups, plotOptions);
    }

    private static Plot heatmap(DelimitedData table, Options options,
            Counter omitted) {
        String columnName = require(options.x, "--x");
        String rowName = require(options.y, "--y");
        String fillName = require(options.fill, "--fill");
        int columnIndex = table.column(columnName);
        int rowIndex = table.column(rowName);
        int fillIndex = table.column(fillName);
        Map<String, Integer> columns = new LinkedHashMap<>();
        Map<String, Integer> rows = new LinkedHashMap<>();
        List<Cell> cells = new ArrayList<>();
        for (String[] row : table.rows()) {
            String column = text(row[columnIndex]);
            String rowLabel = text(row[rowIndex]);
            Double value = number(row[fillIndex], fillName);
            if (column == null || rowLabel == null || value == null) {
                missing(options, omitted, "heatmap mappings");
                continue;
            }
            columns.computeIfAbsent(column, ignored -> columns.size());
            rows.computeIfAbsent(rowLabel, ignored -> rows.size());
            cells.add(new Cell(rowLabel, column, value));
        }
        if (cells.isEmpty()) throw new IllegalArgumentException("heatmap has no complete cells");
        double[][] values = new double[rows.size()][columns.size()];
        for (double[] row : values) Arrays.fill(row, Double.NaN);
        for (Cell cell : cells) {
            int row = rows.get(cell.row()), column = columns.get(cell.column());
            if (Double.isFinite(values[row][column]))
                throw new IllegalArgumentException("duplicate heatmap cell: "
                    + cell.row() + ", " + cell.column());
            values[row][column] = cell.value();
        }
        for (int row = 0; row < values.length; row++)
            for (int column = 0; column < values[row].length; column++)
                if (!Double.isFinite(values[row][column]))
                    throw new IllegalArgumentException("heatmap grid is incomplete at "
                        + new ArrayList<>(rows.keySet()).get(row) + ", "
                        + new ArrayList<>(columns.keySet()).get(column));
        Color[] fallback = {new Color(0x21, 0x66, 0xAC), Color.WHITE,
            new Color(0xB2, 0x18, 0x2B)};
        return DiagramPlots.heatmap(new DiagramPlots.HeatmapData(
            List.copyOf(columns.keySet()), List.copyOf(rows.keySet()), values),
            options.clusterMode(), options.plotOptions(columnName, rowName, fallback));
    }

    private static Plot spider(DelimitedData table, Options options,
            Counter omitted) {
        String axisName = require(options.x, "--x");
        String valueName = require(options.y, "--y");
        int axisIndex = table.column(axisName), valueIndex = table.column(valueName);
        int groupIndex = options.group == null ? -1 : table.column(options.group);
        Map<String, Integer> axes = new LinkedHashMap<>();
        Map<String, Map<String, Double>> selected = new LinkedHashMap<>();
        for (String[] row : table.rows()) {
            String axis = text(row[axisIndex]);
            String group = groupIndex < 0 ? "data" : text(row[groupIndex]);
            Double value = number(row[valueIndex], valueName);
            if (axis == null || group == null || value == null) {
                missing(options, omitted, "spider mappings");
                continue;
            }
            axes.computeIfAbsent(axis, ignored -> axes.size());
            if (selected.computeIfAbsent(group, ignored -> new LinkedHashMap<>())
                    .put(axis, value) != null)
                throw new IllegalArgumentException("duplicate spider axis in series: "
                    + group + ", " + axis);
        }
        Map<String, double[]> series = new LinkedHashMap<>();
        selected.forEach((group, mapped) -> {
            double[] values = new double[axes.size()];
            int index = 0;
            for (String axis : axes.keySet()) {
                Double value = mapped.get(axis);
                if (value == null) throw new IllegalArgumentException(
                    "spider series is missing axis " + axis + ": " + group);
                values[index++] = value;
            }
            series.put(group, values);
        });
        return DiagramPlots.spider(List.copyOf(axes.keySet()), series,
            options.plotOptions(axisName, valueName));
    }

    private static Plot circularBar(DelimitedData table, Options options,
            Counter omitted) {
        String category = require(options.x, "--x");
        String value = require(options.y, "--y");
        Aggregation aggregation = aggregate(table, category, value, options.group,
            false, options, omitted);
        List<DiagramPlots.TreeItem> bars = new ArrayList<>();
        aggregation.groups().forEach((group, summaries) -> {
            for (String label : aggregation.categories()) {
                Summary summary = summaries.get(label);
                if (summary != null)
                    bars.add(new DiagramPlots.TreeItem(label, group, summary.sum()));
            }
        });
        return DiagramPlots.circularBar(bars, options.plotOptions(category, value));
    }

    private static Plot parallel(DelimitedData table, Options options,
            Counter omitted) {
        DiagramPlots.ObservationData data = observations(table, options, omitted);
        return DiagramPlots.parallel(data, options.plotOptions("", ""));
    }

    private static Plot dendrogram(DelimitedData table, Options options,
            Counter omitted) {
        DiagramPlots.ObservationData data = observations(table, options, omitted);
        return DiagramPlots.dendrogram(data, options.orientation(),
            options.plotOptions("", ""));
    }

    private static Plot treemap(DelimitedData table, Options options,
            Counter omitted) {
        String category = require(options.x, "--x");
        String value = require(options.y, "--y");
        Aggregation aggregation = aggregate(table, category, value, options.group,
            false, options, omitted);
        List<DiagramPlots.TreeItem> items = new ArrayList<>();
        aggregation.groups().forEach((group, summaries) -> {
            for (String label : aggregation.categories()) {
                Summary summary = summaries.get(label);
                if (summary != null)
                    items.add(new DiagramPlots.TreeItem(label, group, summary.sum()));
            }
        });
        return DiagramPlots.treemap(items, options.plotOptions("", ""));
    }

    private static DiagramPlots.ObservationData observations(DelimitedData table,
            Options options, Counter omitted) {
        String idName = require(options.id, "--id");
        if (options.columns == null || options.columns.isEmpty())
            throw new IllegalArgumentException("--columns is required");
        int idIndex = table.column(idName);
        int groupIndex = options.group == null ? -1 : table.column(options.group);
        int[] columnIndices = options.columns.stream().mapToInt(table::column).toArray();
        List<String> labels = new ArrayList<>(), groups = new ArrayList<>();
        List<double[]> values = new ArrayList<>();
        Set<String> unique = new LinkedHashSet<>();
        for (String[] row : table.rows()) {
            String label = text(row[idIndex]);
            String group = groupIndex < 0 ? "data" : text(row[groupIndex]);
            double[] numbers = new double[columnIndices.length];
            boolean complete = label != null && group != null;
            for (int column = 0; column < columnIndices.length; column++) {
                Double value = number(row[columnIndices[column]], options.columns.get(column));
                if (value == null) complete = false; else numbers[column] = value;
            }
            if (!complete) {
                missing(options, omitted, "observation mappings");
                continue;
            }
            if (!unique.add(label))
                throw new IllegalArgumentException("duplicate observation ID: " + label);
            labels.add(label);
            groups.add(group);
            values.add(numbers);
        }
        if (values.isEmpty()) throw new IllegalArgumentException(
            "plot has no complete observations");
        return new DiagramPlots.ObservationData(labels, options.columns,
            values.toArray(double[][]::new), groups);
    }

    private static List<DiagramPlots.Edge> edges(DelimitedData table,
            Options options, Counter omitted) {
        String sourceName = require(options.source, "--source");
        String targetName = require(options.target, "--target");
        int sourceIndex = table.column(sourceName), targetIndex = table.column(targetName);
        int weightIndex = options.weight == null ? -1 : table.column(options.weight);
        int groupIndex = options.group == null ? -1 : table.column(options.group);
        List<DiagramPlots.Edge> result = new ArrayList<>();
        for (String[] row : table.rows()) {
            String source = text(row[sourceIndex]), target = text(row[targetIndex]);
            String group = groupIndex < 0 ? "edges" : text(row[groupIndex]);
            Double weight = weightIndex < 0 ? 1.0 : number(row[weightIndex], options.weight);
            if (source == null || target == null || group == null || weight == null) {
                missing(options, omitted, "edge mappings");
                continue;
            }
            result.add(new DiagramPlots.Edge(source, target, weight, group));
        }
        if (result.isEmpty()) throw new IllegalArgumentException("plot has no complete edges");
        return result;
    }

    private static Plot venn(DelimitedData table, Options options,
            Counter omitted) {
        String itemName = require(options.x, "--x");
        String setName = require(options.group, "--group");
        int itemIndex = table.column(itemName), setIndex = table.column(setName);
        Map<String, Set<String>> sets = new LinkedHashMap<>();
        for (String[] row : table.rows()) {
            String item = text(row[itemIndex]), set = text(row[setIndex]);
            if (item == null || set == null) {
                missing(options, omitted, "Venn mappings");
                continue;
            }
            sets.computeIfAbsent(set, ignored -> new LinkedHashSet<>()).add(item);
        }
        if (sets.size() != 2)
            throw new IllegalArgumentException("Venn plot requires exactly two sets in --group");
        var iterator = sets.entrySet().iterator();
        var left = iterator.next();
        var right = iterator.next();
        Set<String> common = new LinkedHashSet<>(left.getValue());
        common.retainAll(right.getValue());
        return DiagramPlots.venn(new DiagramPlots.VennData(left.getKey(),
            right.getKey(), left.getValue().size() - common.size(), common.size(),
            right.getValue().size() - common.size()),
            options.plotOptions("", ""));
    }

    private static Map<String, List<Point>> points(DelimitedData table,
            Options options, boolean withSize, Counter omitted) {
        int x = table.column(require(options.x, "--x"));
        int y = table.column(require(options.y, "--y"));
        int size = withSize ? table.column(require(options.size, "--size")) : -1;
        int group = options.group == null ? -1 : table.column(options.group);
        Map<String, List<Point>> result = new LinkedHashMap<>();
        for (String[] row : table.rows()) {
            Double xValue = number(row[x], options.x);
            Double yValue = number(row[y], options.y);
            Double sizeValue = withSize ? number(row[size], options.size) : 1.0;
            String label = group < 0 ? "data" : text(row[group]);
            if (xValue == null || yValue == null || sizeValue == null || label == null) {
                missing(options, omitted, "mapped plot columns");
                continue;
            }
            result.computeIfAbsent(label, ignored -> new ArrayList<>())
                .add(new Point(xValue, yValue, sizeValue));
        }
        if (result.isEmpty()) throw new IllegalArgumentException("plot has no complete rows");
        return result;
    }

    private static double[] numbers(DelimitedData table, String column,
            Options options, Counter omitted) {
        int index = table.column(column);
        List<Double> result = new ArrayList<>();
        for (String[] row : table.rows()) {
            Double value = number(row[index], column);
            if (value == null) {
                missing(options, omitted, column);
                continue;
            }
            result.add(value);
        }
        if (result.isEmpty()) throw new IllegalArgumentException("plot has no complete rows");
        return result.stream().mapToDouble(Double::doubleValue).toArray();
    }

    private static Aggregation aggregate(DelimitedData table, String category,
            String value, String group, boolean count, Options options,
            Counter omitted) {
        int categoryIndex = table.column(category);
        int valueIndex = value == null ? -1 : table.column(value);
        int groupIndex = group == null ? -1 : table.column(group);
        Set<String> categories = new LinkedHashSet<>();
        Map<String, Map<String, Summary>> groups = new LinkedHashMap<>();
        for (String[] row : table.rows()) {
            String categoryValue = text(row[categoryIndex]);
            String groupValue = groupIndex < 0 ? "data" : text(row[groupIndex]);
            Double number = count ? 1.0 : number(row[valueIndex], value);
            if (categoryValue == null || groupValue == null || number == null) {
                missing(options, omitted, "mapped plot columns");
                continue;
            }
            categories.add(categoryValue);
            groups.computeIfAbsent(groupValue, ignored -> new LinkedHashMap<>())
                .computeIfAbsent(categoryValue, ignored -> new Summary()).add(number);
        }
        if (categories.isEmpty()) throw new IllegalArgumentException("plot has no complete rows");
        return new Aggregation(List.copyOf(categories), groups);
    }

    private static Double number(String value, String column) {
        String normalized = value.trim();
        if (normalized.isEmpty() || normalized.equalsIgnoreCase("NA")
                || normalized.equalsIgnoreCase("N/A")
                || normalized.equalsIgnoreCase("null")
                || normalized.equals(".")) return null;
        try {
            double result = Double.parseDouble(normalized);
            return Double.isFinite(result) ? result : null;
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException("non-numeric value in " + column
                + ": " + value);
        }
    }

    private static String text(String value) {
        String result = value.trim();
        return result.isEmpty() || result.equalsIgnoreCase("NA") ? null : result;
    }

    private static void missing(Options options, Counter omitted, String columns) {
        if (!options.omitMissing)
            throw new IllegalArgumentException("missing value in " + columns
                + " (use --na omit to drop incomplete rows)");
        omitted.value++;
    }

    private static String require(String value, String option) {
        if (value == null || value.isBlank())
            throw new IllegalArgumentException(option + " is required for this plot type");
        return value;
    }

    private static String help() {
        return """
            Usage: java -jar jlinalg-<version>.jar plot --type TYPE --input TABLE --out PLOT
              Mappings: --x COLUMN [--y COLUMN] [--group COLUMN] [--size COLUMN]
              Types: scatter, line, connected-scatter, area, step, dot, histogram,
                density, box, violin, bar, pie, donut, bubble, heatmap, spider,
                circular-bar, parallel, dendrogram, treemap, network, chord,
                edge-bundling, sankey, venn
              Bar: omit --y for counts; use --group COLUMN and --position grouped|stacked
              Box/violin: --x GROUP --y VALUE (or --group GROUP --y VALUE)
              Bubble: --x X --y Y --size SIZE [--group GROUP]
              Heatmap: --x COLUMN --y ROW --fill VALUE
                [--cluster none|rows|columns|both]
              Spider: --x AXIS --y VALUE [--group SERIES]; values must be in [0,1]
              Parallel/dendrogram: --id LABEL --columns X,Y,Z [--group GROUP]
                [--orientation vertical|horizontal]
              Network/chord/edge-bundling/sankey: --source SOURCE --target TARGET
                [--weight WEIGHT] [--group EDGE_GROUP]
              Venn (two sets): --x ITEM --group SET; repeat an item in both sets
              Export: choose .png, .jpg/.jpeg, .bmp, .gif, .tif/.tiff, .svg, .eps or .pdf
              Style: [--title TEXT] [--x-label TEXT] [--y-label TEXT]
                [--width 800] [--height 600]
                [--theme default|ggplot2|minimal|classic]
                [--legend auto|none|right|bottom|inside-ne|inside-nw|inside-se|inside-sw]
                [--no-grid] [--log-x] [--log-y]
                [--x-range MIN,MAX] [--y-range MIN,MAX]
                [--number-pattern PATTERN] [--point-size PIXELS]
                [--point-shape circle|square|diamond|triangle-up|triangle-down|cross|plus|none]
                [--line-width WIDTH] [--line-pattern solid|dashed|dotted|dash-dot]
                [--colors '#0072B2,#D55E00,...'] [--bins 30] [--na omit|error]
            """;
    }

    private record Point(double x, double y, double size) { }
    private record Cell(String row, String column, double value) { }
    private static final class Counter { int value; }

    private static final class Summary {
        private double sum;
        private int count;
        void add(double value) { sum += value; count++; }
        double sum() { return sum; }
        double mean() { return sum / count; }
    }

    private record Aggregation(List<String> categories,
            Map<String, Map<String, Summary>> groups) {
        double[] singleValues(boolean mean) {
            return values(groups.values().iterator().next(), mean);
        }

        Map<String, double[]> series(boolean mean) {
            Map<String, double[]> result = new LinkedHashMap<>();
            groups.forEach((name, summaries) -> result.put(name, values(summaries, mean)));
            return result;
        }

        private double[] values(Map<String, Summary> summaries, boolean mean) {
            double[] result = new double[categories.size()];
            for (int index = 0; index < result.length; index++) {
                Summary summary = summaries.get(categories.get(index));
                result[index] = summary == null ? 0.0 : mean ? summary.mean() : summary.sum();
            }
            return result;
        }
    }

    private static final class Options {
        private String type;
        private Path input;
        private Path output;
        private String x;
        private String y;
        private String group;
        private String size;
        private String fill;
        private String id;
        private List<String> columns;
        private String source;
        private String target;
        private String weight;
        private String cluster = "none";
        private String orientation = "vertical";
        private String title = "";
        private String xLabel;
        private String yLabel;
        private String position = "grouped";
        private PlotTheme theme = PlotTheme.GGPLOT2;
        private LegendPosition legend = LegendPosition.AUTO;
        private int width = 800;
        private int height = 600;
        private int bins = 30;
        private int markerSize = 8;
        private float lineWidth = 2.0f;
        private LinePattern linePattern = LinePattern.SOLID;
        private PointShape pointShape = PointShape.CIRCLE;
        private boolean grid = true;
        private boolean logX;
        private boolean logY;
        private boolean omitMissing = true;
        private double[] xRange;
        private double[] yRange;
        private String numberPattern = "0.####";
        private Color[] colors;
        private boolean help;

        static Options parse(String[] arguments) {
            Options result = new Options();
            for (int index = 0; index < arguments.length; index++) {
                String option = arguments[index];
                switch (option) {
                    case "--help", "-h" -> result.help = true;
                    case "--type" -> result.type = value(arguments, ++index, option)
                        .toLowerCase(Locale.ROOT);
                    case "--input", "--data" -> result.input = Path.of(
                        value(arguments, ++index, option));
                    case "--out", "--output" -> result.output = Path.of(
                        value(arguments, ++index, option));
                    case "--x" -> result.x = value(arguments, ++index, option);
                    case "--y" -> result.y = value(arguments, ++index, option);
                    case "--group" -> result.group = value(arguments, ++index, option);
                    case "--size" -> result.size = value(arguments, ++index, option);
                    case "--fill" -> result.fill = value(arguments, ++index, option);
                    case "--id" -> result.id = value(arguments, ++index, option);
                    case "--columns" -> result.columns = names(
                        value(arguments, ++index, option), option);
                    case "--source" -> result.source = value(arguments, ++index, option);
                    case "--target" -> result.target = value(arguments, ++index, option);
                    case "--weight" -> result.weight = value(arguments, ++index, option);
                    case "--cluster" -> result.cluster = value(arguments, ++index, option)
                        .toLowerCase(Locale.ROOT);
                    case "--orientation" -> result.orientation = value(
                        arguments, ++index, option).toLowerCase(Locale.ROOT);
                    case "--title" -> result.title = value(arguments, ++index, option);
                    case "--x-label" -> result.xLabel = value(arguments, ++index, option);
                    case "--y-label" -> result.yLabel = value(arguments, ++index, option);
                    case "--position" -> result.position = value(arguments, ++index, option)
                        .toLowerCase(Locale.ROOT);
                    case "--theme" -> result.theme = PlotTheme.valueOf(
                        value(arguments, ++index, option).toUpperCase(Locale.ROOT));
                    case "--legend" -> result.legend = legend(
                        value(arguments, ++index, option));
                    case "--width" -> result.width = integer(arguments, ++index, option);
                    case "--height" -> result.height = integer(arguments, ++index, option);
                    case "--bins" -> result.bins = integer(arguments, ++index, option);
                    case "--marker-size", "--point-size", "--dot-size" ->
                        result.markerSize = integer(
                        arguments, ++index, option);
                    case "--line-width" -> result.lineWidth = floating(
                        arguments, ++index, option);
                    case "--line-pattern" -> result.linePattern = linePattern(
                        value(arguments, ++index, option));
                    case "--point-shape", "--dot-pattern" -> result.pointShape = pointShape(
                        value(arguments, ++index, option));
                    case "--no-grid" -> result.grid = false;
                    case "--log-x" -> result.logX = true;
                    case "--log-y" -> result.logY = true;
                    case "--x-range" -> result.xRange = range(
                        value(arguments, ++index, option), option);
                    case "--y-range" -> result.yRange = range(
                        value(arguments, ++index, option), option);
                    case "--number-pattern" -> result.numberPattern = value(
                        arguments, ++index, option);
                    case "--colors" -> result.colors = colors(
                        value(arguments, ++index, option));
                    case "--na" -> result.omitMissing = missingMode(
                        value(arguments, ++index, option));
                    default -> throw new IllegalArgumentException(
                        "unknown plot option: " + option);
                }
            }
            if (result.help) return result;
            require(result.type, "--type");
            if (result.input == null) throw new IllegalArgumentException("--input is required");
            if (result.output == null) throw new IllegalArgumentException("--out is required");
            if (!Set.of("scatter", "line", "connected-scatter", "area", "step",
                    "dot", "histogram", "density", "box", "violin", "bar", "pie",
                    "donut", "bubble", "heatmap", "spider", "radar", "circular-bar",
                    "parallel", "dendrogram", "treemap", "network", "chord",
                    "edge-bundling", "sankey", "venn").contains(result.type))
                throw new IllegalArgumentException("unsupported plot type: " + result.type);
            if (!Set.of("grouped", "stacked").contains(result.position))
                throw new IllegalArgumentException("--position must be grouped or stacked");
            if (result.bins < 1) throw new IllegalArgumentException("--bins must be positive");
            result.clusterMode();
            result.orientation();
            return result;
        }

        PlotOptions plotOptions(String defaultX, String defaultY) {
            return plotOptions(defaultX, defaultY, null);
        }

        PlotOptions plotOptions(String defaultX, String defaultY,
                Color[] fallbackColors) {
            PlotOptions.Builder builder = PlotOptions.builder().size(width, height)
                .title(title).xLabel(xLabel == null ? defaultX : xLabel)
                .yLabel(yLabel == null ? defaultY : yLabel).theme(theme)
                .legend(legend).grid(grid).logarithmicX(logX).logarithmicY(logY)
                .numberPattern(numberPattern).markerSize(markerSize)
                .lineWidth(lineWidth).linePattern(linePattern).pointShape(pointShape);
            if (xRange != null) builder.xRange(xRange[0], xRange[1]);
            if (yRange != null) builder.yRange(yRange[0], yRange[1]);
            if (colors != null) builder.colors(colors);
            else if (fallbackColors != null) builder.colors(fallbackColors);
            return builder.build();
        }

        DiagramPlots.ClusterMode clusterMode() {
            return switch (cluster) {
                case "none" -> DiagramPlots.ClusterMode.NONE;
                case "rows", "vertical" -> DiagramPlots.ClusterMode.ROWS;
                case "columns", "horizontal" -> DiagramPlots.ClusterMode.COLUMNS;
                case "both" -> DiagramPlots.ClusterMode.BOTH;
                default -> throw new IllegalArgumentException(
                    "--cluster must be none, rows, columns, or both");
            };
        }

        DiagramPlots.Orientation orientation() {
            return switch (orientation) {
                case "vertical" -> DiagramPlots.Orientation.VERTICAL;
                case "horizontal" -> DiagramPlots.Orientation.HORIZONTAL;
                default -> throw new IllegalArgumentException(
                    "--orientation must be vertical or horizontal");
            };
        }

        private static String value(String[] arguments, int index, String option) {
            if (index >= arguments.length) throw new IllegalArgumentException(
                option + " needs a value");
            return arguments[index];
        }

        private static int integer(String[] arguments, int index, String option) {
            try { return Integer.parseInt(value(arguments, index, option)); }
            catch (NumberFormatException exception) {
                throw new IllegalArgumentException(option + " needs an integer");
            }
        }

        private static float floating(String[] arguments, int index, String option) {
            try { return Float.parseFloat(value(arguments, index, option)); }
            catch (NumberFormatException exception) {
                throw new IllegalArgumentException(option + " needs a number");
            }
        }

        private static List<String> names(String value, String option) {
            List<String> result = Arrays.stream(value.split(",", -1))
                .map(String::trim).toList();
            if (result.isEmpty() || result.stream().anyMatch(String::isEmpty)
                    || new LinkedHashSet<>(result).size() != result.size())
                throw new IllegalArgumentException(
                    option + " needs unique comma-separated column names");
            return result;
        }

        private static double[] range(String value, String option) {
            String[] fields = value.split(",", -1);
            if (fields.length != 2) throw new IllegalArgumentException(
                option + " needs MIN,MAX");
            try { return new double[] {Double.parseDouble(fields[0]),
                Double.parseDouble(fields[1])}; }
            catch (NumberFormatException exception) {
                throw new IllegalArgumentException(option + " needs numeric bounds");
            }
        }

        private static Color[] colors(String value) {
            try {
                return Arrays.stream(value.split(",", -1)).map(String::trim)
                    .map(Color::decode).toArray(Color[]::new);
            } catch (NumberFormatException exception) {
                throw new IllegalArgumentException("--colors needs comma-separated hex colors");
            }
        }

        private static boolean missingMode(String value) {
            if (value.equalsIgnoreCase("omit")) return true;
            if (value.equalsIgnoreCase("error")) return false;
            throw new IllegalArgumentException("--na must be omit or error");
        }

        private static LegendPosition legend(String value) {
            return switch (value.toLowerCase(Locale.ROOT)) {
                case "auto" -> LegendPosition.AUTO;
                case "none" -> LegendPosition.NONE;
                case "right" -> LegendPosition.RIGHT;
                case "bottom" -> LegendPosition.BOTTOM;
                case "inside-ne" -> LegendPosition.INSIDE_NORTH_EAST;
                case "inside-nw" -> LegendPosition.INSIDE_NORTH_WEST;
                case "inside-se" -> LegendPosition.INSIDE_SOUTH_EAST;
                case "inside-sw" -> LegendPosition.INSIDE_SOUTH_WEST;
                default -> throw new IllegalArgumentException("unknown legend position: " + value);
            };
        }

        private static LinePattern linePattern(String value) {
            return switch (value.toLowerCase(Locale.ROOT)) {
                case "solid" -> LinePattern.SOLID;
                case "dash", "dashed" -> LinePattern.DASHED;
                case "dot", "dotted" -> LinePattern.DOTTED;
                case "dash-dot", "dot-dash" -> LinePattern.DASH_DOT;
                default -> throw new IllegalArgumentException("unknown line pattern: " + value);
            };
        }

        private static PointShape pointShape(String value) {
            return switch (value.toLowerCase(Locale.ROOT)) {
                case "circle" -> PointShape.CIRCLE;
                case "square" -> PointShape.SQUARE;
                case "diamond" -> PointShape.DIAMOND;
                case "triangle", "triangle-up" -> PointShape.TRIANGLE_UP;
                case "triangle-down" -> PointShape.TRIANGLE_DOWN;
                case "cross" -> PointShape.CROSS;
                case "plus" -> PointShape.PLUS;
                case "none" -> PointShape.NONE;
                default -> throw new IllegalArgumentException("unknown point shape: " + value);
            };
        }
    }
}
