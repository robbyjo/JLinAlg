/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import org.knowm.xchart.BoxChart;
import org.knowm.xchart.BoxChartBuilder;
import org.knowm.xchart.BubbleChart;
import org.knowm.xchart.BubbleChartBuilder;
import org.knowm.xchart.CategoryChart;
import org.knowm.xchart.CategoryChartBuilder;
import org.knowm.xchart.CategorySeries;
import org.knowm.xchart.Histogram;
import org.knowm.xchart.PieChart;
import org.knowm.xchart.PieChartBuilder;
import org.knowm.xchart.PieSeries;
import org.knowm.xchart.style.markers.SeriesMarkers;

/** Convenience factories for common statistical charts. */
public final class Plots {
    private Plots() { }

    public static XYPlotBuilder xy() { return new XYPlotBuilder(); }

    public static Plot scatter(double[] x, double[] y, PlotOptions options) {
        return xy().options(options).points("data", x, y).build();
    }

    public static Plot line(double[] x, double[] y, PlotOptions options) {
        return xy().options(options).line("data", x, y).build();
    }

    public static Plot area(double[] x, double[] y, PlotOptions options) {
        return xy().options(options).area("data", x, y).build();
    }

    public static Plot bubble(double[] xValues, double[] yValues, double[] sizes,
            PlotOptions options) {
        return bubble(Map.of("data", new BubbleData(xValues, yValues, sizes)), options);
    }

    /** Coordinates and positive sizes for one bubble-plot series. */
    public record BubbleData(double[] x, double[] y, double[] sizes) { }

    public static Plot bubble(Map<String, BubbleData> groups, PlotOptions options) {
        if (groups == null || groups.isEmpty())
            throw new IllegalArgumentException("bubble groups are required");
        BubbleChart chart = new BubbleChartBuilder().width(options.width())
            .height(options.height()).title(options.title())
            .xAxisTitle(options.xLabel()).yAxisTitle(options.yLabel())
            .theme(PlotSupport.theme(options)).build();
        for (Map.Entry<String, BubbleData> entry : groups.entrySet()) {
            BubbleData data = Objects.requireNonNull(entry.getValue());
            double[] x = PlotSupport.finite(data.x(), "x");
            double[] y = PlotSupport.finite(data.y(), "y");
            double[] bubbles = PlotSupport.finite(data.sizes(), "bubble size");
            PlotSupport.sameLength(x, y, "x and y");
            PlotSupport.sameLength(x, bubbles, "x and bubble size");
            for (double size : bubbles)
                if (size <= 0.0)
                    throw new IllegalArgumentException("bubble sizes must be positive");
            chart.addSeries(entry.getKey(), x, y, bubbles)
                .setShowInLegend(groups.size() > 1);
        }
        PlotSupport.style(chart, options, groups.size() > 1);
        return new Plot(chart);
    }

    public static Plot histogram(double[] values, int bins, PlotOptions options) {
        double[] data = PlotSupport.finite(values, "histogram");
        if (bins < 1) throw new IllegalArgumentException("histogram bins must be positive");
        Histogram histogram = new Histogram(Arrays.stream(data).boxed().toList(), bins);
        CategoryChart chart = categoryChart(options);
        chart.addSeries("count", histogram.getxAxisData(), histogram.getyAxisData())
            .setShowInLegend(false);
        chart.getStyler().setAvailableSpaceFill(1.0);
        PlotSupport.style(chart, options, false);
        return new Plot(chart);
    }

    public static Plot density(double[] values, PlotOptions options) {
        DensityEstimate estimate = DensityEstimate.of(PlotSupport.finite(values, "density"), 256);
        return xy().options(options).area("density", estimate.coordinates(),
            estimate.density()).build();
    }

    public static Plot dotPlot(List<String> categories, double[] values,
            PlotOptions options) {
        List<String> labels = labels(categories, values.length);
        double[] data = PlotSupport.finite(values, "dot plot");
        CategoryChart chart = categoryChart(options);
        CategorySeries series = chart.addSeries("data", labels,
            Arrays.stream(data).boxed().toList());
        series.setChartCategorySeriesRenderStyle(
            CategorySeries.CategorySeriesRenderStyle.Scatter);
        series.setMarker(SeriesMarkers.CIRCLE).setShowInLegend(false);
        PlotSupport.style(chart, options, false);
        return new Plot(chart);
    }

    public static Plot barChart(List<String> categories, double[] values,
            PlotOptions options) {
        List<String> labels = labels(categories, values.length);
        double[] data = PlotSupport.finite(values, "bar");
        CategoryChart chart = categoryChart(options);
        chart.addSeries("data", labels, Arrays.stream(data).boxed().toList())
            .setShowInLegend(false);
        PlotSupport.style(chart, options, false);
        return new Plot(chart);
    }

    public static Plot groupedBarChart(List<String> categories,
            Map<String, double[]> series, PlotOptions options) {
        return multipleBarChart(categories, series, options, false);
    }

    public static Plot stackedBarChart(List<String> categories,
            Map<String, double[]> series, PlotOptions options) {
        return multipleBarChart(categories, series, options, true);
    }

    private static Plot multipleBarChart(List<String> categories,
            Map<String, double[]> values, PlotOptions options, boolean stacked) {
        if (values == null || values.isEmpty())
            throw new IllegalArgumentException("bar series are required");
        List<String> labels = labels(categories, categories.size());
        CategoryChart chart = categoryChart(options);
        for (Map.Entry<String, double[]> entry : values.entrySet()) {
            if (entry.getKey() == null || entry.getKey().isBlank())
                throw new IllegalArgumentException("bar series names are required");
            double[] data = PlotSupport.finite(entry.getValue(), "bar");
            if (data.length != labels.size())
                throw new IllegalArgumentException("all bar series must match category count");
            chart.addSeries(entry.getKey(), labels, Arrays.stream(data).boxed().toList());
        }
        chart.getStyler().setStacked(stacked);
        PlotSupport.style(chart, options, values.size() > 1);
        return new Plot(chart);
    }

    public static Plot boxPlot(Map<String, double[]> groups, PlotOptions options) {
        Map<String, double[]> data = groups(groups, "box plot");
        BoxChart chart = new BoxChartBuilder().width(options.width())
            .height(options.height()).title(options.title())
            .xAxisTitle(options.xLabel()).yAxisTitle(options.yLabel())
            .theme(PlotSupport.theme(options)).build();
        data.forEach(chart::addSeries);
        PlotSupport.style(chart, options, false);
        chart.getStyler().setLegendVisible(false);
        return new Plot(chart);
    }

    public static Plot violinPlot(Map<String, double[]> groups, PlotOptions options) {
        return new Plot(new ViolinChart(groups(groups, "violin plot"), options));
    }

    public static Plot pieChart(Map<String, ? extends Number> slices,
            PlotOptions options) {
        return pie(slices, options, false);
    }

    public static Plot donutChart(Map<String, ? extends Number> slices,
            PlotOptions options) {
        return pie(slices, options, true);
    }

    private static Plot pie(Map<String, ? extends Number> slices,
            PlotOptions options, boolean donut) {
        if (slices == null || slices.isEmpty())
            throw new IllegalArgumentException("pie slices are required");
        PieChart chart = new PieChartBuilder().width(options.width())
            .height(options.height()).title(options.title())
            .theme(PlotSupport.theme(options)).build();
        double total = 0.0;
        for (Map.Entry<String, ? extends Number> entry : slices.entrySet()) {
            String name = entry.getKey();
            double value = entry.getValue() == null ? Double.NaN
                : entry.getValue().doubleValue();
            if (name == null || name.isBlank() || !Double.isFinite(value) || value < 0.0)
                throw new IllegalArgumentException("pie slices need names and finite nonnegative values");
            PieSeries series = chart.addSeries(name, value);
            if (donut) series.setChartPieSeriesRenderStyle(
                PieSeries.PieSeriesRenderStyle.Donut);
            total += value;
        }
        if (!(total > 0.0)) throw new IllegalArgumentException("pie slice total must be positive");
        PlotSupport.style(chart, options, true);
        chart.getStyler().setLabelsVisible(true);
        return new Plot(chart);
    }

    private static CategoryChart categoryChart(PlotOptions options) {
        return new CategoryChartBuilder().width(options.width()).height(options.height())
            .title(options.title()).xAxisTitle(options.xLabel())
            .yAxisTitle(options.yLabel()).theme(PlotSupport.theme(options)).build();
    }

    private static List<String> labels(List<String> values, int expected) {
        Objects.requireNonNull(values, "categories");
        if (values.size() != expected)
            throw new IllegalArgumentException("category and value counts must match");
        List<String> result = new ArrayList<>(values.size());
        for (String value : values) {
            if (value == null || value.isBlank())
                throw new IllegalArgumentException("category labels must be nonblank");
            result.add(value);
        }
        return List.copyOf(result);
    }

    private static Map<String, double[]> groups(Map<String, double[]> source,
            String name) {
        if (source == null || source.isEmpty())
            throw new IllegalArgumentException(name + " groups are required");
        Map<String, double[]> result = new LinkedHashMap<>();
        source.forEach((label, values) -> {
            if (label == null || label.isBlank())
                throw new IllegalArgumentException(name + " group names must be nonblank");
            result.put(label, PlotSupport.finite(values, name));
        });
        return result;
    }

    record DensityEstimate(double[] coordinates, double[] density) {
        static DensityEstimate of(double[] values, int points) {
            double[] sorted = values.clone();
            Arrays.sort(sorted);
            double minimum = sorted[0], maximum = sorted[sorted.length - 1];
            double spread = maximum - minimum;
            if (!(spread > 0.0)) spread = Math.max(1.0, Math.abs(minimum) * 0.1);
            double standardDeviation = standardDeviation(sorted);
            double bandwidth = 1.06 * standardDeviation * Math.pow(sorted.length, -0.2);
            if (!(bandwidth > 0.0)) bandwidth = spread / 20.0;
            minimum -= 3.0 * bandwidth;
            maximum += 3.0 * bandwidth;
            double[] x = new double[points], density = new double[points];
            double normalizer = sorted.length * bandwidth * Math.sqrt(2.0 * Math.PI);
            for (int index = 0; index < points; index++) {
                x[index] = minimum + (maximum - minimum) * index / (points - 1.0);
                double sum = 0.0;
                for (double value : sorted) {
                    double z = (x[index] - value) / bandwidth;
                    sum += Math.exp(-0.5 * z * z);
                }
                density[index] = sum / normalizer;
            }
            return new DensityEstimate(x, density);
        }

        private static double standardDeviation(double[] values) {
            double mean = Arrays.stream(values).average().orElseThrow();
            double sum = 0.0;
            for (double value : values) sum += (value - mean) * (value - mean);
            return values.length > 1 ? Math.sqrt(sum / (values.length - 1.0)) : 0.0;
        }
    }
}
