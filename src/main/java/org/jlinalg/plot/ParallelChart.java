/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.geom.Path2D;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.knowm.xchart.internal.chartpart.IChart;

final class ParallelChart implements IChart {
    private final List<String> labels;
    private final List<String> dimensions;
    private final double[][] values;
    private final List<String> observationGroups;
    private final List<String> groups;
    private final double[] minima;
    private final double[] maxima;
    private final PlotOptions options;

    ParallelChart(DiagramPlots.ObservationData data, PlotOptions options) {
        validate(data, 2);
        labels = List.copyOf(data.labels());
        dimensions = List.copyOf(data.dimensions());
        values = copy(data.values());
        observationGroups = data.groups() == null || data.groups().isEmpty()
            ? java.util.Collections.nCopies(labels.size(), "data")
            : List.copyOf(data.groups());
        if (observationGroups.size() != labels.size())
            throw new IllegalArgumentException("parallel group count must match observations");
        Map<String, Boolean> ordered = new LinkedHashMap<>();
        observationGroups.forEach(group -> ordered.put(group, Boolean.TRUE));
        groups = List.copyOf(ordered.keySet());
        minima = new double[dimensions.size()];
        maxima = new double[dimensions.size()];
        java.util.Arrays.fill(minima, Double.POSITIVE_INFINITY);
        java.util.Arrays.fill(maxima, Double.NEGATIVE_INFINITY);
        for (double[] row : values)
            for (int column = 0; column < row.length; column++) {
                minima[column] = Math.min(minima[column], row[column]);
                maxima[column] = Math.max(maxima[column], row[column]);
            }
        this.options = options;
    }

    @Override public int getWidth() { return options.width(); }
    @Override public int getHeight() { return options.height(); }
    @Override public String getTitle() { return options.title(); }

    @Override
    public void paint(Graphics2D source, int width, int height) {
        Graphics2D graphics = DiagramSupport.begin(source, options, width, height);
        try {
            int top = DiagramSupport.top(options) + 22;
            int left = 65, right = groups.size() > 1 ? 145 : 50, bottom = 58;
            double plotWidth = width - left - right;
            double plotHeight = height - top - bottom;
            double step = plotWidth / (dimensions.size() - 1.0);
            graphics.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
            FontMetrics metrics = graphics.getFontMetrics();
            DecimalFormat format = new DecimalFormat(options.numberPattern());
            for (int dimension = 0; dimension < dimensions.size(); dimension++) {
                int x = (int) Math.round(left + dimension * step);
                graphics.setColor(DiagramSupport.GRID);
                graphics.setStroke(PlotSupport.strokeWidth(1f));
                graphics.drawLine(x, top, x, (int) (top + plotHeight));
                graphics.setColor(DiagramSupport.TEXT);
                String label = DiagramSupport.shorten(dimensions.get(dimension), metrics, 90);
                graphics.drawString(label, x - metrics.stringWidth(label) / 2, top - 8);
                String high = format.format(maxima[dimension]);
                String low = format.format(minima[dimension]);
                graphics.drawString(high, x - metrics.stringWidth(high) / 2, top + 12);
                graphics.drawString(low, x - metrics.stringWidth(low) / 2,
                    (int) (top + plotHeight) + 18);
            }
            for (int observation = 0; observation < values.length; observation++) {
                Color color = DiagramSupport.color(options,
                    groups.indexOf(observationGroups.get(observation)));
                Path2D path = new Path2D.Double();
                for (int dimension = 0; dimension < dimensions.size(); dimension++) {
                    double fraction = maxima[dimension] == minima[dimension] ? 0.5
                        : (values[observation][dimension] - minima[dimension])
                            / (maxima[dimension] - minima[dimension]);
                    double x = left + dimension * step;
                    double y = top + plotHeight * (1.0 - fraction);
                    if (dimension == 0) path.moveTo(x, y); else path.lineTo(x, y);
                }
                graphics.setColor(DiagramSupport.translucent(color, 135));
                graphics.setStroke(PlotSupport.stroke(options));
                graphics.draw(path);
                for (int dimension = 0; dimension < dimensions.size(); dimension++) {
                    double fraction = maxima[dimension] == minima[dimension] ? 0.5
                        : (values[observation][dimension] - minima[dimension])
                            / (maxima[dimension] - minima[dimension]);
                    DiagramSupport.marker(graphics, options.pointShape(),
                        left + dimension * step, top + plotHeight * (1.0 - fraction),
                        options.markerSize(), color);
                }
            }
            DiagramSupport.legend(graphics, options, groups,
                width - right + 24, top + 10);
        } finally {
            graphics.dispose();
        }
    }

    static void validate(DiagramPlots.ObservationData data, int minimumDimensions) {
        if (data == null || data.labels() == null || data.dimensions() == null
                || data.values() == null || data.labels().isEmpty()
                || data.dimensions().size() < minimumDimensions
                || data.values().length != data.labels().size())
            throw new IllegalArgumentException("observation matrix and labels are required");
        int dimensions = data.dimensions().size();
        for (double[] row : data.values()) {
            if (row == null || row.length != dimensions)
                throw new IllegalArgumentException("observation matrix must be rectangular");
            for (double value : row)
                if (!Double.isFinite(value))
                    throw new IllegalArgumentException("observation values must be finite");
        }
    }

    static double[][] copy(double[][] source) {
        double[][] result = new double[source.length][];
        for (int index = 0; index < source.length; index++)
            result[index] = source[index].clone();
        return result;
    }
}
