/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.geom.Path2D;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.knowm.xchart.internal.chartpart.IChart;

final class SpiderChart implements IChart {
    private final List<String> axes;
    private final Map<String, double[]> series;
    private final PlotOptions options;

    SpiderChart(List<String> axisLabels, Map<String, double[]> source,
            PlotOptions options) {
        if (axisLabels == null || axisLabels.size() < 3)
            throw new IllegalArgumentException("spider plot needs at least three axes");
        if (source == null || source.isEmpty())
            throw new IllegalArgumentException("spider plot series are required");
        axes = List.copyOf(axisLabels);
        series = new LinkedHashMap<>();
        source.forEach((name, values) -> {
            if (name == null || name.isBlank() || values == null
                    || values.length != axes.size())
                throw new IllegalArgumentException("spider series must match every named axis");
            double[] copy = values.clone();
            for (double value : copy)
                if (!Double.isFinite(value) || value < 0.0 || value > 1.0)
                    throw new IllegalArgumentException(
                        "spider values must be finite and in [0,1]");
            series.put(name, copy);
        });
        this.options = options;
    }

    @Override public int getWidth() { return options.width(); }
    @Override public int getHeight() { return options.height(); }
    @Override public String getTitle() { return options.title(); }

    @Override
    public void paint(Graphics2D source, int width, int height) {
        Graphics2D graphics = DiagramSupport.begin(source, options, width, height);
        try {
            int top = DiagramSupport.top(options);
            int legendWidth = series.size() > 1 ? 125 : 25;
            double centerX = (width - legendWidth) / 2.0;
            double centerY = top + (height - top - 30) / 2.0;
            double radius = Math.max(20.0, Math.min(width - legendWidth,
                height - top - 30) * 0.36);
            graphics.setStroke(PlotSupport.strokeWidth(1f));
            for (int ring = 1; ring <= 5; ring++) {
                Path2D path = polygon(centerX, centerY, radius * ring / 5.0,
                    uniform(axes.size(), 1.0));
                graphics.setColor(DiagramSupport.GRID);
                graphics.draw(path);
            }
            graphics.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
            FontMetrics metrics = graphics.getFontMetrics();
            for (int axis = 0; axis < axes.size(); axis++) {
                double angle = angle(axis);
                double x = centerX + radius * Math.cos(angle);
                double y = centerY + radius * Math.sin(angle);
                graphics.setColor(DiagramSupport.GRID);
                graphics.drawLine((int) centerX, (int) centerY, (int) x, (int) y);
                String label = axes.get(axis);
                double labelX = centerX + (radius + 15) * Math.cos(angle);
                double labelY = centerY + (radius + 15) * Math.sin(angle);
                graphics.setColor(DiagramSupport.TEXT);
                graphics.drawString(label, (int) labelX
                    - (Math.cos(angle) < -0.2 ? metrics.stringWidth(label)
                        : Math.abs(Math.cos(angle)) < 0.2 ? metrics.stringWidth(label) / 2 : 0),
                    (int) labelY + metrics.getAscent() / 2);
            }
            int index = 0;
            for (Map.Entry<String, double[]> entry : series.entrySet()) {
                Color color = DiagramSupport.color(options, index++);
                Path2D path = polygon(centerX, centerY, radius, entry.getValue());
                graphics.setColor(DiagramSupport.translucent(color, 55));
                graphics.fill(path);
                graphics.setColor(color);
                graphics.setStroke(PlotSupport.stroke(options));
                graphics.draw(path);
                for (int axis = 0; axis < axes.size(); axis++) {
                    double angle = angle(axis);
                    DiagramSupport.marker(graphics, options.pointShape(),
                        centerX + radius * entry.getValue()[axis] * Math.cos(angle),
                        centerY + radius * entry.getValue()[axis] * Math.sin(angle),
                        options.markerSize(), color);
                }
            }
            DiagramSupport.legend(graphics, options,
                new ArrayList<>(series.keySet()), width - legendWidth + 8, top + 15);
        } finally {
            graphics.dispose();
        }
    }

    private Path2D polygon(double centerX, double centerY, double radius,
            double[] values) {
        Path2D path = new Path2D.Double();
        for (int index = 0; index < values.length; index++) {
            double angle = angle(index);
            double x = centerX + radius * values[index] * Math.cos(angle);
            double y = centerY + radius * values[index] * Math.sin(angle);
            if (index == 0) path.moveTo(x, y); else path.lineTo(x, y);
        }
        path.closePath();
        return path;
    }

    private double angle(int index) {
        return -Math.PI / 2.0 + 2.0 * Math.PI * index / axes.size();
    }

    private static double[] uniform(int length, double value) {
        double[] result = new double[length];
        java.util.Arrays.fill(result, value);
        return result;
    }
}
