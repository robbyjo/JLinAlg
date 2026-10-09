/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.AffineTransform;
import java.awt.geom.Path2D;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import org.knowm.xchart.internal.chartpart.IChart;

final class ViolinChart implements IChart {
    private final PlotOptions options;
    private final List<Violin> violins;
    private final double minimum;
    private final double maximum;

    ViolinChart(Map<String, double[]> groups, PlotOptions options) {
        this.options = options;
        if (options.yLogarithmic())
            throw new IllegalArgumentException("violin plots do not support a logarithmic value axis");
        violins = new ArrayList<>();
        double low = Double.POSITIVE_INFINITY, high = Double.NEGATIVE_INFINITY;
        for (Map.Entry<String, double[]> entry : groups.entrySet()) {
            double[] values = entry.getValue().clone();
            Arrays.sort(values);
            Plots.DensityEstimate estimate = Plots.DensityEstimate.of(values, 160);
            violins.add(new Violin(entry.getKey(), values, estimate));
            low = Math.min(low, estimate.coordinates()[0]);
            high = Math.max(high, estimate.coordinates()[estimate.coordinates().length - 1]);
        }
        double span = high - low;
        minimum = options.yMinimum() == null ? low - 0.02 * span : options.yMinimum();
        maximum = options.yMaximum() == null ? high + 0.02 * span : options.yMaximum();
    }

    @Override public int getWidth() { return options.width(); }
    @Override public int getHeight() { return options.height(); }
    @Override public String getTitle() { return options.title(); }

    @Override
    public void paint(Graphics2D graphics, int width, int height) {
        Graphics2D g = (Graphics2D) graphics.create();
        try {
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                RenderingHints.VALUE_ANTIALIAS_ON);
            g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,
                RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
            g.setColor(options.background());
            g.fillRect(0, 0, width, height);
            int left = 82, right = 28, top = options.title().isBlank() ? 28 : 62;
            int bottom = 76;
            double plotWidth = width - left - right;
            double plotHeight = height - top - bottom;
            if (plotWidth < 20 || plotHeight < 20) return;
            drawPanel(g, left, top, plotWidth, plotHeight);
            drawViolins(g, left, top, plotWidth, plotHeight);
            drawLabels(g, width, height, left, top, plotWidth, plotHeight);
        } finally {
            g.dispose();
        }
    }

    private void drawPanel(Graphics2D g, int left, int top, double width,
            double height) {
        Color panel = options.theme() == PlotTheme.GGPLOT2
            ? new Color(0xEB, 0xEB, 0xEB) : options.background();
        g.setColor(panel);
        g.fillRect(left, top, (int) width, (int) height);
        DecimalFormat format = new DecimalFormat(options.numberPattern());
        g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 12));
        FontMetrics metrics = g.getFontMetrics();
        for (int tick = 0; tick <= 5; tick++) {
            double value = minimum + (maximum - minimum) * tick / 5.0;
            int y = coordinate(value, top, height);
            if (options.gridVisible() && options.theme() != PlotTheme.CLASSIC) {
                g.setColor(options.theme() == PlotTheme.GGPLOT2 ? Color.WHITE
                    : new Color(0xE5, 0xE7, 0xEB));
                g.setStroke(new BasicStroke(1f));
                g.drawLine(left, y, (int) (left + width), y);
            }
            String label = format.format(value);
            g.setColor(new Color(0x33, 0x33, 0x33));
            g.drawString(label, left - metrics.stringWidth(label) - 9,
                y + metrics.getAscent() / 2 - 2);
        }
        g.setColor(new Color(0x44, 0x44, 0x44));
        g.drawLine(left, top, left, (int) (top + height));
        g.drawLine(left, (int) (top + height), (int) (left + width),
            (int) (top + height));
    }

    private void drawViolins(Graphics2D g, int left, int top, double width,
            double height) {
        Color[] colors = options.colors();
        double spacing = width / violins.size();
        for (int group = 0; group < violins.size(); group++) {
            Violin violin = violins.get(group);
            double center = left + spacing * (group + 0.5);
            double maximumDensity = Arrays.stream(violin.estimate().density())
                .max().orElse(1.0);
            double scale = spacing * 0.38 / maximumDensity;
            Path2D path = new Path2D.Double();
            double[] values = violin.estimate().coordinates();
            double[] density = violin.estimate().density();
            path.moveTo(center - density[0] * scale,
                coordinate(values[0], top, height));
            for (int index = 1; index < values.length; index++)
                path.lineTo(center - density[index] * scale,
                    coordinate(values[index], top, height));
            for (int index = values.length - 1; index >= 0; index--)
                path.lineTo(center + density[index] * scale,
                    coordinate(values[index], top, height));
            path.closePath();
            Color base = colors[group % colors.length];
            g.setColor(new Color(base.getRed(), base.getGreen(), base.getBlue(), 150));
            g.fill(path);
            g.setColor(base.darker());
            g.setStroke(new BasicStroke(1.4f));
            g.draw(path);
            double q1 = quantile(violin.values(), 0.25);
            double median = quantile(violin.values(), 0.5);
            double q3 = quantile(violin.values(), 0.75);
            g.setColor(new Color(0x22, 0x22, 0x22));
            g.setStroke(new BasicStroke(4f));
            g.drawLine((int) center, coordinate(q1, top, height),
                (int) center, coordinate(q3, top, height));
            g.fillOval((int) center - 4, coordinate(median, top, height) - 4, 8, 8);
        }
    }

    private void drawLabels(Graphics2D g, int canvasWidth, int canvasHeight,
            int left, int top, double width, double height) {
        g.setColor(new Color(0x22, 0x22, 0x22));
        g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 12));
        FontMetrics metrics = g.getFontMetrics();
        double spacing = width / violins.size();
        for (int group = 0; group < violins.size(); group++) {
            String label = violins.get(group).name();
            int x = (int) (left + spacing * (group + 0.5));
            g.drawString(label, x - metrics.stringWidth(label) / 2,
                (int) (top + height) + 23);
        }
        g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 14));
        metrics = g.getFontMetrics();
        g.drawString(options.xLabel(),
            left + (int) (width - metrics.stringWidth(options.xLabel())) / 2,
            canvasHeight - 17);
        AffineTransform transform = g.getTransform();
        g.rotate(-Math.PI / 2.0);
        g.drawString(options.yLabel(),
            -(top + (int) (height + metrics.stringWidth(options.yLabel())) / 2), 20);
        g.setTransform(transform);
        if (!options.title().isBlank()) {
            g.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
            metrics = g.getFontMetrics();
            g.drawString(options.title(),
                (canvasWidth - metrics.stringWidth(options.title())) / 2, 30);
        }
    }

    private int coordinate(double value, int top, double height) {
        return (int) Math.round(top + height * (maximum - value) / (maximum - minimum));
    }

    private static double quantile(double[] sorted, double probability) {
        double position = probability * (sorted.length - 1);
        int low = (int) Math.floor(position), high = (int) Math.ceil(position);
        return sorted[low] + (position - low) * (sorted[high] - sorted[low]);
    }

    private record Violin(String name, double[] values,
            Plots.DensityEstimate estimate) { }
}
