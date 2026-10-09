/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.geom.Arc2D;
import java.awt.geom.Path2D;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.knowm.xchart.internal.chartpart.IChart;

final class CircularBarChart implements IChart {
    private final List<DiagramPlots.TreeItem> bars;
    private final PlotOptions options;
    private final List<String> groups;
    private final double maximum;

    CircularBarChart(List<DiagramPlots.TreeItem> source, PlotOptions options) {
        if (source == null || source.isEmpty())
            throw new IllegalArgumentException("circular bars are required");
        bars = new ArrayList<>();
        Map<String, Boolean> groupOrder = new LinkedHashMap<>();
        double high = 0.0;
        for (DiagramPlots.TreeItem item : source) {
            if (item == null || item.label() == null || item.label().isBlank()
                    || !Double.isFinite(item.value()) || item.value() < 0.0)
                throw new IllegalArgumentException(
                    "circular bars need labels and finite nonnegative values");
            String group = item.group() == null || item.group().isBlank()
                ? "data" : item.group();
            bars.add(new DiagramPlots.TreeItem(item.label(), group, item.value()));
            groupOrder.put(group, Boolean.TRUE);
            high = Math.max(high, item.value());
        }
        if (!(high > 0.0)) throw new IllegalArgumentException(
            "at least one circular bar must be positive");
        maximum = high;
        groups = List.copyOf(groupOrder.keySet());
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
            double centerX = width / 2.0;
            double centerY = top + (height - top) / 2.0;
            double outer = Math.min(width, height - top) * 0.34;
            double inner = outer * 0.22;
            double step = 2.0 * Math.PI / bars.size();
            graphics.setStroke(PlotSupport.strokeWidth(1f));
            graphics.setColor(DiagramSupport.GRID);
            for (int ring = 1; ring <= 4; ring++) {
                double radius = inner + (outer - inner) * ring / 4.0;
                graphics.draw(new Arc2D.Double(centerX - radius, centerY - radius,
                    2 * radius, 2 * radius, 0, 360, Arc2D.OPEN));
            }
            graphics.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 10));
            FontMetrics metrics = graphics.getFontMetrics();
            for (int index = 0; index < bars.size(); index++) {
                DiagramPlots.TreeItem bar = bars.get(index);
                double start = -Math.PI / 2.0 + index * step + step * 0.08;
                double end = start + step * 0.84;
                double radius = inner + (outer - inner) * bar.value() / maximum;
                Color color = DiagramSupport.color(options, groups.indexOf(bar.group()));
                Path2D path = sector(centerX, centerY, inner, radius, start, end);
                graphics.setColor(DiagramSupport.translucent(color, 205));
                graphics.fill(path);
                graphics.setColor(color.darker());
                graphics.setStroke(PlotSupport.stroke(options));
                graphics.draw(path);
                double middle = (start + end) / 2.0;
                String label = DiagramSupport.shorten(bar.label(), metrics, 75);
                double labelRadius = outer + 10;
                double x = centerX + labelRadius * Math.cos(middle);
                double y = centerY + labelRadius * Math.sin(middle);
                graphics.setColor(DiagramSupport.TEXT);
                graphics.drawString(label, (int) x
                    - (Math.cos(middle) < 0 ? metrics.stringWidth(label) : 0),
                    (int) y + metrics.getAscent() / 2);
            }
            DiagramSupport.legend(graphics, options, groups, 16, top + 15);
        } finally {
            graphics.dispose();
        }
    }

    private static Path2D sector(double centerX, double centerY, double inner,
            double outer, double start, double end) {
        Path2D path = new Path2D.Double();
        int segments = 12;
        path.moveTo(centerX + inner * Math.cos(start),
            centerY + inner * Math.sin(start));
        for (int index = 0; index <= segments; index++) {
            double angle = start + (end - start) * index / segments;
            path.lineTo(centerX + outer * Math.cos(angle),
                centerY + outer * Math.sin(angle));
        }
        for (int index = segments; index >= 0; index--) {
            double angle = start + (end - start) * index / segments;
            path.lineTo(centerX + inner * Math.cos(angle),
                centerY + inner * Math.sin(angle));
        }
        path.closePath();
        return path;
    }
}
