/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.knowm.xchart.internal.chartpart.IChart;

final class TreemapChart implements IChart {
    private final Map<String, List<DiagramPlots.TreeItem>> groups;
    private final PlotOptions options;
    private final double total;

    TreemapChart(List<DiagramPlots.TreeItem> source, PlotOptions options) {
        if (source == null || source.isEmpty())
            throw new IllegalArgumentException("treemap items are required");
        groups = new LinkedHashMap<>();
        double sum = 0.0;
        for (DiagramPlots.TreeItem item : source) {
            if (item == null || item.label() == null || item.label().isBlank()
                    || !Double.isFinite(item.value()) || !(item.value() > 0.0))
                throw new IllegalArgumentException(
                    "treemap items need labels and finite positive values");
            String group = item.group() == null || item.group().isBlank()
                ? "data" : item.group();
            groups.computeIfAbsent(group, ignored -> new ArrayList<>()).add(
                new DiagramPlots.TreeItem(item.label(), group, item.value()));
            sum += item.value();
        }
        total = sum;
        this.options = options;
    }

    @Override public int getWidth() { return options.width(); }
    @Override public int getHeight() { return options.height(); }
    @Override public String getTitle() { return options.title(); }

    @Override
    public void paint(Graphics2D source, int width, int height) {
        Graphics2D graphics = DiagramSupport.begin(source, options, width, height);
        try {
            int top = DiagramSupport.top(options) + 5;
            Rect available = new Rect(18, top, width - 36, height - top - 18);
            double offset = available.x();
            int groupIndex = 0;
            for (Map.Entry<String, List<DiagramPlots.TreeItem>> entry : groups.entrySet()) {
                double groupTotal = entry.getValue().stream()
                    .mapToDouble(DiagramPlots.TreeItem::value).sum();
                double groupWidth = available.width() * groupTotal / total;
                Rect groupRect = new Rect(offset, available.y(), groupWidth,
                    available.height());
                drawGroup(graphics, groupRect, entry.getKey(), entry.getValue(),
                    groupIndex++);
                offset += groupWidth;
            }
        } finally {
            graphics.dispose();
        }
    }

    private void drawGroup(Graphics2D graphics, Rect rectangle, String group,
            List<DiagramPlots.TreeItem> items, int groupIndex) {
        Color base = DiagramSupport.color(options, groupIndex);
        graphics.setColor(base.darker());
        graphics.setStroke(PlotSupport.strokeWidth(2f));
        graphics.drawRect((int) rectangle.x(), (int) rectangle.y(),
            (int) rectangle.width(), (int) rectangle.height());
        double totalValue = items.stream().mapToDouble(DiagramPlots.TreeItem::value).sum();
        double offset = rectangle.height() >= rectangle.width()
            ? rectangle.y() : rectangle.x();
        for (int index = 0; index < items.size(); index++) {
            DiagramPlots.TreeItem item = items.get(index);
            double fraction = item.value() / totalValue;
            Rect itemRect;
            if (rectangle.height() >= rectangle.width()) {
                double itemHeight = rectangle.height() * fraction;
                itemRect = new Rect(rectangle.x(), offset, rectangle.width(), itemHeight);
                offset += itemHeight;
            } else {
                double itemWidth = rectangle.width() * fraction;
                itemRect = new Rect(offset, rectangle.y(), itemWidth, rectangle.height());
                offset += itemWidth;
            }
            drawItem(graphics, itemRect, item, base, index, items.size());
        }
        if (groups.size() > 1) {
            graphics.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 12));
            graphics.setColor(DiagramSupport.TEXT);
            graphics.drawString(group, (int) rectangle.x() + 5,
                (int) rectangle.y() + 15);
        }
    }

    private void drawItem(Graphics2D graphics, Rect rectangle,
            DiagramPlots.TreeItem item, Color base, int index, int count) {
        double light = 0.25 + 0.55 * index / Math.max(1.0, count - 1.0);
        Color fill = new Color(
            (int) Math.round(base.getRed() + (255 - base.getRed()) * light),
            (int) Math.round(base.getGreen() + (255 - base.getGreen()) * light),
            (int) Math.round(base.getBlue() + (255 - base.getBlue()) * light));
        graphics.setColor(fill);
        graphics.fillRect((int) rectangle.x() + 1, (int) rectangle.y() + 1,
            Math.max(1, (int) rectangle.width() - 2),
            Math.max(1, (int) rectangle.height() - 2));
        graphics.setColor(options.background());
        graphics.setStroke(PlotSupport.strokeWidth(1f));
        graphics.drawRect((int) rectangle.x(), (int) rectangle.y(),
            (int) rectangle.width(), (int) rectangle.height());
        if (rectangle.width() < 28 || rectangle.height() < 20) return;
        graphics.setColor(DiagramSupport.TEXT);
        graphics.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
        FontMetrics metrics = graphics.getFontMetrics();
        String label = DiagramSupport.shorten(item.label(), metrics,
            Math.max(10, (int) rectangle.width() - 8));
        graphics.drawString(label, (int) rectangle.x() + 4,
            (int) rectangle.y() + metrics.getAscent() + (groups.size() > 1 ? 16 : 2));
        if (rectangle.height() > 38) {
            String value = new java.text.DecimalFormat(options.numberPattern())
                .format(item.value());
            graphics.drawString(value, (int) rectangle.x() + 4,
                (int) rectangle.y() + 2 * metrics.getHeight()
                    + (groups.size() > 1 ? 16 : 2));
        }
    }

    private record Rect(double x, double y, double width, double height) { }
}
