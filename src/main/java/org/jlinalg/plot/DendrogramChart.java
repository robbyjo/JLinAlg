/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.util.List;
import org.knowm.xchart.internal.chartpart.IChart;

final class DendrogramChart implements IChart {
    private final List<String> labels;
    private final HierarchicalClustering.Node tree;
    private final int[] order;
    private final int[] positions;
    private final DiagramPlots.Orientation orientation;
    private final PlotOptions options;

    DendrogramChart(DiagramPlots.ObservationData data,
            DiagramPlots.Orientation orientation, PlotOptions options) {
        ParallelChart.validate(data, 1);
        if (data.labels().size() < 2)
            throw new IllegalArgumentException("dendrogram needs at least two observations");
        labels = List.copyOf(data.labels());
        tree = HierarchicalClustering.cluster(ParallelChart.copy(data.values()));
        order = tree.order();
        positions = new int[order.length];
        for (int index = 0; index < order.length; index++) positions[order[index]] = index;
        this.orientation = orientation;
        this.options = options;
    }

    @Override public int getWidth() { return options.width(); }
    @Override public int getHeight() { return options.height(); }
    @Override public String getTitle() { return options.title(); }

    @Override
    public void paint(Graphics2D source, int width, int height) {
        Graphics2D graphics = DiagramSupport.begin(source, options, width, height);
        try {
            graphics.setColor(DiagramSupport.color(options, 0));
            graphics.setStroke(PlotSupport.stroke(options));
            graphics.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
            int top = DiagramSupport.top(options) + 15;
            if (orientation == DiagramPlots.Orientation.VERTICAL) {
                int left = 42, right = 25, bottom = 105;
                double plotWidth = width - left - right;
                double leafY = height - bottom;
                drawVertical(graphics, tree, left, plotWidth, leafY, top,
                    Math.max(tree.height(), 1e-12));
                FontMetrics metrics = graphics.getFontMetrics();
                AffineTransform original = graphics.getTransform();
                for (int position = 0; position < order.length; position++) {
                    double x = left + (position + 0.5) * plotWidth / order.length;
                    String label = DiagramSupport.shorten(labels.get(order[position]),
                        metrics, 85);
                    graphics.translate(x, leafY + 8);
                    graphics.rotate(-Math.PI / 3.0);
                    graphics.drawString(label, 0, metrics.getAscent());
                    graphics.setTransform(original);
                }
            } else {
                int left = 120, right = 32, bottom = 30;
                double plotHeight = height - top - bottom;
                double leafX = left;
                drawHorizontal(graphics, tree, top, plotHeight, leafX,
                    width - right, Math.max(tree.height(), 1e-12));
                FontMetrics metrics = graphics.getFontMetrics();
                for (int position = 0; position < order.length; position++) {
                    double y = top + (position + 0.5) * plotHeight / order.length;
                    String label = DiagramSupport.shorten(labels.get(order[position]),
                        metrics, left - 12);
                    graphics.drawString(label, left - 8 - metrics.stringWidth(label),
                        (int) y + metrics.getAscent() / 2 - 2);
                }
            }
        } finally {
            graphics.dispose();
        }
    }

    private double drawVertical(Graphics2D graphics,
            HierarchicalClustering.Node node, double left, double width,
            double leafY, double rootY, double maximumHeight) {
        if (node.isLeaf())
            return left + (positions[node.leaf()] + 0.5) * width / order.length;
        double firstX = drawVertical(graphics, node.left(), left, width,
            leafY, rootY, maximumHeight);
        double secondX = drawVertical(graphics, node.right(), left, width,
            leafY, rootY, maximumHeight);
        double y = leafY - (leafY - rootY) * node.height() / maximumHeight;
        double firstY = node.left().isLeaf() ? leafY
            : leafY - (leafY - rootY) * node.left().height() / maximumHeight;
        double secondY = node.right().isLeaf() ? leafY
            : leafY - (leafY - rootY) * node.right().height() / maximumHeight;
        graphics.drawLine((int) firstX, (int) firstY, (int) firstX, (int) y);
        graphics.drawLine((int) secondX, (int) secondY, (int) secondX, (int) y);
        graphics.drawLine((int) firstX, (int) y, (int) secondX, (int) y);
        return (firstX + secondX) / 2.0;
    }

    private double drawHorizontal(Graphics2D graphics,
            HierarchicalClustering.Node node, double top, double height,
            double leafX, double rootX, double maximumHeight) {
        if (node.isLeaf())
            return top + (positions[node.leaf()] + 0.5) * height / order.length;
        double firstY = drawHorizontal(graphics, node.left(), top, height,
            leafX, rootX, maximumHeight);
        double secondY = drawHorizontal(graphics, node.right(), top, height,
            leafX, rootX, maximumHeight);
        double x = leafX + (rootX - leafX) * node.height() / maximumHeight;
        double firstX = node.left().isLeaf() ? leafX
            : leafX + (rootX - leafX) * node.left().height() / maximumHeight;
        double secondX = node.right().isLeaf() ? leafX
            : leafX + (rootX - leafX) * node.right().height() / maximumHeight;
        graphics.drawLine((int) firstX, (int) firstY, (int) x, (int) firstY);
        graphics.drawLine((int) secondX, (int) secondY, (int) x, (int) secondY);
        graphics.drawLine((int) x, (int) firstY, (int) x, (int) secondY);
        return (firstY + secondY) / 2.0;
    }
}
