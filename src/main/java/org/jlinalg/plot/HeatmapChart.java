/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.util.Arrays;
import java.util.List;
import org.knowm.xchart.internal.chartpart.IChart;

final class HeatmapChart implements IChart {
    private final List<String> columns;
    private final List<String> rows;
    private final double[][] values;
    private final PlotOptions options;
    private final HierarchicalClustering.Node rowTree;
    private final HierarchicalClustering.Node columnTree;
    private final int[] rowOrder;
    private final int[] columnOrder;
    private final double minimum;
    private final double maximum;

    HeatmapChart(DiagramPlots.HeatmapData data,
            DiagramPlots.ClusterMode clustering, PlotOptions options) {
        if (data == null || data.columns() == null || data.rows() == null
                || data.values() == null || data.rows().isEmpty()
                || data.columns().isEmpty())
            throw new IllegalArgumentException("heatmap grid is required");
        this.options = options;
        columns = List.copyOf(data.columns());
        rows = List.copyOf(data.rows());
        if (data.values().length != rows.size())
            throw new IllegalArgumentException("heatmap row count does not match labels");
        values = new double[rows.size()][columns.size()];
        double low = Double.POSITIVE_INFINITY, high = Double.NEGATIVE_INFINITY;
        for (int row = 0; row < values.length; row++) {
            if (data.values()[row] == null
                    || data.values()[row].length != columns.size())
                throw new IllegalArgumentException("heatmap grid must be complete and rectangular");
            for (int column = 0; column < columns.size(); column++) {
                double value = data.values()[row][column];
                if (!Double.isFinite(value))
                    throw new IllegalArgumentException("heatmap values must be finite");
                values[row][column] = value;
                low = Math.min(low, value);
                high = Math.max(high, value);
            }
        }
        boolean clusterRows = clustering == DiagramPlots.ClusterMode.ROWS
            || clustering == DiagramPlots.ClusterMode.BOTH;
        boolean clusterColumns = clustering == DiagramPlots.ClusterMode.COLUMNS
            || clustering == DiagramPlots.ClusterMode.BOTH;
        rowTree = clusterRows ? HierarchicalClustering.cluster(values) : null;
        columnTree = clusterColumns ? HierarchicalClustering.cluster(transpose(values)) : null;
        rowOrder = rowTree == null ? sequence(rows.size()) : rowTree.order();
        columnOrder = columnTree == null ? sequence(columns.size()) : columnTree.order();
        minimum = options.yMinimum() == null ? low : options.yMinimum();
        maximum = options.yMaximum() == null ? high : options.yMaximum();
        if (minimum > high || maximum < low)
            throw new IllegalArgumentException("heatmap color range excludes all values");
    }

    @Override public int getWidth() { return options.width(); }
    @Override public int getHeight() { return options.height(); }
    @Override public String getTitle() { return options.title(); }

    @Override
    public void paint(Graphics2D source, int width, int height) {
        Graphics2D graphics = DiagramSupport.begin(source, options, width, height);
        try {
            int top = DiagramSupport.top(options) + (columnTree == null ? 8 : 95);
            int left = rowTree == null ? 58 : 150;
            int right = 135;
            int bottom = 105;
            double gridWidth = Math.max(20.0, width - left - right);
            double gridHeight = Math.max(20.0, height - top - bottom);
            double cellWidth = gridWidth / columns.size();
            double cellHeight = gridHeight / rows.size();
            for (int displayedRow = 0; displayedRow < rows.size(); displayedRow++)
                for (int displayedColumn = 0; displayedColumn < columns.size(); displayedColumn++) {
                    double value = values[rowOrder[displayedRow]][columnOrder[displayedColumn]];
                    double fraction = maximum == minimum ? 0.5
                        : (value - minimum) / (maximum - minimum);
                    graphics.setColor(DiagramSupport.gradient(options, fraction));
                    int x1 = (int) Math.round(left + displayedColumn * cellWidth);
                    int x2 = (int) Math.round(left + (displayedColumn + 1) * cellWidth);
                    int y1 = (int) Math.round(top + displayedRow * cellHeight);
                    int y2 = (int) Math.round(top + (displayedRow + 1) * cellHeight);
                    graphics.fillRect(x1, y1, Math.max(1, x2 - x1),
                        Math.max(1, y2 - y1));
                }
            graphics.setColor(new Color(0x55, 0x55, 0x55));
            graphics.drawRect(left, top, (int) gridWidth, (int) gridHeight);
            drawLabels(graphics, left, top, gridWidth, gridHeight,
                cellWidth, cellHeight);
            if (rowTree != null) drawRowTree(graphics, rowTree, left - 8,
                top, gridHeight, cellHeight);
            if (columnTree != null) drawColumnTree(graphics, columnTree,
                left, top - 8, gridWidth, cellWidth);
            drawScale(graphics, left + (int) gridWidth + 72, top,
                Math.min(180, (int) gridHeight));
        } finally {
            graphics.dispose();
        }
    }

    private void drawLabels(Graphics2D graphics, int left, int top,
            double gridWidth, double gridHeight, double cellWidth,
            double cellHeight) {
        graphics.setColor(DiagramSupport.TEXT);
        graphics.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
        FontMetrics metrics = graphics.getFontMetrics();
        for (int displayedRow = 0; displayedRow < rows.size(); displayedRow++) {
            String label = DiagramSupport.shorten(rows.get(rowOrder[displayedRow]),
                metrics, 65);
            graphics.drawString(label, left + (int) gridWidth + 7,
                (int) (top + (displayedRow + 0.5) * cellHeight
                    + metrics.getAscent() / 2.0 - 2));
        }
        AffineTransform original = graphics.getTransform();
        for (int displayedColumn = 0; displayedColumn < columns.size(); displayedColumn++) {
            String label = DiagramSupport.shorten(
                columns.get(columnOrder[displayedColumn]), metrics, 85);
            double x = left + (displayedColumn + 0.5) * cellWidth;
            double y = top + gridHeight + 7;
            graphics.translate(x, y);
            graphics.rotate(-Math.PI / 3.0);
            graphics.drawString(label, 0, metrics.getAscent());
            graphics.setTransform(original);
        }
        graphics.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 13));
        metrics = graphics.getFontMetrics();
        graphics.drawString(options.xLabel(),
            left + (int) (gridWidth - metrics.stringWidth(options.xLabel())) / 2,
            top + (int) gridHeight + 88);
        original = graphics.getTransform();
        graphics.rotate(-Math.PI / 2.0);
        graphics.drawString(options.yLabel(),
            -(top + (int) (gridHeight + metrics.stringWidth(options.yLabel())) / 2), 18);
        graphics.setTransform(original);
    }

    private void drawScale(Graphics2D graphics, int x, int y, int height) {
        int barWidth = 16;
        for (int offset = 0; offset < height; offset++) {
            graphics.setColor(DiagramSupport.gradient(options,
                1.0 - offset / Math.max(1.0, height - 1.0)));
            graphics.drawLine(x, y + offset, x + barWidth, y + offset);
        }
        graphics.setColor(DiagramSupport.TEXT);
        graphics.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
        java.text.DecimalFormat format = new java.text.DecimalFormat(options.numberPattern());
        graphics.drawString(format.format(maximum), x + barWidth + 5, y + 5);
        graphics.drawString(format.format(minimum), x + barWidth + 5, y + height);
        graphics.drawRect(x, y, barWidth, height);
    }

    private void drawRowTree(Graphics2D graphics, HierarchicalClustering.Node tree,
            double leafX, double top, double height, double cellHeight) {
        int[] positions = inverse(rowOrder);
        graphics.setColor(DiagramSupport.TEXT);
        graphics.setStroke(PlotSupport.stroke(options));
        drawRowNode(graphics, tree, leafX, leafX - 82, top, cellHeight,
            positions, Math.max(tree.height(), 1e-12));
    }

    private double drawRowNode(Graphics2D graphics,
            HierarchicalClustering.Node node, double leafX, double rootX,
            double top, double cellHeight, int[] positions, double maximumHeight) {
        if (node.isLeaf()) return top + (positions[node.leaf()] + 0.5) * cellHeight;
        double firstY = drawRowNode(graphics, node.left(), leafX, rootX,
            top, cellHeight, positions, maximumHeight);
        double secondY = drawRowNode(graphics, node.right(), leafX, rootX,
            top, cellHeight, positions, maximumHeight);
        double x = leafX - (leafX - rootX) * node.height() / maximumHeight;
        double firstX = node.left().isLeaf() ? leafX
            : leafX - (leafX - rootX) * node.left().height() / maximumHeight;
        double secondX = node.right().isLeaf() ? leafX
            : leafX - (leafX - rootX) * node.right().height() / maximumHeight;
        graphics.drawLine((int) firstX, (int) firstY, (int) x, (int) firstY);
        graphics.drawLine((int) secondX, (int) secondY, (int) x, (int) secondY);
        graphics.drawLine((int) x, (int) firstY, (int) x, (int) secondY);
        return (firstY + secondY) / 2.0;
    }

    private void drawColumnTree(Graphics2D graphics,
            HierarchicalClustering.Node tree, double left, double leafY,
            double width, double cellWidth) {
        int[] positions = inverse(columnOrder);
        graphics.setColor(DiagramSupport.TEXT);
        graphics.setStroke(PlotSupport.stroke(options));
        drawColumnNode(graphics, tree, leafY, leafY - 72, left, cellWidth,
            positions, Math.max(tree.height(), 1e-12));
    }

    private double drawColumnNode(Graphics2D graphics,
            HierarchicalClustering.Node node, double leafY, double rootY,
            double left, double cellWidth, int[] positions, double maximumHeight) {
        if (node.isLeaf()) return left + (positions[node.leaf()] + 0.5) * cellWidth;
        double firstX = drawColumnNode(graphics, node.left(), leafY, rootY,
            left, cellWidth, positions, maximumHeight);
        double secondX = drawColumnNode(graphics, node.right(), leafY, rootY,
            left, cellWidth, positions, maximumHeight);
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

    private static double[][] transpose(double[][] source) {
        double[][] result = new double[source[0].length][source.length];
        for (int row = 0; row < source.length; row++)
            for (int column = 0; column < source[row].length; column++)
                result[column][row] = source[row][column];
        return result;
    }

    private static int[] sequence(int length) {
        int[] result = new int[length];
        for (int index = 0; index < length; index++) result[index] = index;
        return result;
    }

    private static int[] inverse(int[] order) {
        int[] result = new int[order.length];
        for (int position = 0; position < order.length; position++)
            result[order[position]] = position;
        return result;
    }
}
