/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.geom.Arc2D;
import java.awt.geom.CubicCurve2D;
import java.awt.geom.QuadCurve2D;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.knowm.xchart.internal.chartpart.IChart;

final class NetworkDiagram implements IChart {
    enum Mode { NETWORK, CHORD, EDGE_BUNDLE }

    private final List<DiagramPlots.Edge> edges;
    private final List<String> nodes;
    private final Map<String, Integer> nodeIndex;
    private final List<String> groups;
    private final Mode mode;
    private final PlotOptions options;
    private final double maximumWeight;

    NetworkDiagram(List<DiagramPlots.Edge> source, Mode mode,
            PlotOptions options) {
        if (source == null || source.isEmpty())
            throw new IllegalArgumentException("diagram edges are required");
        edges = new ArrayList<>();
        nodeIndex = new LinkedHashMap<>();
        Map<String, Boolean> groupOrder = new LinkedHashMap<>();
        double maximum = 0.0;
        for (DiagramPlots.Edge edge : source) {
            if (edge == null || edge.source() == null || edge.source().isBlank()
                    || edge.target() == null || edge.target().isBlank()
                    || !Double.isFinite(edge.weight()) || !(edge.weight() > 0.0))
                throw new IllegalArgumentException(
                    "edges need source, target and a finite positive weight");
            String group = edge.group() == null || edge.group().isBlank()
                ? "edges" : edge.group();
            edges.add(new DiagramPlots.Edge(edge.source(), edge.target(),
                edge.weight(), group));
            nodeIndex.computeIfAbsent(edge.source(), ignored -> nodeIndex.size());
            nodeIndex.computeIfAbsent(edge.target(), ignored -> nodeIndex.size());
            groupOrder.put(group, Boolean.TRUE);
            maximum = Math.max(maximum, edge.weight());
        }
        nodes = List.copyOf(nodeIndex.keySet());
        groups = List.copyOf(groupOrder.keySet());
        maximumWeight = maximum;
        this.mode = mode;
        this.options = options;
    }

    @Override public int getWidth() { return options.width(); }
    @Override public int getHeight() { return options.height(); }
    @Override public String getTitle() { return options.title(); }

    @Override
    public void paint(Graphics2D source, int width, int height) {
        Graphics2D graphics = DiagramSupport.begin(source, options, width, height);
        try {
            if (mode == Mode.NETWORK) drawNetwork(graphics, width, height);
            else drawRadial(graphics, width, height);
        } finally {
            graphics.dispose();
        }
    }

    private void drawNetwork(Graphics2D graphics, int width, int height) {
        int top = DiagramSupport.top(options);
        int legend = groups.size() > 1 ? 120 : 15;
        double plotWidth = width - legend - 40.0;
        double plotHeight = height - top - 28.0;
        double[][] layout = forceLayout();
        double[] x = new double[nodes.size()], y = new double[nodes.size()];
        for (int index = 0; index < nodes.size(); index++) {
            x[index] = 20 + (layout[index][0] + 1.0) * plotWidth / 2.0;
            y[index] = top + 10 + (layout[index][1] + 1.0) * plotHeight / 2.0;
        }
        for (DiagramPlots.Edge edge : edges) {
            int first = nodeIndex.get(edge.source()), second = nodeIndex.get(edge.target());
            Color color = DiagramSupport.color(options, groups.indexOf(edge.group()));
            graphics.setColor(DiagramSupport.translucent(color, 125));
            graphics.setStroke(PlotSupport.stroke(options,
                scaledWidth(edge.weight())));
            graphics.drawLine((int) x[first], (int) y[first],
                (int) x[second], (int) y[second]);
        }
        graphics.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
        FontMetrics metrics = graphics.getFontMetrics();
        for (int index = 0; index < nodes.size(); index++) {
            Color color = DiagramSupport.color(options, index);
            DiagramSupport.marker(graphics, options.pointShape(), x[index], y[index],
                Math.max(7, options.markerSize() + 3), color);
            graphics.setColor(DiagramSupport.TEXT);
            graphics.drawString(nodes.get(index), (int) x[index] + 6,
                (int) y[index] - 5 + metrics.getAscent() / 2);
        }
        DiagramSupport.legend(graphics, options, groups,
            width - legend + 8, top + 15);
    }

    private void drawRadial(Graphics2D graphics, int width, int height) {
        int top = DiagramSupport.top(options);
        int legend = groups.size() > 1 ? 115 : 10;
        double centerX = (width - legend) / 2.0;
        double centerY = top + (height - top) / 2.0;
        double radius = Math.min(width - legend, height - top) * 0.36;
        double[] x = new double[nodes.size()], y = new double[nodes.size()];
        for (int index = 0; index < nodes.size(); index++) {
            double angle = angle(index);
            x[index] = centerX + radius * Math.cos(angle);
            y[index] = centerY + radius * Math.sin(angle);
        }
        if (mode == Mode.CHORD) drawNodeArcs(graphics, centerX, centerY, radius);
        for (DiagramPlots.Edge edge : edges) {
            int first = nodeIndex.get(edge.source()), second = nodeIndex.get(edge.target());
            Color color = mode == Mode.CHORD
                ? DiagramSupport.color(options, first)
                : DiagramSupport.color(options, groups.indexOf(edge.group()));
            graphics.setColor(DiagramSupport.translucent(color,
                mode == Mode.CHORD ? 105 : 115));
            graphics.setStroke(PlotSupport.stroke(options,
                scaledWidth(edge.weight())));
            if (mode == Mode.CHORD) {
                graphics.draw(new CubicCurve2D.Double(x[first], y[first], centerX,
                    centerY, centerX, centerY, x[second], y[second]));
            } else {
                int group = groups.indexOf(edge.group());
                double controlAngle = -Math.PI / 2.0
                    + 2.0 * Math.PI * group / Math.max(1, groups.size());
                double controlRadius = groups.size() == 1 ? 0.0 : radius * 0.18;
                double controlX = centerX + controlRadius * Math.cos(controlAngle);
                double controlY = centerY + controlRadius * Math.sin(controlAngle);
                graphics.draw(new QuadCurve2D.Double(x[first], y[first],
                    controlX, controlY, x[second], y[second]));
            }
        }
        drawRadialLabels(graphics, centerX, centerY, radius + 13);
        DiagramSupport.legend(graphics, options, groups,
            width - legend + 4, top + 14);
    }

    private void drawNodeArcs(Graphics2D graphics, double centerX,
            double centerY, double radius) {
        double extent = 360.0 / nodes.size() * 0.82;
        for (int index = 0; index < nodes.size(); index++) {
            double start = 90.0 - 360.0 * index / nodes.size() - extent / 2.0;
            graphics.setColor(DiagramSupport.color(options, index));
            graphics.setStroke(PlotSupport.strokeWidth(Math.max(5f,
                options.lineWidth() * 2.5f)));
            graphics.draw(new Arc2D.Double(centerX - radius, centerY - radius,
                2 * radius, 2 * radius, start, -extent, Arc2D.OPEN));
        }
    }

    private void drawRadialLabels(Graphics2D graphics, double centerX,
            double centerY, double radius) {
        graphics.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
        FontMetrics metrics = graphics.getFontMetrics();
        graphics.setColor(DiagramSupport.TEXT);
        for (int index = 0; index < nodes.size(); index++) {
            double angle = angle(index);
            String label = nodes.get(index);
            double x = centerX + radius * Math.cos(angle);
            double y = centerY + radius * Math.sin(angle);
            graphics.drawString(label, (int) x
                - (Math.cos(angle) < 0 ? metrics.stringWidth(label) : 0),
                (int) y + metrics.getAscent() / 2);
        }
    }

    private double angle(int index) {
        return -Math.PI / 2.0 + 2.0 * Math.PI * index / nodes.size();
    }

    private float scaledWidth(double weight) {
        return Math.max(0.6f, options.lineWidth()
            * (float) (0.55 + 1.7 * Math.sqrt(weight / maximumWeight)));
    }

    private double[][] forceLayout() {
        int count = nodes.size();
        double[][] positions = new double[count][2];
        for (int index = 0; index < count; index++) {
            positions[index][0] = 0.72 * Math.cos(2.0 * Math.PI * index / count);
            positions[index][1] = 0.72 * Math.sin(2.0 * Math.PI * index / count);
        }
        double ideal = Math.sqrt(4.0 / Math.max(1, count));
        for (int iteration = 0; iteration < 220; iteration++) {
            double[][] displacement = new double[count][2];
            for (int first = 0; first < count; first++)
                for (int second = first + 1; second < count; second++) {
                    double dx = positions[first][0] - positions[second][0];
                    double dy = positions[first][1] - positions[second][1];
                    double distance = Math.max(1e-4, Math.hypot(dx, dy));
                    double force = ideal * ideal / distance;
                    displacement[first][0] += dx / distance * force;
                    displacement[first][1] += dy / distance * force;
                    displacement[second][0] -= dx / distance * force;
                    displacement[second][1] -= dy / distance * force;
                }
            for (DiagramPlots.Edge edge : edges) {
                int first = nodeIndex.get(edge.source()), second = nodeIndex.get(edge.target());
                double dx = positions[first][0] - positions[second][0];
                double dy = positions[first][1] - positions[second][1];
                double distance = Math.max(1e-4, Math.hypot(dx, dy));
                double force = distance * distance / ideal
                    * Math.min(2.0, edge.weight() / maximumWeight + 0.5);
                displacement[first][0] -= dx / distance * force;
                displacement[first][1] -= dy / distance * force;
                displacement[second][0] += dx / distance * force;
                displacement[second][1] += dy / distance * force;
            }
            double temperature = 0.10 * (1.0 - iteration / 220.0);
            for (int index = 0; index < count; index++) {
                displacement[index][0] -= positions[index][0] * 0.18;
                displacement[index][1] -= positions[index][1] * 0.18;
                double length = Math.max(1e-9, Math.hypot(displacement[index][0],
                    displacement[index][1]));
                positions[index][0] += displacement[index][0] / length
                    * Math.min(length, temperature);
                positions[index][1] += displacement[index][1] / length
                    * Math.min(length, temperature);
                positions[index][0] = Math.max(-1.8,
                    Math.min(1.8, positions[index][0]));
                positions[index][1] = Math.max(-1.8,
                    Math.min(1.8, positions[index][1]));
            }
        }
        double meanX = 0.0, meanY = 0.0;
        for (double[] position : positions) {
            meanX += position[0];
            meanY += position[1];
        }
        meanX /= count;
        meanY /= count;
        double scale = 0.0;
        for (double[] position : positions)
            scale = Math.max(scale, Math.max(Math.abs(position[0] - meanX),
                Math.abs(position[1] - meanY)));
        scale = Math.max(scale, 1e-9) / 0.82;
        for (double[] position : positions) {
            position[0] = (position[0] - meanX) / scale;
            position[1] = (position[1] - meanY) / scale;
        }
        return positions;
    }
}
