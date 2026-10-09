/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.awt.Color;
import java.awt.BasicStroke;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.geom.CubicCurve2D;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.knowm.xchart.internal.chartpart.IChart;

final class SankeyChart implements IChart {
    private final List<DiagramPlots.Edge> edges;
    private final List<String> nodes;
    private final Map<String, Integer> nodeIndex;
    private final List<String> groups;
    private final int[] depths;
    private final int maximumDepth;
    private final double[] nodeValues;
    private final PlotOptions options;

    SankeyChart(List<DiagramPlots.Edge> source, PlotOptions options) {
        if (source == null || source.isEmpty())
            throw new IllegalArgumentException("Sankey edges are required");
        edges = new ArrayList<>();
        nodeIndex = new LinkedHashMap<>();
        Map<String, Boolean> groupOrder = new LinkedHashMap<>();
        for (DiagramPlots.Edge edge : source) {
            if (edge == null || edge.source() == null || edge.source().isBlank()
                    || edge.target() == null || edge.target().isBlank()
                    || edge.source().equals(edge.target()) || !Double.isFinite(edge.weight())
                    || !(edge.weight() > 0.0))
                throw new IllegalArgumentException(
                    "Sankey edges need distinct nodes and finite positive weights");
            String group = edge.group() == null || edge.group().isBlank()
                ? "flow" : edge.group();
            edges.add(new DiagramPlots.Edge(edge.source(), edge.target(),
                edge.weight(), group));
            nodeIndex.computeIfAbsent(edge.source(), ignored -> nodeIndex.size());
            nodeIndex.computeIfAbsent(edge.target(), ignored -> nodeIndex.size());
            groupOrder.put(group, Boolean.TRUE);
        }
        nodes = List.copyOf(nodeIndex.keySet());
        groups = List.copyOf(groupOrder.keySet());
        depths = depths();
        maximumDepth = Arrays.stream(depths).max().orElse(0);
        nodeValues = values();
        this.options = options;
    }

    @Override public int getWidth() { return options.width(); }
    @Override public int getHeight() { return options.height(); }
    @Override public String getTitle() { return options.title(); }

    @Override
    public void paint(Graphics2D source, int width, int height) {
        Graphics2D graphics = DiagramSupport.begin(source, options, width, height);
        try {
            int top = DiagramSupport.top(options) + 12;
            int left = 45, right = groups.size() > 1 ? 180 : 85, bottom = 28;
            double plotWidth = width - left - right;
            double plotHeight = height - top - bottom;
            List<List<Integer>> layers = new ArrayList<>();
            for (int depth = 0; depth <= maximumDepth; depth++) layers.add(new ArrayList<>());
            for (int node = 0; node < nodes.size(); node++) layers.get(depths[node]).add(node);
            double maximumLayerValue = layers.stream().mapToDouble(layer -> layer.stream()
                .mapToDouble(node -> nodeValues[node]).sum()).max().orElse(1.0);
            int maximumLayerNodes = layers.stream().mapToInt(List::size).max().orElse(1);
            double gap = 14.0;
            double scale = Math.max(0.5,
                (plotHeight - gap * Math.max(0, maximumLayerNodes - 1))
                    / maximumLayerValue);
            double[] x = new double[nodes.size()], y = new double[nodes.size()];
            double nodeWidth = 16.0;
            for (int depth = 0; depth < layers.size(); depth++) {
                List<Integer> layer = layers.get(depth);
                double occupied = layer.stream().mapToDouble(node -> nodeValues[node] * scale)
                    .sum() + gap * Math.max(0, layer.size() - 1);
                double offset = top + (plotHeight - occupied) / 2.0;
                for (int node : layer) {
                    x[node] = left + (maximumDepth == 0 ? plotWidth / 2.0
                        : plotWidth * depth / maximumDepth);
                    y[node] = offset;
                    offset += nodeValues[node] * scale + gap;
                }
            }
            double[] sourceOffset = new double[nodes.size()];
            double[] targetOffset = new double[nodes.size()];
            for (DiagramPlots.Edge edge : edges) {
                int first = nodeIndex.get(edge.source()), second = nodeIndex.get(edge.target());
                double thickness = Math.max(1.0, edge.weight() * scale);
                double startY = y[first] + sourceOffset[first] + thickness / 2.0;
                double endY = y[second] + targetOffset[second] + thickness / 2.0;
                sourceOffset[first] += thickness;
                targetOffset[second] += thickness;
                double startX = x[first] + nodeWidth;
                double endX = x[second];
                double control = (startX + endX) / 2.0;
                Color color = DiagramSupport.color(options, groups.indexOf(edge.group()));
                graphics.setColor(DiagramSupport.translucent(color, 115));
                graphics.setStroke(flowStroke((float) thickness));
                graphics.draw(new CubicCurve2D.Double(startX, startY, control,
                    startY, control, endY, endX, endY));
            }
            graphics.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
            FontMetrics metrics = graphics.getFontMetrics();
            for (int node = 0; node < nodes.size(); node++) {
                Color color = DiagramSupport.color(options, node);
                graphics.setColor(color.darker());
                graphics.fillRoundRect((int) x[node], (int) y[node], (int) nodeWidth,
                    Math.max(2, (int) Math.round(nodeValues[node] * scale)), 3, 3);
                graphics.setColor(DiagramSupport.TEXT);
                int labelX = depths[node] == maximumDepth
                    ? (int) x[node] - metrics.stringWidth(nodes.get(node)) - 6
                    : (int) x[node] + 20;
                graphics.drawString(nodes.get(node), labelX,
                    (int) (y[node] + Math.min(14.0,
                        nodeValues[node] * scale / 2.0 + 4.0)));
            }
            DiagramSupport.legend(graphics, options, groups,
                width - right + 42, top + 10);
        } finally {
            graphics.dispose();
        }
    }

    private int[] depths() {
        int count = nodes.size();
        int[] indegree = new int[count], depth = new int[count];
        List<List<Integer>> outgoing = new ArrayList<>();
        for (int node = 0; node < count; node++) outgoing.add(new ArrayList<>());
        for (DiagramPlots.Edge edge : edges) {
            int first = nodeIndex.get(edge.source()), second = nodeIndex.get(edge.target());
            outgoing.get(first).add(second);
            indegree[second]++;
        }
        ArrayDeque<Integer> queue = new ArrayDeque<>();
        for (int node = 0; node < count; node++) if (indegree[node] == 0) queue.add(node);
        int processed = 0;
        while (!queue.isEmpty()) {
            int node = queue.remove();
            processed++;
            for (int target : outgoing.get(node)) {
                depth[target] = Math.max(depth[target], depth[node] + 1);
                if (--indegree[target] == 0) queue.add(target);
            }
        }
        if (processed != count)
            throw new IllegalArgumentException("Sankey input must be a directed acyclic graph");
        return depth;
    }

    private BasicStroke flowStroke(float width) {
        BasicStroke base = PlotSupport.stroke(options, width);
        return new BasicStroke(width, BasicStroke.CAP_BUTT,
            BasicStroke.JOIN_ROUND, base.getMiterLimit(), base.getDashArray(),
            base.getDashPhase());
    }

    private double[] values() {
        double[] incoming = new double[nodes.size()], outgoing = new double[nodes.size()];
        for (DiagramPlots.Edge edge : edges) {
            outgoing[nodeIndex.get(edge.source())] += edge.weight();
            incoming[nodeIndex.get(edge.target())] += edge.weight();
        }
        double[] result = new double[nodes.size()];
        for (int node = 0; node < result.length; node++)
            result[node] = Math.max(incoming[node], outgoing[node]);
        return result;
    }
}
