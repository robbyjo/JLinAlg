/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.markers.SeriesMarkers;

/** Builder for layered Cartesian point, line, step and area plots. */
public final class XYPlotBuilder {
    private final List<Layer> layers = new ArrayList<>();
    private PlotOptions options = PlotOptions.defaults();

    public XYPlotBuilder options(PlotOptions value) {
        options = Objects.requireNonNull(value);
        return this;
    }

    public XYPlotBuilder points(String name, double[] x, double[] y) {
        return add(name, x, y, Geometry.POINTS, false, null);
    }

    public XYPlotBuilder points(String name, double[] x, double[] y, Color color) {
        return add(name, x, y, Geometry.POINTS, false, color);
    }

    public XYPlotBuilder line(String name, double[] x, double[] y) {
        return add(name, x, y, Geometry.LINE, false, null);
    }

    public XYPlotBuilder line(String name, double[] x, double[] y, Color color) {
        return add(name, x, y, Geometry.LINE, false, color);
    }

    public XYPlotBuilder smoothLine(String name, double[] x, double[] y) {
        return add(name, x, y, Geometry.LINE, true, null);
    }

    public XYPlotBuilder lineAndPoints(String name, double[] x, double[] y) {
        return add(name, x, y, Geometry.LINE_POINTS, false, null);
    }

    public XYPlotBuilder area(String name, double[] x, double[] y) {
        return add(name, x, y, Geometry.AREA, false, null);
    }

    public XYPlotBuilder step(String name, double[] x, double[] y) {
        return add(name, x, y, Geometry.STEP, false, null);
    }

    private XYPlotBuilder add(String name, double[] xValues, double[] yValues,
            Geometry geometry, boolean smooth, Color color) {
        if (name == null || name.isBlank())
            throw new IllegalArgumentException("series name is required");
        double[] x = PlotSupport.finite(xValues, "x");
        double[] y = PlotSupport.finite(yValues, "y");
        PlotSupport.sameLength(x, y, "x and y");
        layers.add(new Layer(name, x, y, geometry, smooth, color));
        return this;
    }

    public Plot build() {
        if (layers.isEmpty()) throw new IllegalStateException("at least one plot layer is required");
        XYChart chart = new XYChartBuilder().width(options.width())
            .height(options.height()).title(options.title())
            .xAxisTitle(options.xLabel()).yAxisTitle(options.yLabel())
            .theme(PlotSupport.theme(options)).build();
        for (Layer layer : layers) {
            XYSeries series = chart.addSeries(layer.name(), layer.x(), layer.y());
            series.setXYSeriesRenderStyle(switch (layer.geometry()) {
                case POINTS -> XYSeries.XYSeriesRenderStyle.Scatter;
                case LINE, LINE_POINTS -> XYSeries.XYSeriesRenderStyle.Line;
                case AREA -> XYSeries.XYSeriesRenderStyle.Area;
                case STEP -> XYSeries.XYSeriesRenderStyle.Step;
            });
            if (layer.geometry() == Geometry.LINE || layer.geometry() == Geometry.AREA
                    || layer.geometry() == Geometry.STEP)
                series.setMarker(SeriesMarkers.NONE);
            else series.setMarker(PlotSupport.marker(options.pointShape()));
            series.setLineStyle(PlotSupport.stroke(options));
            if (layer.color() != null) {
                series.setLineColor(layer.color());
                series.setMarkerColor(layer.color());
                series.setFillColor(new Color(layer.color().getRed(), layer.color().getGreen(),
                    layer.color().getBlue(), 120));
            }
            series.setSmooth(layer.smooth());
        }
        PlotSupport.style(chart, options, layers.size() > 1);
        return new Plot(chart);
    }

    private enum Geometry { POINTS, LINE, LINE_POINTS, AREA, STEP }

    private record Layer(String name, double[] x, double[] y, Geometry geometry,
            boolean smooth, Color color) { }
}
