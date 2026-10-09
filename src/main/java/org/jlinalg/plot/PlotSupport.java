/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.awt.Color;
import java.awt.BasicStroke;
import org.knowm.xchart.internal.chartpart.AxesChart;
import org.knowm.xchart.internal.chartpart.Chart;
import org.knowm.xchart.style.AxesChartStyler;
import org.knowm.xchart.style.Styler;
import org.knowm.xchart.style.markers.Marker;
import org.knowm.xchart.style.markers.SeriesMarkers;

final class PlotSupport {
    private PlotSupport() { }

    static Styler.ChartTheme theme(PlotOptions options) {
        return options.theme() == PlotTheme.GGPLOT2
            ? Styler.ChartTheme.GGPlot2 : Styler.ChartTheme.XChart;
    }

    static void style(Chart<?, ?> chart, PlotOptions options,
            boolean defaultLegend) {
        Styler styler = chart.getStyler();
        styler.setChartBackgroundColor(options.background());
        styler.setSeriesColors(options.colors());
        BasicStroke stroke = stroke(options);
        BasicStroke[] lines = new BasicStroke[options.colors().length];
        java.util.Arrays.fill(lines, stroke);
        styler.setSeriesLines(lines);
        Marker marker = marker(options.pointShape());
        Marker[] markers = new Marker[options.colors().length];
        java.util.Arrays.fill(markers, marker);
        styler.setSeriesMarkers(markers);
        styler.setMarkerSize(options.markerSize());
        styler.setChartTitleVisible(!options.title().isBlank());
        styler.setDecimalPattern(options.numberPattern());
        LegendPosition position = options.legendPosition();
        styler.setLegendVisible(position == LegendPosition.AUTO
            ? defaultLegend : position != LegendPosition.NONE);
        if (position != LegendPosition.AUTO && position != LegendPosition.NONE)
            styler.setLegendPosition(legend(position));
        if (options.theme() == PlotTheme.MINIMAL || options.theme() == PlotTheme.CLASSIC) {
            styler.setPlotBackgroundColor(options.background());
            styler.setPlotBorderVisible(false);
        }
        if (chart instanceof AxesChart<?, ?> axes) styleAxes(axes, options);
    }

    private static void styleAxes(AxesChart<?, ?> chart, PlotOptions options) {
        AxesChartStyler styler = chart.getStyler();
        styler.setPlotGridLinesVisible(options.gridVisible()
            && options.theme() != PlotTheme.CLASSIC);
        if (options.theme() == PlotTheme.MINIMAL) {
            styler.setPlotGridLinesColor(new Color(0xE5, 0xE7, 0xEB));
            styler.setAxisTicksLineVisible(false);
        }
        styler.setXAxisLogarithmic(options.xLogarithmic());
        styler.setYAxisLogarithmic(options.yLogarithmic());
        if (options.xMinimum() != null) {
            styler.setXAxisMin(options.xMinimum());
            styler.setXAxisMax(options.xMaximum());
        }
        if (options.yMinimum() != null) {
            styler.setYAxisMin(options.yMinimum());
            styler.setYAxisMax(options.yMaximum());
        }
        styler.setXAxisDecimalPattern(options.numberPattern());
        styler.setYAxisDecimalPattern(options.numberPattern());
    }

    private static Styler.LegendPosition legend(LegendPosition position) {
        return switch (position) {
            case RIGHT -> Styler.LegendPosition.OutsideE;
            case BOTTOM -> Styler.LegendPosition.OutsideS;
            case INSIDE_NORTH_WEST -> Styler.LegendPosition.InsideNW;
            case INSIDE_NORTH_EAST -> Styler.LegendPosition.InsideNE;
            case INSIDE_SOUTH_WEST -> Styler.LegendPosition.InsideSW;
            case INSIDE_SOUTH_EAST -> Styler.LegendPosition.InsideSE;
            case AUTO, NONE -> throw new IllegalArgumentException(
                "legend position has no backend placement: " + position);
        };
    }

    static double[] finite(double[] values, String name) {
        if (values == null || values.length == 0)
            throw new IllegalArgumentException(name + " values are required");
        double[] copy = values.clone();
        for (double value : copy)
            if (!Double.isFinite(value))
                throw new IllegalArgumentException(name + " values must be finite");
        return copy;
    }

    static BasicStroke stroke(PlotOptions options) {
        return stroke(options, options.lineWidth());
    }

    static BasicStroke stroke(PlotOptions options, float width) {
        float[] dash = switch (options.linePattern()) {
            case SOLID -> null;
            case DASHED -> new float[] {6f * width, 4f * width};
            case DOTTED -> new float[] {width, 3f * width};
            case DASH_DOT -> new float[] {6f * width, 3f * width,
                width, 3f * width};
        };
        return dash == null ? new BasicStroke(width, BasicStroke.CAP_ROUND,
            BasicStroke.JOIN_ROUND) : new BasicStroke(width,
            BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 10f, dash, 0f);
    }

    static BasicStroke strokeWidth(float width) {
        return new BasicStroke(width, BasicStroke.CAP_ROUND,
            BasicStroke.JOIN_ROUND);
    }

    static Marker marker(PointShape shape) {
        return switch (shape) {
            case CIRCLE -> SeriesMarkers.CIRCLE;
            case SQUARE -> SeriesMarkers.SQUARE;
            case DIAMOND -> SeriesMarkers.DIAMOND;
            case TRIANGLE_UP -> SeriesMarkers.TRIANGLE_UP;
            case TRIANGLE_DOWN -> SeriesMarkers.TRIANGLE_DOWN;
            case CROSS -> SeriesMarkers.CROSS;
            case PLUS -> SeriesMarkers.PLUS;
            case NONE -> SeriesMarkers.NONE;
        };
    }

    static void sameLength(double[] first, double[] second, String names) {
        if (first.length != second.length)
            throw new IllegalArgumentException(names + " must have equal lengths");
    }
}
