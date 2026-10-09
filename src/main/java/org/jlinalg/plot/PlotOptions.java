/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.awt.Color;
import java.util.Arrays;
import java.util.Objects;

/** Shared dimensions, labels, scales, colors and theme for a plot. */
public final class PlotOptions {
    private static final Color[] ACCESSIBLE_PALETTE = {
        new Color(0x00, 0x72, 0xB2), new Color(0xD5, 0x5E, 0x00),
        new Color(0x00, 0x9E, 0x73), new Color(0xCC, 0x79, 0xA7),
        new Color(0xE6, 0x9F, 0x00), new Color(0x56, 0xB4, 0xE9),
        new Color(0xF0, 0xE4, 0x42), new Color(0x00, 0x00, 0x00)
    };

    private final int width;
    private final int height;
    private final String title;
    private final String xLabel;
    private final String yLabel;
    private final PlotTheme theme;
    private final LegendPosition legendPosition;
    private final boolean gridVisible;
    private final boolean xLogarithmic;
    private final boolean yLogarithmic;
    private final Double xMinimum;
    private final Double xMaximum;
    private final Double yMinimum;
    private final Double yMaximum;
    private final String numberPattern;
    private final int markerSize;
    private final float lineWidth;
    private final LinePattern linePattern;
    private final PointShape pointShape;
    private final Color background;
    private final Color[] colors;

    private PlotOptions(Builder builder) {
        width = builder.width;
        height = builder.height;
        title = builder.title;
        xLabel = builder.xLabel;
        yLabel = builder.yLabel;
        theme = builder.theme;
        legendPosition = builder.legendPosition;
        gridVisible = builder.gridVisible;
        xLogarithmic = builder.xLogarithmic;
        yLogarithmic = builder.yLogarithmic;
        xMinimum = builder.xMinimum;
        xMaximum = builder.xMaximum;
        yMinimum = builder.yMinimum;
        yMaximum = builder.yMaximum;
        numberPattern = builder.numberPattern;
        markerSize = builder.markerSize;
        lineWidth = builder.lineWidth;
        linePattern = builder.linePattern;
        pointShape = builder.pointShape;
        background = builder.background;
        colors = builder.colors.clone();
    }

    public static Builder builder() { return new Builder(); }

    public static PlotOptions defaults() { return builder().build(); }

    public int width() { return width; }
    public int height() { return height; }
    public String title() { return title; }
    public String xLabel() { return xLabel; }
    public String yLabel() { return yLabel; }
    public PlotTheme theme() { return theme; }
    public LegendPosition legendPosition() { return legendPosition; }
    public boolean gridVisible() { return gridVisible; }
    public boolean xLogarithmic() { return xLogarithmic; }
    public boolean yLogarithmic() { return yLogarithmic; }
    public Double xMinimum() { return xMinimum; }
    public Double xMaximum() { return xMaximum; }
    public Double yMinimum() { return yMinimum; }
    public Double yMaximum() { return yMaximum; }
    public String numberPattern() { return numberPattern; }
    public int markerSize() { return markerSize; }
    public float lineWidth() { return lineWidth; }
    public LinePattern linePattern() { return linePattern; }
    public PointShape pointShape() { return pointShape; }
    public Color background() { return background; }
    public Color[] colors() { return colors.clone(); }

    /** Fluent builder for immutable plot options. */
    public static final class Builder {
        private int width = 800;
        private int height = 600;
        private String title = "";
        private String xLabel = "";
        private String yLabel = "";
        private PlotTheme theme = PlotTheme.GGPLOT2;
        private LegendPosition legendPosition = LegendPosition.AUTO;
        private boolean gridVisible = true;
        private boolean xLogarithmic;
        private boolean yLogarithmic;
        private Double xMinimum;
        private Double xMaximum;
        private Double yMinimum;
        private Double yMaximum;
        private String numberPattern = "0.####";
        private int markerSize = 8;
        private float lineWidth = 2.0f;
        private LinePattern linePattern = LinePattern.SOLID;
        private PointShape pointShape = PointShape.CIRCLE;
        private Color background = Color.WHITE;
        private Color[] colors = ACCESSIBLE_PALETTE.clone();

        private Builder() { }

        public Builder size(int valueWidth, int valueHeight) {
            width = valueWidth;
            height = valueHeight;
            return this;
        }

        public Builder title(String value) { title = Objects.requireNonNull(value); return this; }
        public Builder xLabel(String value) { xLabel = Objects.requireNonNull(value); return this; }
        public Builder yLabel(String value) { yLabel = Objects.requireNonNull(value); return this; }
        public Builder theme(PlotTheme value) { theme = Objects.requireNonNull(value); return this; }
        public Builder legend(LegendPosition value) { legendPosition = Objects.requireNonNull(value); return this; }
        public Builder grid(boolean value) { gridVisible = value; return this; }
        public Builder logarithmicX(boolean value) { xLogarithmic = value; return this; }
        public Builder logarithmicY(boolean value) { yLogarithmic = value; return this; }

        public Builder xRange(double minimum, double maximum) {
            xMinimum = minimum;
            xMaximum = maximum;
            return this;
        }

        public Builder yRange(double minimum, double maximum) {
            yMinimum = minimum;
            yMaximum = maximum;
            return this;
        }

        public Builder numberPattern(String value) {
            numberPattern = Objects.requireNonNull(value);
            return this;
        }

        public Builder markerSize(int value) { markerSize = value; return this; }
        public Builder lineWidth(float value) { lineWidth = value; return this; }
        public Builder linePattern(LinePattern value) {
            linePattern = Objects.requireNonNull(value);
            return this;
        }
        public Builder pointShape(PointShape value) {
            pointShape = Objects.requireNonNull(value);
            return this;
        }
        public Builder background(Color value) { background = Objects.requireNonNull(value); return this; }

        public Builder colors(Color... values) {
            Objects.requireNonNull(values);
            colors = Arrays.stream(values).map(Objects::requireNonNull).toArray(Color[]::new);
            return this;
        }

        public PlotOptions build() {
            if (width < 100 || height < 100)
                throw new IllegalArgumentException("plot dimensions must each be at least 100 pixels");
            if (markerSize < 1) throw new IllegalArgumentException("marker size must be positive");
            if (!(lineWidth > 0.0f) || !Float.isFinite(lineWidth))
                throw new IllegalArgumentException("line width must be finite and positive");
            if (colors.length == 0) throw new IllegalArgumentException("at least one plot color is required");
            checkRange("x", xMinimum, xMaximum, xLogarithmic);
            checkRange("y", yMinimum, yMaximum, yLogarithmic);
            return new PlotOptions(this);
        }

        private static void checkRange(String name, Double minimum, Double maximum,
                boolean logarithmic) {
            if ((minimum == null) != (maximum == null))
                throw new IllegalArgumentException(name + " range needs both bounds");
            if (minimum != null && (!Double.isFinite(minimum) || !Double.isFinite(maximum)
                    || minimum >= maximum))
                throw new IllegalArgumentException(name + " range must be finite and increasing");
            if (logarithmic && minimum != null && minimum <= 0.0)
                throw new IllegalArgumentException(name + " logarithmic range must be positive");
        }
    }
}
