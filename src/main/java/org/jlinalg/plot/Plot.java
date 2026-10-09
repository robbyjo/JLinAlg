/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.ChartEncoder;
import org.knowm.xchart.internal.chartpart.IChart;

/** A rendered plot that can be painted or exported to raster and vector formats. */
public final class Plot {
    private final IChart chart;

    Plot(IChart chart) {
        this.chart = chart;
    }

    public int width() { return chart.getWidth(); }
    public int height() { return chart.getHeight(); }
    public String title() { return chart.getTitle(); }

    /** Paints the plot into an existing Java2D surface. */
    public void paint(Graphics2D graphics, int width, int height) {
        if (graphics == null) throw new IllegalArgumentException("graphics surface is required");
        if (width < 1 || height < 1) throw new IllegalArgumentException("paint dimensions must be positive");
        chart.paint(graphics, width, height);
    }

    /** Produces an in-memory image at the plot's configured dimensions. */
    public BufferedImage image() {
        return BitmapEncoder.getBufferedImage(chart);
    }

    /** Encodes this plot in memory. */
    public byte[] render(PlotFormat format) throws IOException {
        if (format == null) throw new IllegalArgumentException("plot format is required");
        return ChartEncoder.getBytes(chart, format.extension());
    }

    /** Saves using the destination extension to select the output format. */
    public void save(Path destination) throws IOException {
        save(destination, PlotFormat.fromPath(destination));
    }

    /** Saves this plot without silently replacing an existing file. */
    public void save(Path destination, PlotFormat format) throws IOException {
        if (destination == null) throw new IllegalArgumentException("plot destination is required");
        Path absolute = destination.toAbsolutePath().normalize();
        if (Files.exists(absolute))
            throw new IOException("plot output exists; choose a new path: " + destination);
        Path parent = absolute.getParent();
        if (parent != null) Files.createDirectories(parent);
        Files.write(absolute, render(format));
    }
}
