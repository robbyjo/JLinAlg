/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Path2D;
import java.util.List;

final class DiagramSupport {
    static final Color TEXT = new Color(0x2B, 0x2B, 0x2B);
    static final Color GRID = new Color(0xD1, 0xD5, 0xDB);

    private DiagramSupport() { }

    static Graphics2D begin(Graphics2D source, PlotOptions options,
            int width, int height) {
        Graphics2D graphics = (Graphics2D) source.create();
        graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
            RenderingHints.VALUE_ANTIALIAS_ON);
        graphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,
            RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        graphics.setColor(options.background());
        graphics.fillRect(0, 0, width, height);
        if (!options.title().isBlank()) {
            graphics.setColor(TEXT);
            graphics.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
            FontMetrics metrics = graphics.getFontMetrics();
            graphics.drawString(options.title(),
                (width - metrics.stringWidth(options.title())) / 2, 29);
        }
        return graphics;
    }

    static int top(PlotOptions options) { return options.title().isBlank() ? 24 : 52; }

    static Color color(PlotOptions options, int index) {
        Color[] colors = options.colors();
        return colors[Math.floorMod(index, colors.length)];
    }

    static Color translucent(Color color, int alpha) {
        return new Color(color.getRed(), color.getGreen(), color.getBlue(), alpha);
    }

    static Color gradient(PlotOptions options, double fraction) {
        Color[] colors = options.colors();
        double clipped = Math.max(0.0, Math.min(1.0, fraction));
        if (colors.length == 1) return colors[0];
        double position = clipped * (colors.length - 1);
        int low = Math.min((int) Math.floor(position), colors.length - 2);
        return interpolate(colors[low], colors[low + 1], position - low);
    }

    private static Color interpolate(Color first, Color second, double fraction) {
        return new Color(
            channel(first.getRed(), second.getRed(), fraction),
            channel(first.getGreen(), second.getGreen(), fraction),
            channel(first.getBlue(), second.getBlue(), fraction));
    }

    private static int channel(int first, int second, double fraction) {
        return (int) Math.round(first + fraction * (second - first));
    }

    static void marker(Graphics2D graphics, PointShape shape, double x, double y,
            int size, Color color) {
        if (shape == PointShape.NONE) return;
        int radius = Math.max(2, size / 2);
        graphics.setColor(color);
        switch (shape) {
            case CIRCLE -> graphics.fillOval((int) x - radius, (int) y - radius,
                2 * radius, 2 * radius);
            case SQUARE -> graphics.fillRect((int) x - radius, (int) y - radius,
                2 * radius, 2 * radius);
            case DIAMOND -> {
                Path2D path = new Path2D.Double();
                path.moveTo(x, y - radius);
                path.lineTo(x + radius, y);
                path.lineTo(x, y + radius);
                path.lineTo(x - radius, y);
                path.closePath();
                graphics.fill(path);
            }
            case TRIANGLE_UP, TRIANGLE_DOWN -> {
                double sign = shape == PointShape.TRIANGLE_UP ? -1.0 : 1.0;
                Path2D path = new Path2D.Double();
                path.moveTo(x, y + sign * radius);
                path.lineTo(x + radius, y - sign * radius);
                path.lineTo(x - radius, y - sign * radius);
                path.closePath();
                graphics.fill(path);
            }
            case CROSS -> {
                graphics.setStroke(PlotSupport.strokeWidth(Math.max(1f, size / 5f)));
                graphics.drawLine((int) x - radius, (int) y - radius,
                    (int) x + radius, (int) y + radius);
                graphics.drawLine((int) x - radius, (int) y + radius,
                    (int) x + radius, (int) y - radius);
            }
            case PLUS -> {
                graphics.setStroke(PlotSupport.strokeWidth(Math.max(1f, size / 5f)));
                graphics.drawLine((int) x - radius, (int) y,
                    (int) x + radius, (int) y);
                graphics.drawLine((int) x, (int) y - radius,
                    (int) x, (int) y + radius);
            }
            case NONE -> { }
        }
    }

    static void legend(Graphics2D graphics, PlotOptions options,
            List<String> labels, int x, int y) {
        if (labels.size() < 2 || options.legendPosition() == LegendPosition.NONE) return;
        graphics.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 12));
        int offset = 0;
        for (int index = 0; index < labels.size(); index++) {
            Color color = color(options, index);
            graphics.setColor(color);
            graphics.fillRect(x, y + offset - 9, 14, 10);
            graphics.setColor(TEXT);
            graphics.drawString(labels.get(index), x + 20, y + offset);
            offset += 18;
        }
    }

    static String shorten(String text, FontMetrics metrics, int width) {
        if (metrics.stringWidth(text) <= width) return text;
        String suffix = "…";
        int limit = text.length();
        while (limit > 0 && metrics.stringWidth(text.substring(0, limit) + suffix) > width)
            limit--;
        return text.substring(0, limit) + suffix;
    }
}
