/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.plot;

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.geom.Ellipse2D;
import org.knowm.xchart.internal.chartpart.IChart;

/** Deterministic, area-neutral two-set Venn diagram. */
final class VennChart implements IChart {
    private final DiagramPlots.VennData data;
    private final PlotOptions options;

    VennChart(DiagramPlots.VennData data, PlotOptions options) {
        if (data == null || data.leftLabel() == null || data.leftLabel().isBlank()
                || data.rightLabel() == null || data.rightLabel().isBlank())
            throw new IllegalArgumentException("Venn sets need nonblank labels");
        if (data.leftOnly() < 0 || data.overlap() < 0 || data.rightOnly() < 0)
            throw new IllegalArgumentException("Venn region counts must be nonnegative");
        if ((long) data.leftOnly() + data.overlap() + data.rightOnly() == 0)
            throw new IllegalArgumentException("at least one Venn region must be nonempty");
        this.data = data;
        this.options = options;
    }

    @Override public int getWidth() { return options.width(); }
    @Override public int getHeight() { return options.height(); }
    @Override public String getTitle() { return options.title(); }

    @Override
    public void paint(Graphics2D source, int width, int height) {
        Graphics2D graphics = DiagramSupport.begin(source, options, width, height);
        try {
            int top = DiagramSupport.top(options) + 25;
            int usableHeight = Math.max(120, height - top - 75);
            double diameter = Math.min(width * 0.52, usableHeight * 0.92);
            double overlap = diameter * 0.38;
            double totalWidth = 2.0 * diameter - overlap;
            double leftX = (width - totalWidth) / 2.0;
            double rightX = leftX + diameter - overlap;
            double y = top + (usableHeight - diameter) / 2.0;
            Ellipse2D left = new Ellipse2D.Double(leftX, y, diameter, diameter);
            Ellipse2D right = new Ellipse2D.Double(rightX, y, diameter, diameter);

            Color leftColor = DiagramSupport.color(options, 0);
            Color rightColor = DiagramSupport.color(options, 1);
            graphics.setComposite(AlphaComposite.SrcOver.derive(0.42f));
            graphics.setColor(leftColor);
            graphics.fill(left);
            graphics.setColor(rightColor);
            graphics.fill(right);
            graphics.setComposite(AlphaComposite.SrcOver);
            graphics.setStroke(PlotSupport.strokeWidth(2f));
            graphics.setColor(leftColor.darker());
            graphics.draw(left);
            graphics.setColor(rightColor.darker());
            graphics.draw(right);

            graphics.setColor(DiagramSupport.TEXT);
            graphics.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
            drawCentered(graphics, Integer.toString(data.leftOnly()),
                leftX + (diameter - overlap) * 0.48, y + diameter * 0.52);
            drawCentered(graphics, Integer.toString(data.overlap()),
                leftX + diameter - overlap / 2.0, y + diameter * 0.52);
            drawCentered(graphics, Integer.toString(data.rightOnly()),
                rightX + overlap + (diameter - overlap) * 0.52,
                y + diameter * 0.52);

            graphics.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 13));
            drawCentered(graphics, data.leftLabel(), leftX + diameter * 0.38,
                y + diameter + 28);
            drawCentered(graphics, data.rightLabel(), rightX + diameter * 0.62,
                y + diameter + 28);
            graphics.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 11));
            drawCentered(graphics, "Set sizes: " + data.leftSize() + " and "
                + data.rightSize(), width / 2.0, y + diameter + 49);
        } finally {
            graphics.dispose();
        }
    }

    private static void drawCentered(Graphics2D graphics, String text,
            double centerX, double baseline) {
        FontMetrics metrics = graphics.getFontMetrics();
        graphics.drawString(text, (float) (centerX - metrics.stringWidth(text) / 2.0),
            (float) baseline);
    }
}
