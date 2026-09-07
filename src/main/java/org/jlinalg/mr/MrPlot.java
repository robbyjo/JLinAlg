/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

/** Dependency-free SVG scatter plot for summary-data MR. */
public final class MrPlot {
    private MrPlot() { }

    public static void writeSvg(Path destination,
            List<HarmonizedInstrument> instruments, MrEstimate estimate)
            throws IOException {
        if (destination == null || instruments == null || instruments.isEmpty()
                || estimate == null) throw new IllegalArgumentException("plot inputs are required");
        double xmin = Double.POSITIVE_INFINITY, xmax = Double.NEGATIVE_INFINITY;
        double ymin = Double.POSITIVE_INFINITY, ymax = Double.NEGATIVE_INFINITY;
        for (HarmonizedInstrument value : instruments) {
            xmin = Math.min(xmin, value.exposureEffect()); xmax = Math.max(xmax, value.exposureEffect());
            ymin = Math.min(ymin, value.outcomeEffect()); ymax = Math.max(ymax, value.outcomeEffect());
        }
        double xspan = Math.max(1e-12, xmax - xmin), yspan = Math.max(1e-12, ymax - ymin);
        xmin -= 0.08 * xspan; xmax += 0.08 * xspan; ymin -= 0.08 * yspan; ymax += 0.08 * yspan;
        StringBuilder svg = new StringBuilder("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"720\" height=\"520\" viewBox=\"0 0 720 520\">\n");
        svg.append("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n<g transform=\"translate(70,25)\">");
        svg.append("<rect x=\"0\" y=\"0\" width=\"610\" height=\"430\" fill=\"#fafafa\" stroke=\"#555\"/>\n");
        for (HarmonizedInstrument value : instruments) {
            double x = 610.0 * (value.exposureEffect() - xmin) / (xmax - xmin);
            double y = 430.0 * (1.0 - (value.outcomeEffect() - ymin) / (ymax - ymin));
            svg.append("<circle cx=\"").append(x).append("\" cy=\"").append(y)
                .append("\" r=\"4\" fill=\"#2563eb\"><title>")
                .append(escape(value.variantId())).append("</title></circle>\n");
        }
        double y0 = 430.0 * (1.0 - (estimate.estimate() * xmin - ymin) / (ymax - ymin));
        double y1 = 430.0 * (1.0 - (estimate.estimate() * xmax - ymin) / (ymax - ymin));
        svg.append("<line x1=\"0\" y1=\"").append(y0).append("\" x2=\"610\" y2=\"").append(y1)
            .append("\" stroke=\"#dc2626\" stroke-width=\"2\"/>\n</g>")
            .append("<text x=\"360\" y=\"505\" text-anchor=\"middle\">Exposure effect</text>")
            .append("<text transform=\"translate(16 240) rotate(-90)\" text-anchor=\"middle\">Outcome effect</text></svg>\n");
        Path parent = destination.toAbsolutePath().normalize().getParent();
        if (parent != null) Files.createDirectories(parent);
        Files.writeString(destination, svg.toString());
    }

    private static String escape(String value) {
        return value.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace("\"", "&quot;");
    }
}
