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
        if (estimate != null && (estimate.method() == MrMethod.MR_EGGER
                || estimate.method() == MrMethod.MR_EGGER_GENERALIZED))
            throw new IllegalArgumentException("Egger plotting requires the full MrEggerResult including its intercept");
        writeSvg(destination, instruments, estimate, 0.0, false);
    }

    /** Plots the fitted intercept and exposure-increasing orientation used by Egger. */
    public static void writeSvg(Path destination,
            List<HarmonizedInstrument> instruments, MrEggerResult egger) throws IOException {
        if (egger == null) throw new IllegalArgumentException("Egger result is required");
        writeSvg(destination, instruments, egger.slope(), egger.intercept(), true);
    }

    private static void writeSvg(Path destination, List<HarmonizedInstrument> instruments,
            MrEstimate estimate, double intercept, boolean orient) throws IOException {
        if (destination == null || instruments == null || instruments.isEmpty()
                || estimate == null) throw new IllegalArgumentException("plot inputs are required");
        double xmin = Double.POSITIVE_INFINITY, xmax = Double.NEGATIVE_INFINITY;
        double ymin = Double.POSITIVE_INFINITY, ymax = Double.NEGATIVE_INFINITY;
        for (HarmonizedInstrument value : instruments) {
            double sign = orient ? Math.copySign(1.0, value.exposureEffect()) : 1.0;
            double x = sign * value.exposureEffect(), y = sign * value.outcomeEffect();
            if (!Double.isFinite(x) || !Double.isFinite(y)) throw new IllegalArgumentException("plot effects must be finite");
            xmin = Math.min(xmin, x); xmax = Math.max(xmax, x);
            ymin = Math.min(ymin, y); ymax = Math.max(ymax, y);
        }
        double xspan = Math.max(1e-12, xmax - xmin);
        xmin -= 0.08 * xspan; xmax += 0.08 * xspan;
        ymin = Math.min(ymin, Math.min(intercept + estimate.estimate() * xmin, intercept + estimate.estimate() * xmax));
        ymax = Math.max(ymax, Math.max(intercept + estimate.estimate() * xmin, intercept + estimate.estimate() * xmax));
        double yspan = Math.max(1e-12, ymax - ymin);
        ymin -= 0.08 * yspan; ymax += 0.08 * yspan;
        StringBuilder svg = new StringBuilder("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"720\" height=\"520\" viewBox=\"0 0 720 520\">\n");
        svg.append("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n<g transform=\"translate(70,25)\">");
        svg.append("<rect x=\"0\" y=\"0\" width=\"610\" height=\"430\" fill=\"#fafafa\" stroke=\"#555\"/>\n");
        for (HarmonizedInstrument value : instruments) {
            double sign = orient ? Math.copySign(1.0, value.exposureEffect()) : 1.0;
            double x = 610.0 * (sign * value.exposureEffect() - xmin) / (xmax - xmin);
            double y = 430.0 * (1.0 - (sign * value.outcomeEffect() - ymin) / (ymax - ymin));
            svg.append("<circle cx=\"").append(x).append("\" cy=\"").append(y)
                .append("\" r=\"4\" fill=\"#2563eb\"><title>")
                .append(escape(value.variantId())).append("</title></circle>\n");
        }
        double y0 = 430.0 * (1.0 - (intercept + estimate.estimate() * xmin - ymin) / (ymax - ymin));
        double y1 = 430.0 * (1.0 - (intercept + estimate.estimate() * xmax - ymin) / (ymax - ymin));
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
