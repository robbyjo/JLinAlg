/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import jdistlib.Normal;

/** Dependency-free SVG forest plot for a multivariate MR effect matrix. */
public final class MultivariateMrPlot {
    private MultivariateMrPlot() { }

    public static void writeForestSvg(Path path, MultivariateMrResult result,
            double confidenceLevel) throws IOException {
        if (path == null || result == null || !(confidenceLevel > 0
                && confidenceLevel < 1))
            throw new IllegalArgumentException(
                "path, result, and confidence level in (0,1) are required");
        double critical = Normal.quantile(.5 + confidenceLevel / 2,
            0, 1, true, false);
        double[] beta = result.beta(), se = result.standardErrors();
        double minimum = 0, maximum = 0;
        for (int index = 0; index < beta.length; index++) {
            minimum = Math.min(minimum, beta[index] - critical * se[index]);
            maximum = Math.max(maximum, beta[index] + critical * se[index]);
        }
        if (!(maximum > minimum)) { minimum -= 1; maximum += 1; }
        double padding = .08 * (maximum - minimum);
        minimum -= padding; maximum += padding;
        int rows = beta.length, width = 900, left = 285, right = 45;
        int top = 70, rowHeight = 34, bottom = 65;
        int height = top + rows * rowHeight + bottom;
        double plotWidth = width - left - right;
        StringBuilder svg = new StringBuilder();
        svg.append("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"")
            .append(width).append("\" height=\"").append(height)
            .append("\" viewBox=\"0 0 ").append(width).append(' ')
            .append(height).append("\" role=\"img\" aria-labelledby=\"title desc\">")
            .append("<title id=\"title\">Multivariate Mendelian randomization forest plot</title>")
            .append("<desc id=\"desc\">Exposure by outcome effects with confidence intervals</desc>")
            .append("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>")
            .append("<style>text{font-family:system-ui,sans-serif;fill:#172033}.label{font-size:13px}.axis{font-size:12px;fill:#566176}.ci{stroke:#176b87;stroke-width:2}.point{fill:#d34f2f}</style>");
        double zero = coordinate(0, minimum, maximum, left, plotWidth);
        svg.append("<line x1=\"").append(zero).append("\" x2=\"")
            .append(zero).append("\" y1=\"").append(top - 18)
            .append("\" y2=\"").append(top + rows * rowHeight)
            .append("\" stroke=\"#9aa4b2\" stroke-dasharray=\"4 4\"/>");
        int index = 0;
        for (int outcome = 0; outcome < result.outcomeNames().size(); outcome++)
            for (int exposure = 0; exposure < result.exposureNames().size(); exposure++) {
                int y = top + index * rowHeight;
                double low = beta[index] - critical * se[index];
                double high = beta[index] + critical * se[index];
                svg.append("<text class=\"label\" x=\"15\" y=\"")
                    .append(y + 4).append("\">")
                    .append(escape(result.exposureNames().get(exposure)))
                    .append(" &#x2192; ")
                    .append(escape(result.outcomeNames().get(outcome)))
                    .append("</text><line class=\"ci\" x1=\"")
                    .append(coordinate(low, minimum, maximum, left, plotWidth))
                    .append("\" x2=\"")
                    .append(coordinate(high, minimum, maximum, left, plotWidth))
                    .append("\" y1=\"").append(y).append("\" y2=\"")
                    .append(y).append("\"/><circle class=\"point\" cx=\"")
                    .append(coordinate(beta[index], minimum, maximum, left, plotWidth))
                    .append("\" cy=\"").append(y).append("\" r=\"4\"/>");
                index++;
            }
        svg.append("<line x1=\"").append(left).append("\" x2=\"")
            .append(width - right).append("\" y1=\"")
            .append(top + rows * rowHeight + 5).append("\" y2=\"")
            .append(top + rows * rowHeight + 5).append("\" stroke=\"#172033\"/>")
            .append("<text class=\"axis\" x=\"").append(left)
            .append("\" y=\"").append(height - 25).append("\">")
            .append(String.format(java.util.Locale.ROOT, "%.4g", minimum))
            .append("</text><text class=\"axis\" text-anchor=\"end\" x=\"")
            .append(width - right).append("\" y=\"").append(height - 25)
            .append("\">").append(String.format(java.util.Locale.ROOT,
                "%.4g", maximum)).append("</text></svg>");
        Files.writeString(path, svg.toString(),
            java.nio.file.StandardOpenOption.CREATE_NEW);
    }

    private static double coordinate(double value, double minimum,
            double maximum, int left, double width) {
        return left + (value - minimum) / (maximum - minimum) * width;
    }

    private static String escape(String value) {
        return value.replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace("\"", "&quot;");
    }
}
