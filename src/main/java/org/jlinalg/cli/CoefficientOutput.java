/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/** Phenotype-only fixed-effect output with BH over non-intercept terms. */
final class CoefficientOutput {
    private CoefficientOutput() { }

    static long write(
            Path path, boolean overwrite, List<String> names,
            double[] beta, double[] standardErrors, double[] statistics,
            double[] degreesOfFreedom, double[] pValues,
            String statisticType, String dfMethod,
            double[] transformedEffects, String transformedEffectName)
            throws IOException {
        try (ExternalBh output = new ExternalBh(path, overwrite)) {
            output.writeHeader(List.of("status", "term", "beta",
                "standard_error", "statistic", "statistic_type",
                "df_numerator", "df_denominator", "df_method",
                "partial_r2", "partial_r2_method", transformedEffectName,
                "p_value"));
            for (int index = 0; index < names.size(); index++) {
                double df = degreesOfFreedom[index];
                double squared = statistics[index] * statistics[index];
                double partial = statisticType.equals("t")
                        && Double.isFinite(df) && df > 0
                    ? squared / (squared + df) : Double.NaN;
                List<String> fields = new ArrayList<>();
                fields.add("ok");
                fields.add(names.get(index));
                fields.add(number(beta[index]));
                fields.add(number(standardErrors[index]));
                fields.add(number(statistics[index]));
                fields.add(statisticType);
                fields.add("1");
                fields.add(number(df));
                fields.add(dfMethod);
                fields.add(number(partial));
                fields.add(Double.isFinite(partial)
                    ? "test-statistic" : "");
                fields.add(transformedEffects == null ? ""
                    : number(transformedEffects[index]));
                fields.add(number(pValues[index]));
                double adjusted = names.get(index).equals("(Intercept)")
                    ? Double.NaN : pValues[index];
                output.write(fields, adjusted);
            }
            long tests = output.tests();
            output.finish();
            return tests;
        }
    }

    private static String number(double value) {
        return Double.isFinite(value) ? Double.toString(value) : "";
    }
}
