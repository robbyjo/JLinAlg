/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.util.ArrayList;
import java.util.List;

/** Steiger filtering using effect/SE/sample-size derived marginal R-squared. */
public final class SteigerFiltering {
    private SteigerFiltering() { }

    public static SteigerResult analyze(
            List<HarmonizedInstrument> instruments,
            double[] exposureSampleSizes,
            double[] outcomeSampleSizes) {
        List<HarmonizedInstrument> values = MendelianRandomization.validated(instruments, 1);
        if (exposureSampleSizes == null || outcomeSampleSizes == null
                || exposureSampleSizes.length != values.size()
                || outcomeSampleSizes.length != values.size()) {
            throw new IllegalArgumentException("one exposure/outcome sample size is required per instrument");
        }
        List<String> retained = new ArrayList<>();
        List<String> reversed = new ArrayList<>();
        double exposureTotal = 0.0;
        double outcomeTotal = 0.0;
        for (int index = 0; index < values.size(); index++) {
            HarmonizedInstrument value = values.get(index);
            double exposure = rSquared(value.exposureEffect(),
                value.exposureStandardError(), exposureSampleSizes[index]);
            double outcome = rSquared(value.outcomeEffect(),
                value.outcomeStandardError(), outcomeSampleSizes[index]);
            exposureTotal += exposure;
            outcomeTotal += outcome;
            (exposure > outcome ? retained : reversed).add(value.variantId());
        }
        return new SteigerResult(retained, reversed, exposureTotal, outcomeTotal,
            exposureTotal > outcomeTotal);
    }

    private static double rSquared(double beta, double standardError, double sampleSize) {
        if (!(sampleSize > 2.0) || !Double.isFinite(sampleSize)) {
            throw new IllegalArgumentException("sample sizes must be finite and exceed two");
        }
        double f = beta * beta / (standardError * standardError);
        return f / (f + sampleSize - 2.0);
    }
}
