/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.util.ArrayList;
import java.util.List;
import jdistlib.ChiSquare;

/** Fast analytic PRESSO-style residual outlier diagnostics. */
public final class MrPresso {
    private MrPresso() { }

    public static MrPressoResult analyze(
            List<HarmonizedInstrument> instruments, double familyWiseAlpha) {
        List<HarmonizedInstrument> values = MendelianRandomization.validated(instruments, 4);
        if (!(familyWiseAlpha > 0.0 && familyWiseAlpha < 1.0)) {
            throw new IllegalArgumentException("alpha must lie in (0,1)");
        }
        MrEstimate raw = MendelianRandomization.ivw(values, false, 0.95);
        List<String> outliers = new ArrayList<>();
        List<HarmonizedInstrument> retained = new ArrayList<>();
        double q = 0.0;
        double[] ratios = new double[values.size()];
        for (int index = 0; index < values.size(); index++) {
            ratios[index] = values.get(index).outcomeEffect()
                / values.get(index).exposureEffect();
        }
        java.util.Arrays.sort(ratios);
        double robustCenter = ratios.length % 2 == 0
            ? 0.5 * (ratios[ratios.length / 2 - 1] + ratios[ratios.length / 2])
            : ratios[ratios.length / 2];
        for (int omitted = 0; omitted < values.size(); omitted++) {
            HarmonizedInstrument value = values.get(omitted);
            double residual = value.outcomeEffect()
                - robustCenter * value.exposureEffect();
            double studentized = residual / value.outcomeStandardError();
            q += studentized * studentized;
            double adjustedP = Math.min(1.0, values.size()
                * MendelianRandomization.normalPValue(studentized));
            if (adjustedP < familyWiseAlpha) outliers.add(value.variantId());
            else retained.add(value);
        }
        MrEstimate corrected = retained.size() >= 2
            ? MendelianRandomization.ivw(retained, false, 0.95) : raw;
        double distortion = raw.estimate() == 0.0 ? Double.NaN
            : 100.0 * (corrected.estimate() - raw.estimate()) / Math.abs(raw.estimate());
        return new MrPressoResult(raw, corrected, outliers, q,
            ChiSquare.cumulative(q, values.size() - 1.0, false, false), distortion);
    }
}
