/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/** Exhaustive small-order AICc selection for nonseasonal ARIMA. */
public final class AutomaticArima {
    private AutomaticArima() { }

    public static ArimaSelectionResult select(
            double[] series, int maximumAr, int maximumDifferences, int maximumMa) {
        if (maximumAr < 0 || maximumDifferences < 0 || maximumMa < 0
                || maximumAr + maximumMa > 10) {
            throw new IllegalArgumentException("automatic-order bounds are invalid or excessive");
        }
        List<Fit> fits = new ArrayList<>();
        for (int differences = 0; differences <= maximumDifferences; differences++) {
            for (int ar = 0; ar <= maximumAr; ar++) {
                for (int ma = 0; ma <= maximumMa; ma++) {
                    ArimaOrder order = ArimaOrder.arima(ar, differences, ma);
                    ArimaOptions options = ArimaOptions.builder()
                        .includeMean(differences == 0)
                        .includeDrift(differences == 1)
                        .optimizationStarts(5)
                        .build();
                    try {
                        ArimaResult result = Arima.fit(series, order, options);
                        fits.add(new Fit(result, new ArimaCandidate(
                            order, result.aicc(), result.bic(), result.converged())));
                    } catch (IllegalArgumentException ignored) {
                        // Orders invalid for the available series length are omitted.
                    }
                }
            }
        }
        Fit best = fits.stream().filter(value -> value.result().converged())
            .min(Comparator.comparingDouble(value -> value.result().aicc()))
            .orElseGet(() -> fits.stream()
                .min(Comparator.comparingDouble(value -> value.result().aicc()))
                .orElseThrow(() -> new IllegalArgumentException("no ARIMA candidate could be fitted")));
        return new ArimaSelectionResult(best.result(),
            fits.stream().map(Fit::candidate).toList());
    }

    private record Fit(ArimaResult result, ArimaCandidate candidate) { }
}
