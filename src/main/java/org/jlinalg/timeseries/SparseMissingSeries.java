/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;

/** Exact stationary ARMA likelihood on sparse missing observation patterns. */
public final class SparseMissingSeries {
    private SparseMissingSeries() { }

    public static Result fit(double[] series, ArimaOrder order, boolean includeMean,
                             BackendPolicy backendPolicy) {
        if (series == null || order == null || backendPolicy == null) throw new IllegalArgumentException("missing-series inputs are invalid");
        int missing = 0; for (double value : series) if (!Double.isFinite(value)) missing++;
        ExactArmaResult fit = ExactArma.fitPanel(List.of(series), order, includeMean, backendPolicy);
        return new Result(fit, missing, missing == 0 ? "dense Toeplitz" : "sparse observed-pattern Cholesky");
    }

    public record Result(ExactArmaResult fit, int missingCount, String likelihoodPath) { }
}
