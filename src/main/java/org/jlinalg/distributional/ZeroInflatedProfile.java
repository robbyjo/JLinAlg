/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.distributional;

import java.util.Arrays;

/** Profile marginal-Laplace log likelihood for one optimizer-scale parameter. */
public record ZeroInflatedProfile(
        String parameter,
        double[] values,
        double[] logLikelihoods,
        boolean[] converged) {
    public ZeroInflatedProfile {
        if (parameter == null || parameter.isBlank() || values == null
                || logLikelihoods == null || converged == null
                || values.length == 0 || values.length != logLikelihoods.length
                || values.length != converged.length)
            throw new IllegalArgumentException("profile arrays must match");
        values = values.clone();
        logLikelihoods = logLikelihoods.clone();
        converged = converged.clone();
    }

    @Override public double[] values() { return values.clone(); }
    @Override public double[] logLikelihoods() { return logLikelihoods.clone(); }
    @Override public boolean[] converged() { return converged.clone(); }

    /** Likelihood-ratio support interval over the supplied ordered grid. */
    public double[] supportInterval(double cutoff) {
        if (!(cutoff > 0.0) || !Double.isFinite(cutoff))
            throw new IllegalArgumentException("cutoff must be finite and positive");
        double maximum = Arrays.stream(logLikelihoods).max().orElseThrow();
        double lower = Double.NaN;
        double upper = Double.NaN;
        for (int index = 0; index < values.length; index++)
            if (converged[index] && 2.0 * (maximum - logLikelihoods[index])
                    <= cutoff) {
                if (Double.isNaN(lower)) lower = values[index];
                upper = values[index];
            }
        return new double[] {lower, upper};
    }
}
