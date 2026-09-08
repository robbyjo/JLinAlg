/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import org.jlinalg.compute.BackendPolicy;

/** Joint Gaussian RAM maximum likelihood for observed and latent variables. */
public final class Sem {
    private Sem() { }
    public static SemFitResult fit(double[][] data, SemModel model) {
        return fit(data, model, SemOptions.defaults(), BackendPolicy.PREFERRED);
    }
    public static SemFitResult fit(double[][] data, SemModel model, SemOptions options,
            BackendPolicy backendPolicy) {
        return RamFit.fit(data, model, options, backendPolicy, false);
    }
    /** Covariance is row-major and uses the ML divisor N. */
    public static SemFitResult fitCovariance(double[] covariance, int observations,
            SemModel model, SemOptions options, BackendPolicy backendPolicy) {
        if (model == null || model.hasMeanStructure())
            throw new IllegalArgumentException("use fitMoments for a model with intercepts");
        return RamFit.fitMoments(covariance, new double[model.variables().size()],
            observations, model, options, backendPolicy);
    }
    public static SemFitResult fitMoments(double[] covariance, double[] means, int observations,
            SemModel model, SemOptions options, BackendPolicy backendPolicy) {
        return RamFit.fitMoments(covariance, means, observations, model, options, backendPolicy);
    }
}
