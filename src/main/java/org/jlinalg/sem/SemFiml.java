/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.model.MissingDataPolicy;

/** Direct constrained Gaussian observed-data maximum likelihood over missingness patterns. */
public final class SemFiml {
    private SemFiml() { }
    public static SemFimlResult fit(double[][] data, SemModel model) {
        return fit(data, model, 10000, 1e-8, BackendPolicy.PREFERRED);
    }
    /** The legacy iteration limit now bounds objective evaluations, with a minimum budget of 100. */
    public static SemFimlResult fit(double[][] data, SemModel model, int maximumIterations,
            double tolerance, BackendPolicy backendPolicy) {
        if(maximumIterations<1)throw new IllegalArgumentException("positive evaluation limit required");
        SemFitResult fit=RamFit.fit(data,model,new SemOptions(Math.max(100,maximumIterations),
            tolerance,MissingDataPolicy.ERROR),backendPolicy,true);
        return new SemFimlResult(fit,fit.impliedMeans(),fit.impliedCovariance(),
            fit.state().patterns().size(),fit.functionEvaluations(),fit.converged());
    }
}
