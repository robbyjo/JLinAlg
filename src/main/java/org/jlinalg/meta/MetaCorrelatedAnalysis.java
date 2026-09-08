/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.Arrays;
import java.util.List;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;

/** GLS pooling with known sampling covariance and optional independent heterogeneity. */
public final class MetaCorrelatedAnalysis {
    private MetaCorrelatedAnalysis() { }

    public static MetaCorrelatedResult fit(List<MetaStudy> studies,
            double[][] samplingCovariance, MetaAnalysisOptions options,
            BackendPolicy backendPolicy) {
        MetaMath.Data data = MetaMath.data(studies);
        if (options == null || backendPolicy == null) throw new IllegalArgumentException("options and backend required");
        int n = studies.size();
        double[] v = MetaGls.matrix(samplingCovariance, n), x = new double[n];
        Arrays.fill(x, 1);
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            double tau = MetaGls.tau(data.effects(), x, 1, v, options, context.backend());
            MetaGls.Fit fit = MetaGls.fit(data.effects(), x, 1, MetaGls.addDiagonal(v, n, tau), context.backend());
            double variance = fit.bread()[0] * MetaAnalysis.inferenceScale(options, fit.q(), n - 1);
            double margin = MetaAnalysis.critical(options, n - 1) * Math.sqrt(variance);
            double q = MetaGls.fit(data.effects(), x, 1, v, context.backend()).q();
            return new MetaCorrelatedResult(studies.stream().map(MetaStudy::name).toList(),
                MetaAnalysis.statistics(options, fit.beta(), new double[] {Math.sqrt(variance)}, n - 1),
                new double[] {variance}, tau, q, n - 1, MetaAnalysis.normalize(fit.wx()), context.provenance(),
                fit.beta()[0] - margin, fit.beta()[0] + margin);
        }
    }
}
