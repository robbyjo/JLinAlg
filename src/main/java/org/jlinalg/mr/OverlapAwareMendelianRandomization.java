/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.util.List;

/** Iterative IVW using per-instrument exposure/outcome sampling covariance. */
public final class OverlapAwareMendelianRandomization {
    private OverlapAwareMendelianRandomization() { }

    public static OverlapAwareMrResult ivw(
            List<HarmonizedInstrument> instruments, double[] samplingCovariance) {
        List<HarmonizedInstrument> values = MendelianRandomization.validated(instruments, 2);
        if (samplingCovariance == null || samplingCovariance.length != values.size()) {
            throw new IllegalArgumentException("one sampling covariance is required per instrument");
        }
        double beta = MendelianRandomization.ivw(values, false, 0.95).estimate();
        boolean converged = false;
        int iteration = 0;
        double information = 0.0;
        for (iteration = 1; iteration <= 100; iteration++) {
            double numerator = 0.0;
            information = 0.0;
            for (int index = 0; index < values.size(); index++) {
                HarmonizedInstrument value = values.get(index);
                double covariance = samplingCovariance[index];
                double bound = value.exposureStandardError() * value.outcomeStandardError();
                if (!Double.isFinite(covariance) || Math.abs(covariance) > bound + 1e-12) {
                    throw new IllegalArgumentException("sampling covariance violates Cauchy-Schwarz");
                }
                double variance = value.outcomeStandardError() * value.outcomeStandardError()
                    + beta * beta * value.exposureStandardError()
                        * value.exposureStandardError()
                    - 2.0 * beta * covariance;
                if (!(variance > 0.0)) throw new IllegalArgumentException("nonpositive ratio variance");
                double weight = 1.0 / variance;
                numerator += weight * value.exposureEffect() * value.outcomeEffect();
                information += weight * value.exposureEffect() * value.exposureEffect();
            }
            double candidate = numerator / information;
            if (Math.abs(candidate - beta) <= 1e-10 * (1.0 + Math.abs(beta))) {
                beta = candidate;
                converged = true;
                break;
            }
            beta = candidate;
        }
        MrEstimate estimate = MendelianRandomization.estimate(
            MrMethod.OVERLAP_AWARE_IVW, beta, Math.sqrt(1.0 / information),
            0.95, Double.NaN, 0, 1.0, values.size());
        return new OverlapAwareMrResult(estimate, iteration, converged);
    }
}
