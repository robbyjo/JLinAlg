/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.util.List;

/** Robust summary-data estimators accounting for exposure measurement error. */
public final class RobustMendelianRandomization {
    private static final double HUBER = 1.345;
    private RobustMendelianRandomization() { }

    /** Fits a Huber robust adjusted-profile-score estimator. */
    public static MrRapsResult raps(List<HarmonizedInstrument> instruments) {
        List<HarmonizedInstrument> values = MendelianRandomization.validated(instruments, 3);
        double beta = MendelianRandomization.ivw(values, false, 0.95).estimate();
        double tau2 = 0.0;
        boolean converged = false;
        int iteration = 0;
        for (iteration = 1; iteration <= 100; iteration++) {
            double score = score(values, beta, tau2);
            double step = Math.max(1e-6, 1e-5 * (1.0 + Math.abs(beta)));
            double derivative = (score(values, beta + step, tau2)
                - score(values, beta - step, tau2)) / (2.0 * step);
            if (!Double.isFinite(derivative) || Math.abs(derivative) < 1e-12) break;
            double candidate = beta - score / derivative;
            double limit = 2.0 * (1.0 + Math.abs(beta));
            candidate = Math.max(beta - limit, Math.min(beta + limit, candidate));
            double newTau = overdispersion(values, candidate);
            if (Math.abs(candidate - beta) <= 1e-9 * (1.0 + Math.abs(beta))
                    && Math.abs(newTau - tau2) <= 1e-9 * (1.0 + tau2)) {
                beta = candidate;
                tau2 = newTau;
                converged = true;
                break;
            }
            beta = candidate;
            tau2 = newTau;
        }
        double h = Math.max(1e-6, 1e-5 * (1.0 + Math.abs(beta)));
        double derivative = (score(values, beta + h, tau2)
            - score(values, beta - h, tau2)) / (2.0 * h);
        double meat = 0.0;
        for (HarmonizedInstrument value : values) {
            double variance = value.outcomeStandardError() * value.outcomeStandardError()
                + beta * beta * value.exposureStandardError() * value.exposureStandardError()
                + tau2;
            double residual = value.outcomeEffect() - beta * value.exposureEffect();
            double psi = huber(residual / Math.sqrt(variance));
            double contribution = value.exposureEffect() * psi / Math.sqrt(variance);
            meat += contribution * contribution;
        }
        double standardError = Math.sqrt(meat / (derivative * derivative));
        MrEstimate estimate = MendelianRandomization.estimate(MrMethod.MR_RAPS,
            beta, standardError, 0.95, Double.NaN, 0, 1.0, values.size());
        return new MrRapsResult(estimate, tau2, iteration, converged);
    }

    private static double score(
            List<HarmonizedInstrument> values, double beta, double tau2) {
        double result = 0.0;
        for (HarmonizedInstrument value : values) {
            double seX2 = value.exposureStandardError() * value.exposureStandardError();
            double variance = value.outcomeStandardError() * value.outcomeStandardError()
                + beta * beta * seX2 + tau2;
            double residual = value.outcomeEffect() - beta * value.exposureEffect();
            double standardized = residual / Math.sqrt(variance);
            result += huber(standardized)
                * (value.exposureEffect() / Math.sqrt(variance)
                    + beta * seX2 * residual / (variance * Math.sqrt(variance)));
        }
        return result;
    }

    private static double overdispersion(
            List<HarmonizedInstrument> values, double beta) {
        double excess = 0.0;
        for (HarmonizedInstrument value : values) {
            double residual = value.outcomeEffect() - beta * value.exposureEffect();
            double known = value.outcomeStandardError() * value.outcomeStandardError()
                + beta * beta * value.exposureStandardError() * value.exposureStandardError();
            excess += Math.max(0.0, residual * residual - known);
        }
        return excess / values.size();
    }

    private static double huber(double value) {
        return Math.max(-HUBER, Math.min(HUBER, value));
    }
}
