/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.util.HashSet;
import java.util.List;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.internal.LeastSquaresSolver;

/** Weighted multivariable IVW and multivariable MR-Egger. */
public final class MultivariableMendelianRandomization {
    private MultivariableMendelianRandomization() { }

    public static MultivariableMrResult ivw(
            List<MultivariableInstrument> instruments, List<String> exposureNames) {
        return fit(instruments, exposureNames, false, BackendPolicy.PREFERRED);
    }

    public static MultivariableMrResult egger(
            List<MultivariableInstrument> instruments, List<String> exposureNames) {
        return fit(instruments, exposureNames, true, BackendPolicy.PREFERRED);
    }

    public static MultivariableMrResult fit(
            List<MultivariableInstrument> instruments,
            List<String> exposureNames,
            boolean intercept,
            BackendPolicy backendPolicy) {
        int exposures = validate(instruments, exposureNames);
        int rows = instruments.size();
        int columns = exposures + (intercept ? 1 : 0);
        if (rows <= columns) {
            throw new IllegalArgumentException("more instruments than fitted coefficients are required");
        }
        double[] design = new double[rows * columns];
        double[] outcome = new double[rows];
        for (int row = 0; row < rows; row++) {
            MultivariableInstrument value = instruments.get(row);
            double scale = 1.0 / value.outcomeStandardError();
            if (intercept) design[row * columns] = scale;
            double[] effects = value.exposureEffects();
            for (int exposure = 0; exposure < exposures; exposure++) {
                design[row * columns + exposure + (intercept ? 1 : 0)] =
                    effects[exposure] * scale;
            }
            outcome[row] = value.outcomeEffect() * scale;
        }
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            LeastSquaresSolver.Solution solution = LeastSquaresSolver.solve(
                design, outcome, rows, columns, false, backend);
            double[] fitted = org.jlinalg.internal.MatrixOps.multiply(
                backend, design, rows, columns, solution.coefficients());
            double q = 0.0;
            for (int row = 0; row < rows; row++) {
                double residual = outcome[row] - fitted[row];
                q += residual * residual;
            }
            int df = rows - columns;
            double dispersion = Math.max(1.0, q / df);
            double[] covariance = solution.unscaledCovariance().clone();
            for (int index = 0; index < covariance.length; index++) covariance[index] *= dispersion;
            int offset = intercept ? 1 : 0;
            double[] beta = new double[exposures];
            double[] se = new double[exposures];
            double[] exposureCovariance = new double[exposures * exposures];
            for (int first = 0; first < exposures; first++) {
                beta[first] = solution.coefficients()[first + offset];
                se[first] = Math.sqrt(Math.max(0.0,
                    covariance[(first + offset) * columns + first + offset]));
                for (int second = 0; second < exposures; second++) {
                    exposureCovariance[first * exposures + second] =
                        covariance[(first + offset) * columns + second + offset];
                }
            }
            double[] conditionalF = conditionalStrength(instruments, exposures);
            return new MultivariableMrResult(exposureNames,
                AssociationStatistics.normal(beta, se), exposureCovariance,
                conditionalF,
                intercept ? solution.coefficients()[0] : 0.0,
                intercept ? Math.sqrt(Math.max(0.0, covariance[0])) : Double.NaN,
                q, df);
        }
    }

    private static int validate(
            List<MultivariableInstrument> instruments, List<String> names) {
        if (instruments == null || instruments.size() < 3 || names == null || names.isEmpty()) {
            throw new IllegalArgumentException("instruments and exposure names are required");
        }
        int exposures = names.size();
        HashSet<String> variants = new HashSet<>();
        for (MultivariableInstrument value : instruments) {
            if (value == null || value.variantId() == null || !variants.add(value.variantId())
                    || value.exposureEffects() == null
                    || value.exposureEffects().length != exposures
                    || value.exposureStandardErrors().length != exposures
                    || !Double.isFinite(value.outcomeEffect())
                    || !(value.outcomeStandardError() > 0.0)) {
                throw new IllegalArgumentException("invalid multivariable instrument");
            }
            for (int index = 0; index < exposures; index++) {
                if (!Double.isFinite(value.exposureEffects()[index])
                        || !(value.exposureStandardErrors()[index] > 0.0)) {
                    throw new IllegalArgumentException("invalid exposure association");
                }
            }
        }
        return exposures;
    }

    private static double[] conditionalStrength(
            List<MultivariableInstrument> instruments, int exposures) {
        double[] result = new double[exposures];
        for (int exposure = 0; exposure < exposures; exposure++) {
            double sum = 0.0;
            for (MultivariableInstrument value : instruments) {
                double z = value.exposureEffects()[exposure]
                    / value.exposureStandardErrors()[exposure];
                sum += z * z;
            }
            result[exposure] = sum / instruments.size();
        }
        return result;
    }
}
