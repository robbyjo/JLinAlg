/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.util.HashSet;
import java.util.List;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.CholeskyFactor;
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
            double[] marginalF = marginalStrength(instruments, exposures);
            return new MultivariableMrResult(exposureNames,
                AssociationStatistics.normal(beta, se), exposureCovariance,
                marginalF, null,
                intercept ? solution.coefficients()[0] : 0.0,
                intercept ? Math.sqrt(Math.max(0.0, covariance[0])) : Double.NaN,
                q, df);
        }
    }

    /**
     * Generalized multivariable IVW/Egger with a known instrument-by-instrument
     * outcome covariance and per-instrument cross-exposure sampling covariance.
     * The latter is used for genuine conditional-strength diagnostics.
     */
    public static MultivariableMrResult generalizedFit(
            List<MultivariableInstrument> instruments,
            List<String> exposureNames, boolean intercept,
            double[][] outcomeCovariance,
            List<double[][]> exposureSamplingCovariances,
            BackendPolicy backendPolicy) {
        int exposures = validate(instruments, exposureNames);
        int rows = instruments.size();
        int columns = exposures + (intercept ? 1 : 0);
        if (rows <= columns || backendPolicy == null)
            throw new IllegalArgumentException("more instruments than coefficients and a backend are required");
        double[] covariance = org.jlinalg.internal.MatrixOps.rowMajor(
            outcomeCovariance, rows);
        if (outcomeCovariance[0].length != rows)
            throw new IllegalArgumentException("outcome covariance must be square");
        double[] design = new double[rows * columns];
        double[] outcome = new double[rows];
        for (int row = 0; row < rows; row++) {
            if (intercept) design[row * columns] = 1.0;
            double[] effects = instruments.get(row).exposureEffects();
            for (int exposure = 0; exposure < exposures; exposure++)
                design[row * columns + exposure + (intercept ? 1 : 0)] = effects[exposure];
            outcome[row] = instruments.get(row).outcomeEffect();
            for (int column = 0; column < rows; column++) {
                double a = covariance[row * rows + column];
                double b = covariance[column * rows + row];
                if (!Double.isFinite(a) || Math.abs(a - b) > 1e-12 * Math.max(1.0, Math.max(Math.abs(a), Math.abs(b))))
                    throw new IllegalArgumentException("outcome covariance must be finite and symmetric");
            }
        }
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            CholeskyFactor vf;
            try { vf = backend.dpotrf(covariance, rows); }
            catch (IllegalArgumentException | IllegalStateException failure) {
                throw new IllegalArgumentException("outcome covariance must be positive definite", failure);
            }
            double[] inverseDesign = vf.solve(design, columns);
            double[] inverseOutcome = vf.solve(outcome);
            double[] information = new double[columns * columns];
            double[] rhs = new double[columns];
            for (int a = 0; a < columns; a++) for (int row = 0; row < rows; row++) {
                rhs[a] += design[row * columns + a] * inverseOutcome[row];
                for (int b = 0; b < columns; b++)
                    information[a * columns + b] += design[row * columns + a] * inverseDesign[row * columns + b];
            }
            CholeskyFactor factor;
            try { factor = backend.dpotrf(information, columns); }
            catch (IllegalArgumentException | IllegalStateException failure) {
                throw new IllegalArgumentException("generalized multivariable design is rank deficient", failure);
            }
            double[] coefficients = factor.solve(rhs);
            double[] coefficientCovariance = factor.solve(
                org.jlinalg.internal.MatrixOps.identity(columns), columns);
            double[] residual = outcome.clone();
            for (int row = 0; row < rows; row++) for (int column = 0; column < columns; column++)
                residual[row] -= design[row * columns + column] * coefficients[column];
            double[] inverseResidual = vf.solve(residual);
            double q = 0;
            for (int row = 0; row < rows; row++) q += residual[row] * inverseResidual[row];
            int df = rows - columns;
            double dispersion = Math.max(1.0, q / df);
            int offset = intercept ? 1 : 0;
            double[] beta = new double[exposures], se = new double[exposures];
            double[] effectCovariance = new double[exposures * exposures];
            for (int a = 0; a < exposures; a++) {
                beta[a] = coefficients[a + offset];
                se[a] = Math.sqrt(dispersion * coefficientCovariance[(a + offset) * columns + a + offset]);
                for (int b = 0; b < exposures; b++)
                    effectCovariance[a * exposures + b] = dispersion * coefficientCovariance[(a + offset) * columns + b + offset];
            }
            ConditionalStrengthResult strength = conditionalStrength(
                instruments, exposureNames, exposureSamplingCovariances,
                backendPolicy);
            return new MultivariableMrResult(exposureNames,
                AssociationStatistics.normal(beta, se), effectCovariance,
                marginalStrength(instruments, exposures), strength.fStatistics(),
                intercept ? coefficients[0] : 0,
                intercept ? Math.sqrt(dispersion * coefficientCovariance[0]) : Double.NaN,
                q, df);
        }
    }

    /** Conditional instrument strength using full within-instrument exposure covariance. */
    public static ConditionalStrengthResult conditionalStrength(
            List<MultivariableInstrument> instruments,
            List<String> exposureNames,
            List<double[][]> exposureSamplingCovariances,
            BackendPolicy backendPolicy) {
        int exposures = validate(instruments, exposureNames);
        int variants = instruments.size();
        if (exposureSamplingCovariances == null
                || exposureSamplingCovariances.size() != variants
                || backendPolicy == null)
            throw new IllegalArgumentException("one exposure sampling covariance matrix and a backend are required per instrument");
        double[][][] covariance = new double[variants][exposures][exposures];
        double[][] correlations = new double[variants][exposures * exposures];
        for (int variant = 0; variant < variants; variant++) {
            double[][] supplied = exposureSamplingCovariances.get(variant);
            if (supplied == null || supplied.length != exposures)
                throw new IllegalArgumentException("exposure sampling covariance dimensions are invalid");
            // Validate every row before inspecting transposed entries.
            for (int a = 0; a < exposures; a++)
                if (supplied[a] == null || supplied[a].length != exposures)
                    throw new IllegalArgumentException("exposure sampling covariance dimensions are invalid");
            double[] errors = instruments.get(variant).exposureStandardErrors();
            for (int a = 0; a < exposures; a++) {
                for (int b = 0; b < exposures; b++) {
                    double value = supplied[a][b];
                    // Work in correlation units so validation is independent
                    // of exposure units, including very small covariances.
                    double large = Math.max(errors[a], errors[b]);
                    double small = Math.min(errors[a], errors[b]);
                    double standardized = value / large / small;
                    double reverse = supplied[b][a] / large / small;
                    if (!Double.isFinite(standardized) || !Double.isFinite(reverse)
                            || Math.abs(standardized - reverse) > 1e-12 * Math.max(1.0, Math.abs(standardized)))
                        throw new IllegalArgumentException("exposure sampling covariance must be finite and symmetric");
                    covariance[variant][a][b] = .5*value + .5*supplied[b][a];
                    correlations[variant][a * exposures + b] = .5*standardized + .5*reverse;
                }
                if (Math.abs(correlations[variant][a * exposures + a] - 1) > 1e-8)
                    throw new IllegalArgumentException("exposure covariance diagonal must match reported standard errors");
            }
        }
        double[] f = new double[exposures], q = new double[exposures];
        int degrees = variants - exposures + 1;
        if (degrees < 1) throw new IllegalArgumentException("conditional strength has no residual degrees of freedom");
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            for (double[] correlation : correlations)
                for (double eigenvalue : backend.dsyev(correlation, exposures).eigenvalues())
                    if (!Double.isFinite(eigenvalue) || eigenvalue < -1e-12 * exposures)
                        throw new IllegalArgumentException("exposure sampling covariance must be positive semidefinite");
            for (int target = 0; target < exposures; target++) {
                int predictors = exposures - 1;
                double[] gamma = new double[predictors];
                for (int iteration = 0; iteration < 100; iteration++) {
                    double[] information = new double[predictors * predictors];
                    double[] rhs = new double[predictors];
                    for (int variant = 0; variant < variants; variant++) {
                        double variance = covariance[variant][target][target];
                        for (int a = 0, ia = 0; a < exposures; a++) if (a != target) {
                            variance -= 2 * gamma[ia] * covariance[variant][target][a];
                            for (int b = 0, ib = 0; b < exposures; b++) if (b != target)
                                variance += gamma[ia] * gamma[ib++] * covariance[variant][a][b];
                            ia++;
                        }
                        if (!(variance > 0) || !Double.isFinite(variance))
                            throw new IllegalArgumentException("conditional exposure residual variance is not positive");
                        double weight = 1 / variance;
                        double y = instruments.get(variant).exposureEffects()[target];
                        for (int a = 0, ia = 0; a < exposures; a++) if (a != target) {
                            double xa = instruments.get(variant).exposureEffects()[a];
                            rhs[ia] += weight * xa * y;
                            for (int b = 0, ib = 0; b < exposures; b++) if (b != target)
                                information[ia * predictors + ib++] += weight * xa * instruments.get(variant).exposureEffects()[b];
                            ia++;
                        }
                    }
                    double[] next = predictors == 0 ? new double[0]
                        : backend.dpotrf(information, predictors).solve(rhs);
                    double change = 0;
                    for (int i = 0; i < predictors; i++) change = Math.max(change, Math.abs(next[i] - gamma[i]));
                    gamma = next;
                    if (change <= 1e-10 * (1 + java.util.Arrays.stream(gamma).map(Math::abs).max().orElse(0))) break;
                    if (iteration == 99) throw new ArithmeticException("conditional-strength weighting did not converge");
                }
                double statistic = 0;
                for (int variant = 0; variant < variants; variant++) {
                    double residual = instruments.get(variant).exposureEffects()[target];
                    double variance = covariance[variant][target][target];
                    for (int a = 0, ia = 0; a < exposures; a++) if (a != target) {
                        residual -= gamma[ia] * instruments.get(variant).exposureEffects()[a];
                        variance -= 2 * gamma[ia] * covariance[variant][target][a];
                        for (int b = 0, ib = 0; b < exposures; b++) if (b != target)
                            variance += gamma[ia] * gamma[ib++] * covariance[variant][a][b];
                        ia++;
                    }
                    if (!(variance > 0) || !Double.isFinite(variance))
                        throw new IllegalArgumentException("conditional exposure residual variance is not positive");
                    statistic += residual * residual / variance;
                }
                q[target] = statistic;
                f[target] = statistic / degrees;
            }
        }
        return new ConditionalStrengthResult(exposureNames, f, q, degrees);
    }

    private static int validate(
            List<MultivariableInstrument> instruments, List<String> names) {
        if (instruments == null || instruments.size() < 3 || names == null || names.isEmpty()) {
            throw new IllegalArgumentException("instruments and exposure names are required");
        }
        int exposures = names.size();
        if (names.stream().anyMatch(n -> n == null || n.isBlank())
                || new HashSet<>(names).size() != exposures)
            throw new IllegalArgumentException("exposure names must be unique and nonblank");
        HashSet<String> variants = new HashSet<>();
        for (MultivariableInstrument value : instruments) {
            if (value == null || value.variantId() == null || !variants.add(value.variantId())
                    || value.exposureEffects() == null
                    || value.exposureEffects().length != exposures
                    || value.exposureStandardErrors().length != exposures
                    || !Double.isFinite(value.outcomeEffect())
                    || !(value.outcomeStandardError() > 0.0)
                    || !Double.isFinite(value.outcomeStandardError())) {
                throw new IllegalArgumentException("invalid multivariable instrument");
            }
            double[] effects = value.exposureEffects(), errors = value.exposureStandardErrors();
            for (int index = 0; index < exposures; index++) {
                if (!Double.isFinite(effects[index])
                        || !(errors[index] > 0.0) || !Double.isFinite(errors[index])) {
                    throw new IllegalArgumentException("invalid exposure association");
                }
            }
        }
        return exposures;
    }

    private static double[] marginalStrength(
            List<MultivariableInstrument> instruments, int exposures) {
        double[] result = new double[exposures];
        for (MultivariableInstrument value : instruments) {
            double[] effects = value.exposureEffects(), errors = value.exposureStandardErrors();
            for (int exposure = 0; exposure < exposures; exposure++) {
                double z = effects[exposure] / errors[exposure];
                result[exposure] += z * z / instruments.size();
            }
        }
        return result;
    }
}
