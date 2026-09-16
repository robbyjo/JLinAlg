/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import jdistlib.ChiSquare;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.internal.MatrixOps;

/**
 * Summary-data multivariate MR-IVW for one or more exposures and correlated
 * outcomes. This is deliberately distinct from multivariable MR, which has
 * multiple exposures but only one outcome.
 */
public final class MultivariateMendelianRandomization {
    private MultivariateMendelianRandomization() { }

    public static MultivariateMrResult ivw(
            List<MultivariateInstrument> instruments,
            List<String> exposureNames, List<String> outcomeNames,
            double[][] outcomeCorrelation) {
        return fit(instruments, exposureNames, outcomeNames,
            outcomeCorrelation, false, BackendPolicy.PREFERRED);
    }

    public static MultivariateMrResult fit(
            List<MultivariateInstrument> instruments,
            List<String> exposureNames, List<String> outcomeNames,
            double[][] outcomeCorrelation, boolean multiplicativeRandomEffects,
            BackendPolicy backendPolicy) {
        Dimensions dimensions = validate(instruments, exposureNames,
            outcomeNames, outcomeCorrelation, backendPolicy);
        int variants = dimensions.variants();
        int exposures = dimensions.exposures();
        int outcomes = dimensions.outcomes();
        int parameters = exposures * outcomes;
        double[] correlation = MatrixOps.rowMajor(outcomeCorrelation, outcomes);
        double[] information = new double[parameters * parameters];
        double[] rhs = new double[parameters];

        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            requirePositiveDefinite(correlation, outcomes, backend,
                "outcome correlation");
            for (MultivariateInstrument instrument : instruments) {
                double[] inverse = inverseOutcomeCovariance(
                    instrument, correlation, outcomes, backend);
                double[] x = instrument.exposureEffects();
                double[] y = instrument.outcomeEffects();
                for (int a = 0; a < outcomes; a++) {
                    double weightedOutcome = 0;
                    for (int b = 0; b < outcomes; b++)
                        weightedOutcome += inverse[a * outcomes + b] * y[b];
                    for (int j = 0; j < exposures; j++) {
                        int left = a * exposures + j;
                        rhs[left] += x[j] * weightedOutcome;
                        for (int b = 0; b < outcomes; b++) {
                            double weight = inverse[a * outcomes + b];
                            for (int k = 0; k < exposures; k++) {
                                int right = b * exposures + k;
                                information[left * parameters + right] +=
                                    x[j] * weight * x[k];
                            }
                        }
                    }
                }
            }
            CholeskyFactor factor;
            try {
                factor = backend.dpotrf(information, parameters);
            } catch (IllegalArgumentException | IllegalStateException failure) {
                throw new IllegalArgumentException(
                    "multivariate MR design is rank deficient", failure);
            }
            double[] beta = factor.solve(rhs);
            double[] covariance = factor.solve(
                MatrixOps.identity(parameters), parameters);
            double q = heterogeneity(instruments, beta, correlation,
                exposures, outcomes, backend);
            int qDf = outcomes * (variants - exposures);
            double dispersion = multiplicativeRandomEffects
                ? Math.max(1.0, q / qDf) : 1.0;
            if (dispersion != 1.0)
                for (int index = 0; index < covariance.length; index++)
                    covariance[index] *= dispersion;
            double[] standardErrors = new double[parameters];
            for (int index = 0; index < parameters; index++)
                standardErrors[index] = Math.sqrt(Math.max(0,
                    covariance[index * parameters + index]));
            AssociationStatistics statistics =
                AssociationStatistics.normal(beta, standardErrors);
            MultivariateMrJointTest overall = wald("all_effects", beta,
                covariance, range(parameters), backend);
            List<MultivariateMrJointTest> exposureTests = new ArrayList<>();
            for (int exposure = 0; exposure < exposures; exposure++) {
                int[] indices = new int[outcomes];
                for (int outcome = 0; outcome < outcomes; outcome++)
                    indices[outcome] = outcome * exposures + exposure;
                exposureTests.add(wald(exposureNames.get(exposure), beta,
                    covariance, indices, backend));
            }
            List<MultivariateMrJointTest> outcomeTests = new ArrayList<>();
            for (int outcome = 0; outcome < outcomes; outcome++) {
                int[] indices = new int[exposures];
                for (int exposure = 0; exposure < exposures; exposure++)
                    indices[exposure] = outcome * exposures + exposure;
                outcomeTests.add(wald(outcomeNames.get(outcome), beta,
                    covariance, indices, backend));
            }
            return new MultivariateMrResult(exposureNames, outcomeNames,
                statistics, covariance, marginalStrength(instruments, exposures),
                overall, exposureTests, outcomeTests, q, qDf,
                ChiSquare.cumulative(q, qDf, false, false), dispersion,
                multiplicativeRandomEffects);
        }
    }

    private static double heterogeneity(List<MultivariateInstrument> instruments,
            double[] beta, double[] correlation, int exposures, int outcomes,
            ComputeBackend backend) {
        double q = 0;
        for (MultivariateInstrument instrument : instruments) {
            double[] residual = instrument.outcomeEffects();
            double[] x = instrument.exposureEffects();
            for (int outcome = 0; outcome < outcomes; outcome++)
                for (int exposure = 0; exposure < exposures; exposure++)
                    residual[outcome] -= x[exposure]
                        * beta[outcome * exposures + exposure];
            double[] inverse = inverseOutcomeCovariance(
                instrument, correlation, outcomes, backend);
            double[] solved = MatrixOps.multiply(
                backend, inverse, outcomes, outcomes, residual);
            q += backend.ddot(outcomes, residual, 0, 1, solved, 0, 1);
        }
        return Math.max(0, q);
    }

    static double[] inverseOutcomeCovariance(MultivariateInstrument instrument,
            double[] correlation, int outcomes, ComputeBackend backend) {
        double[] errors = instrument.outcomeStandardErrors();
        double[] covariance = new double[outcomes * outcomes];
        for (int row = 0; row < outcomes; row++)
            for (int column = 0; column < outcomes; column++)
                covariance[row * outcomes + column] = errors[row]
                    * correlation[row * outcomes + column] * errors[column];
        try {
            return backend.dpotrf(covariance, outcomes).solve(
                MatrixOps.identity(outcomes), outcomes);
        } catch (IllegalArgumentException | IllegalStateException failure) {
            throw new IllegalArgumentException(
                "per-instrument outcome covariance is not positive definite",
                failure);
        }
    }

    private static MultivariateMrJointTest wald(String name, double[] beta,
            double[] covariance, int[] indices, ComputeBackend backend) {
        int dimension = beta.length, size = indices.length;
        double[] selected = new double[size];
        double[] selectedCovariance = new double[size * size];
        for (int row = 0; row < size; row++) {
            selected[row] = beta[indices[row]];
            for (int column = 0; column < size; column++)
                selectedCovariance[row * size + column] = covariance[
                    indices[row] * dimension + indices[column]];
        }
        double[] solved;
        try {
            solved = backend.dpotrf(selectedCovariance, size).solve(selected);
        } catch (IllegalArgumentException | IllegalStateException failure) {
            throw new IllegalArgumentException(
                "joint MR coefficient covariance is singular", failure);
        }
        double statistic = Math.max(0, backend.ddot(
            size, selected, 0, 1, solved, 0, 1));
        return new MultivariateMrJointTest(name, statistic, size,
            ChiSquare.cumulative(statistic, size, false, false));
    }

    private static int[] range(int size) {
        int[] result = new int[size];
        for (int index = 0; index < size; index++) result[index] = index;
        return result;
    }

    private static Dimensions validate(List<MultivariateInstrument> instruments,
            List<String> exposureNames, List<String> outcomeNames,
            double[][] outcomeCorrelation, BackendPolicy backendPolicy) {
        if (instruments == null || exposureNames == null || outcomeNames == null
                || exposureNames.isEmpty() || outcomeNames.size() < 2
                || backendPolicy == null)
            throw new IllegalArgumentException(
                "instruments, exposures, at least two outcomes, and a backend are required");
        int variants = instruments.size();
        int exposures = exposureNames.size(), outcomes = outcomeNames.size();
        if (variants <= exposures)
            throw new IllegalArgumentException(
                "more instruments than exposures are required");
        requireNames(exposureNames, "exposure");
        requireNames(outcomeNames, "outcome");
        if (outcomeCorrelation == null || outcomeCorrelation.length != outcomes)
            throw new IllegalArgumentException(
                "outcome correlation dimensions must match outcome names");
        HashSet<String> variantsSeen = new HashSet<>();
        for (int row = 0; row < outcomes; row++)
            if (outcomeCorrelation[row] == null
                    || outcomeCorrelation[row].length != outcomes)
                throw new IllegalArgumentException(
                    "outcome correlation must be square");
        for (int row = 0; row < outcomes; row++) {
            for (int column = 0; column < outcomes; column++) {
                double value = outcomeCorrelation[row][column];
                if (!Double.isFinite(value)
                        || Math.abs(value - outcomeCorrelation[column][row])
                            > 1e-12 * Math.max(1, Math.abs(value)))
                    throw new IllegalArgumentException(
                        "outcome correlation must be finite and symmetric");
            }
            if (Math.abs(outcomeCorrelation[row][row] - 1) > 1e-10)
                throw new IllegalArgumentException(
                    "outcome correlation diagonal must equal one");
        }
        for (MultivariateInstrument instrument : instruments) {
            if (instrument == null || !variantsSeen.add(instrument.variantId())
                    || instrument.exposureEffects().length != exposures
                    || instrument.exposureStandardErrors().length != exposures
                    || instrument.outcomeEffects().length != outcomes
                    || instrument.outcomeStandardErrors().length != outcomes)
                throw new IllegalArgumentException(
                    "multivariate instrument dimensions or variant IDs are invalid");
            requireAssociations(instrument.exposureEffects(),
                instrument.exposureStandardErrors(), "exposure");
            requireAssociations(instrument.outcomeEffects(),
                instrument.outcomeStandardErrors(), "outcome");
        }
        return new Dimensions(variants, exposures, outcomes);
    }

    private static void requireNames(List<String> names, String kind) {
        if (names.stream().anyMatch(name -> name == null || name.isBlank())
                || new HashSet<>(names).size() != names.size())
            throw new IllegalArgumentException(
                kind + " names must be unique and nonblank");
    }

    private static void requireAssociations(double[] effects, double[] errors,
            String kind) {
        for (int index = 0; index < effects.length; index++)
            if (!Double.isFinite(effects[index]) || !(errors[index] > 0)
                    || !Double.isFinite(errors[index]))
                throw new IllegalArgumentException(
                    kind + " associations must be finite with positive standard errors");
    }

    private static void requirePositiveDefinite(double[] matrix, int dimension,
            ComputeBackend backend, String name) {
        try {
            backend.dpotrf(matrix, dimension);
        } catch (IllegalArgumentException | IllegalStateException failure) {
            throw new IllegalArgumentException(name + " must be positive definite",
                failure);
        }
    }

    private static double[] marginalStrength(
            List<MultivariateInstrument> instruments, int exposures) {
        double[] result = new double[exposures];
        for (MultivariateInstrument instrument : instruments) {
            double[] effects = instrument.exposureEffects();
            double[] errors = instrument.exposureStandardErrors();
            for (int exposure = 0; exposure < exposures; exposure++) {
                double z = effects[exposure] / errors[exposure];
                result[exposure] += z * z / instruments.size();
            }
        }
        return result;
    }

    private record Dimensions(int variants, int exposures, int outcomes) { }
}
