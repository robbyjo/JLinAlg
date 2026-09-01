/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

import java.util.List;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTranspose;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;

/** Generalized summary-data MR for instruments correlated through LD. */
public final class CorrelatedMendelianRandomization {
    private CorrelatedMendelianRandomization() {
    }

    /** Fits generalized IVW with the preferred backend. */
    public static CorrelatedMrEstimate ivw(
            List<HarmonizedInstrument> instruments,
            double[][] alleleAlignedCorrelation,
            boolean multiplicativeRandomEffects) {
        return ivw(instruments, alleleAlignedCorrelation,
            multiplicativeRandomEffects, 0.95, BackendPolicy.PREFERRED);
    }

    /** Fits generalized IVW using an allele-aligned LD correlation matrix. */
    public static CorrelatedMrEstimate ivw(
            List<HarmonizedInstrument> instruments,
            double[][] alleleAlignedCorrelation,
            boolean multiplicativeRandomEffects,
            double confidenceLevel,
            BackendPolicy backendPolicy) {
        List<HarmonizedInstrument> values =
            MendelianRandomization.validated(instruments, 2);
        MendelianRandomization.validateConfidence(confidenceLevel);
        double[] correlation = correlation(
            alleleAlignedCorrelation, values.size());
        if (backendPolicy == null) {
            throw new IllegalArgumentException("backendPolicy is required");
        }
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            double[] covariance = outcomeCovariance(values, correlation, false);
            CholeskyFactor factor = positiveDefiniteFactor(
                covariance, values.size(), backend);
            double[] exposure = exposure(values, false);
            double[] outcome = outcome(values, false);
            double[] inverseExposure = factor.solve(exposure);
            double[] inverseOutcome = factor.solve(outcome);
            double denominator = backend.ddot(values.size(), exposure, 0, 1,
                inverseExposure, 0, 1);
            double numerator = backend.ddot(values.size(), exposure, 0, 1,
                inverseOutcome, 0, 1);
            if (!(denominator > 0.0) || !Double.isFinite(denominator)) {
                throw new IllegalArgumentException(
                    "generalized IVW information is not finite and positive");
            }
            double estimate = numerator / denominator;
            double[] residual = new double[values.size()];
            for (int index = 0; index < values.size(); index++) {
                residual[index] = outcome[index] - estimate * exposure[index];
            }
            double[] inverseResidual = factor.solve(residual);
            double q = backend.ddot(values.size(), residual, 0, 1,
                inverseResidual, 0, 1);
            int degreesOfFreedom = values.size() - 1;
            double dispersion = multiplicativeRandomEffects
                ? Math.max(1.0, q / degreesOfFreedom) : 1.0;
            double standardError = Math.sqrt(dispersion / denominator);
            MrMethod method = multiplicativeRandomEffects
                ? MrMethod.IVW_GENERALIZED_MULTIPLICATIVE_RANDOM
                : MrMethod.IVW_GENERALIZED_FIXED;
            MrEstimate result = MendelianRandomization.estimate(
                method, estimate, standardError,
                confidenceLevel, q, degreesOfFreedom, dispersion, values.size());
            return new CorrelatedMrEstimate(result, context.provenance());
        }
    }

    /** Fits generalized MR-Egger with the preferred backend. */
    public static CorrelatedMrEggerResult egger(
            List<HarmonizedInstrument> instruments,
            double[][] alleleAlignedCorrelation) {
        return egger(instruments, alleleAlignedCorrelation,
            0.95, BackendPolicy.PREFERRED);
    }

    /**
     * Fits generalized MR-Egger. Correlation signs are changed internally when
     * variants are reoriented to exposure-increasing alleles.
     */
    public static CorrelatedMrEggerResult egger(
            List<HarmonizedInstrument> instruments,
            double[][] alleleAlignedCorrelation,
            double confidenceLevel,
            BackendPolicy backendPolicy) {
        List<HarmonizedInstrument> values =
            MendelianRandomization.validated(instruments, 3);
        MendelianRandomization.validateConfidence(confidenceLevel);
        double[] correlation = correlation(
            alleleAlignedCorrelation, values.size());
        if (backendPolicy == null) {
            throw new IllegalArgumentException("backendPolicy is required");
        }
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            int size = values.size();
            double[] covariance = outcomeCovariance(values, correlation, true);
            CholeskyFactor covarianceFactor = positiveDefiniteFactor(
                covariance, size, backend);
            double[] design = new double[size * 2];
            double[] orientedOutcome = outcome(values, true);
            for (int index = 0; index < size; index++) {
                design[index * 2] = 1.0;
                design[index * 2 + 1] = Math.abs(
                    values.get(index).exposureEffect());
            }
            double[] inverseDesign = covarianceFactor.solve(design, 2);
            double[] information = MatrixOps.transposeMultiply(
                backend, design, size, 2, inverseDesign, 2);
            CholeskyFactor informationFactor;
            try {
                informationFactor = backend.dpotrf(information, 2);
            } catch (IllegalArgumentException | IllegalStateException exception) {
                throw new IllegalArgumentException(
                    "generalized MR-Egger design is rank deficient", exception);
            }
            double[] inverseOutcome = covarianceFactor.solve(orientedOutcome);
            double[] rightSide = new double[2];
            backend.dgemv(MatrixTranspose.TRANSPOSE, size, 2,
                1.0, design, inverseOutcome, 0.0, rightSide);
            double[] coefficients = informationFactor.solve(rightSide);
            double[] fitted = MatrixOps.multiply(
                backend, design, size, 2, coefficients);
            double[] residual = MatrixOps.subtract(orientedOutcome, fitted);
            double[] inverseResidual = covarianceFactor.solve(residual);
            double q = backend.ddot(size, residual, 0, 1,
                inverseResidual, 0, 1);
            int degreesOfFreedom = size - 2;
            double dispersion = Math.max(1.0, q / degreesOfFreedom);
            double[] coefficientCovariance = informationFactor.solve(
                MatrixOps.identity(2), 2);
            double interceptStandardError = Math.sqrt(
                dispersion * coefficientCovariance[0]);
            double slopeStandardError = Math.sqrt(
                dispersion * coefficientCovariance[3]);
            MrEstimate slope = MendelianRandomization.estimate(
                MrMethod.MR_EGGER_GENERALIZED, coefficients[1],
                slopeStandardError, confidenceLevel, q, degreesOfFreedom,
                dispersion, size);
            double interceptStatistic = coefficients[0] / interceptStandardError;
            double critical = MendelianRandomization.normalCritical(confidenceLevel);
            MrEggerResult result = new MrEggerResult(
                slope, coefficients[0], interceptStandardError,
                interceptStatistic,
                MendelianRandomization.normalPValue(interceptStatistic),
                coefficients[0] - critical * interceptStandardError,
                coefficients[0] + critical * interceptStandardError,
                Double.NaN);
            return new CorrelatedMrEggerResult(result, context.provenance());
        }
    }

    private static double[] correlation(double[][] matrix, int size) {
        double[] result = MatrixOps.rowMajor(matrix, size);
        if (matrix[0].length != size) {
            throw new IllegalArgumentException(
                "LD correlation must be square and match instrument count");
        }
        for (int row = 0; row < size; row++) {
            if (Math.abs(result[row * size + row] - 1.0) > 1e-10) {
                throw new IllegalArgumentException(
                    "LD correlation diagonal must equal one");
            }
            for (int column = 0; column < row; column++) {
                double first = result[row * size + column];
                double second = result[column * size + row];
                if (Math.abs(first - second) > 1e-12
                        * Math.max(1.0, Math.max(Math.abs(first), Math.abs(second)))) {
                    throw new IllegalArgumentException(
                        "LD correlation must be symmetric");
                }
                if (Math.abs(first) > 1.0 + 1e-12) {
                    throw new IllegalArgumentException(
                        "LD correlations must lie in [-1, 1]");
                }
            }
        }
        return result;
    }

    private static double[] outcomeCovariance(
            List<HarmonizedInstrument> instruments,
            double[] correlation,
            boolean orientToPositiveExposure) {
        int size = instruments.size();
        double[] result = new double[size * size];
        for (int row = 0; row < size; row++) {
            HarmonizedInstrument rowInstrument = instruments.get(row);
            double rowSign = orientToPositiveExposure
                ? Math.copySign(1.0, rowInstrument.exposureEffect()) : 1.0;
            for (int column = 0; column < size; column++) {
                HarmonizedInstrument columnInstrument = instruments.get(column);
                double columnSign = orientToPositiveExposure
                    ? Math.copySign(1.0, columnInstrument.exposureEffect()) : 1.0;
                result[row * size + column] = correlation[row * size + column]
                    * rowInstrument.outcomeStandardError()
                    * columnInstrument.outcomeStandardError()
                    * rowSign * columnSign;
            }
        }
        return result;
    }

    private static double[] exposure(
            List<HarmonizedInstrument> instruments, boolean orient) {
        double[] result = new double[instruments.size()];
        for (int index = 0; index < result.length; index++) {
            result[index] = orient
                ? Math.abs(instruments.get(index).exposureEffect())
                : instruments.get(index).exposureEffect();
        }
        return result;
    }

    private static double[] outcome(
            List<HarmonizedInstrument> instruments, boolean orient) {
        double[] result = new double[instruments.size()];
        for (int index = 0; index < result.length; index++) {
            HarmonizedInstrument instrument = instruments.get(index);
            result[index] = orient
                ? Math.copySign(1.0, instrument.exposureEffect())
                    * instrument.outcomeEffect()
                : instrument.outcomeEffect();
        }
        return result;
    }

    private static CholeskyFactor positiveDefiniteFactor(
            double[] covariance, int size, ComputeBackend backend) {
        try {
            return backend.dpotrf(covariance, size);
        } catch (IllegalArgumentException | IllegalStateException exception) {
            throw new IllegalArgumentException(
                "LD-implied outcome covariance must be positive definite",
                exception);
        }
    }
}
