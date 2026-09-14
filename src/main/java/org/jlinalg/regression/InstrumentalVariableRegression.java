/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.regression;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import jdistlib.F;
import jdistlib.Gamma;
import jdistlib.Normal;
import jdistlib.T;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.SingularValueDecomposition;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.inference.StatisticDistribution;
import org.jlinalg.internal.LeastSquaresSolver;
import org.jlinalg.internal.LeastSquaresSolver.Solution;
import org.jlinalg.internal.MatrixOps;

/**
 * Individual-level linear instrumental-variable regression by two-stage least
 * squares (2SLS).
 *
 * <p>The structural design is {@code [exogenous, endogenous]} and the instrument
 * design is {@code [exogenous, excludedInstruments]}. Callers include an
 * intercept explicitly when wanted. First-stage projections use the shared
 * equilibrated pivoted-QR solver; SVD supplies a scale-neutral identification
 * condition diagnostic. The selected covariance always uses structural
 * residuals {@code y - X beta}, never residuals from an ordinary regression on
 * first-stage fitted predictors.</p>
 */
public final class InstrumentalVariableRegression {
    private InstrumentalVariableRegression() { }

    /** Fits conventional homoskedastic 2SLS on the preferred backend. */
    public static InstrumentalVariableResult fit(
            double[] response,
            double[][] exogenous,
            double[][] endogenous,
            double[][] excludedInstruments) {
        return fit(response, exogenous, endogenous, excludedInstruments,
            null, InstrumentalVariableOptions.defaults(),
            BackendPolicy.PREFERRED);
    }

    /** Fits 2SLS with caller-selected covariance and backend. */
    public static InstrumentalVariableResult fit(
            double[] response,
            double[][] exogenous,
            double[][] endogenous,
            double[][] excludedInstruments,
            InstrumentalVariableOptions options,
            BackendPolicy backendPolicy) {
        return fit(response, exogenous, endogenous, excludedInstruments,
            null, options, backendPolicy);
    }

    /** Fits 2SLS, optionally using cluster-robust inference. */
    public static InstrumentalVariableResult fit(
            double[] response,
            double[][] exogenous,
            double[][] endogenous,
            double[][] excludedInstruments,
            int[] clusters,
            InstrumentalVariableOptions options,
            BackendPolicy backendPolicy) {
        if (response == null || response.length < 1) {
            throw new IllegalArgumentException("response is required");
        }
        if (options == null || backendPolicy == null) {
            throw new IllegalArgumentException("options and backend policy are required");
        }
        MatrixOps.requireFinite(response, "response");
        int observations = response.length;
        MatrixData exogenousData = matrix(exogenous, observations,
            "exogenous design", true);
        MatrixData endogenousData = matrix(endogenous, observations,
            "endogenous design", false);
        MatrixData excludedData = matrix(excludedInstruments, observations,
            "excluded instruments", false);
        int exogenousColumns = exogenousData.columns();
        int endogenousColumns = endogenousData.columns();
        int excludedColumns = excludedData.columns();
        int parameters = exogenousColumns + endogenousColumns;
        int instruments = exogenousColumns + excludedColumns;
        if (excludedColumns < endogenousColumns) {
            throw new IllegalArgumentException(
                "2SLS is underidentified: excluded instruments must be at least "
                + "the number of endogenous regressors");
        }
        if (observations <= parameters || observations <= instruments) {
            throw new IllegalArgumentException(
                "2SLS requires residual degrees of freedom in both structural and first-stage models");
        }
        if (options.covariance().requiresClusters() != (clusters != null)) {
            throw new IllegalArgumentException(options.covariance().requiresClusters()
                ? "cluster covariance requires one cluster identifier per observation"
                : "cluster identifiers require a cluster covariance estimator");
        }
        int clusterCount = clusters == null ? 0
            : validateClusters(clusters, observations, parameters, instruments);

        double[] structural = concatenate(exogenousData.values(),
            exogenousColumns, endogenousData.values(), endogenousColumns,
            observations);
        double[] instrumentDesign = concatenate(exogenousData.values(),
            exogenousColumns, excludedData.values(), excludedColumns,
            observations);

        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            requireFullRank(structural, observations, parameters,
                "structural design", backend);

            Solution[] firstStages = new Solution[endogenousColumns];
            double[] instrumentedEndogenous = new double[observations * endogenousColumns];
            double[] instrumentedStructural = new double[observations * parameters];
            for (int row = 0; row < observations; row++) {
                System.arraycopy(exogenousData.values(), row * exogenousColumns,
                    instrumentedStructural, row * parameters, exogenousColumns);
            }
            for (int column = 0; column < endogenousColumns; column++) {
                double[] target = column(endogenousData.values(), observations,
                    endogenousColumns, column);
                try {
                    firstStages[column] = LeastSquaresSolver.solve(
                        instrumentDesign, target, observations, instruments,
                        false, backend);
                } catch (IllegalArgumentException exception) {
                    throw new IllegalArgumentException(
                        "instrument design is rank deficient; remove redundant exogenous or instrument columns",
                        exception);
                }
                double[] fitted = MatrixOps.multiply(backend, instrumentDesign,
                    observations, instruments,
                    firstStages[column].coefficients());
                for (int row = 0; row < observations; row++) {
                    instrumentedEndogenous[row * endogenousColumns + column] = fitted[row];
                    instrumentedStructural[row * parameters
                        + exogenousColumns + column] = fitted[row];
                }
            }

            double conditionNumber = identificationConditionNumber(structural,
                instrumentedStructural, observations, parameters, backend);
            Solution secondStage;
            try {
                secondStage = LeastSquaresSolver.solve(instrumentedStructural,
                    response, observations, parameters, false, backend);
            } catch (IllegalArgumentException exception) {
                throw new IllegalArgumentException(
                    "projected structural design is rank deficient; endogenous regressors "
                    + "are not identified by the supplied instruments", exception);
            }
            double[] coefficients = secondStage.coefficients();
            double[] fitted = MatrixOps.multiply(backend, structural,
                observations, parameters, coefficients);
            double[] residuals = MatrixOps.subtract(response, fitted);
            double residualSumOfSquares = sumSquares(residuals);
            int residualDegreesOfFreedom = observations - parameters;
            double residualVariance = residualSumOfSquares
                / residualDegreesOfFreedom;

            double[] covariance = covariance(secondStage.unscaledCovariance(),
                instrumentedStructural, residuals, observations, parameters,
                clusters, clusterCount, options.covariance());
            double inferenceDegreesOfFreedom = switch (options.covariance()) {
                case HOMOSKEDASTIC -> residualDegreesOfFreedom;
                case HC0, HC1 -> Double.POSITIVE_INFINITY;
                case CLUSTER_CR0, CLUSTER_CR1 -> clusterCount - parameters;
            };
            if (options.covariance() == InstrumentalVariableCovariance.HOMOSKEDASTIC) {
                covariance = secondStage.unscaledCovariance().clone();
                scale(covariance, residualVariance);
            }
            double[] standardErrors = diagonalStandardErrors(covariance,
                parameters);
            AssociationStatistics association = association(coefficients,
                standardErrors, inferenceDegreesOfFreedom,
                options.covariance());
            double[][] intervals = intervals(coefficients, standardErrors,
                inferenceDegreesOfFreedom, options.confidenceLevel());
            List<InstrumentStrengthDiagnostic> diagnostics = diagnostics(
                endogenousData.values(), endogenousColumns,
                exogenousData.values(), exogenousColumns,
                instrumentDesign, instruments, firstStages,
                instrumentedEndogenous, observations, clusters, clusterCount,
                options.covariance(), backend);
            return new InstrumentalVariableResult(coefficients, covariance,
                standardErrors, intervals[0], intervals[1], fitted, residuals,
                instrumentedEndogenous, observations, exogenousColumns,
                endogenousColumns, excludedColumns, secondStage.rank(),
                firstStages[0].rank(), residualDegreesOfFreedom,
                inferenceDegreesOfFreedom, clusterCount,
                residualSumOfSquares, residualVariance, conditionNumber,
                diagnostics, association, options.covariance(),
                context.provenance());
        }
    }

    private static List<InstrumentStrengthDiagnostic> diagnostics(
            double[] endogenous,
            int endogenousColumns,
            double[] exogenous,
            int exogenousColumns,
            double[] instrumentDesign,
            int instruments,
            Solution[] firstStages,
            double[] instrumentedEndogenous,
            int observations,
            int[] clusters,
            int clusterCount,
            InstrumentalVariableCovariance estimator,
            ComputeBackend backend) {
        List<InstrumentStrengthDiagnostic> result = new ArrayList<>(
            endogenousColumns);
        int excludedColumns = instruments - exogenousColumns;
        for (int column = 0; column < endogenousColumns; column++) {
            double[] target = column(endogenous, observations,
                endogenousColumns, column);
            double[] unrestrictedFitted = column(instrumentedEndogenous,
                observations, endogenousColumns, column);
            double[] unrestrictedResiduals = MatrixOps.subtract(target,
                unrestrictedFitted);
            double unrestrictedRss = sumSquares(unrestrictedResiduals);
            double restrictedRss;
            if (exogenousColumns == 0) {
                restrictedRss = sumSquares(target);
            } else {
                Solution restricted = LeastSquaresSolver.solve(exogenous,
                    target, observations, exogenousColumns, false, backend);
                restrictedRss = sumSquares(MatrixOps.subtract(target,
                    MatrixOps.multiply(backend, exogenous, observations,
                        exogenousColumns, restricted.coefficients())));
            }
            double improvement = Math.max(0.0,
                restrictedRss - unrestrictedRss);
            double residualTolerance = 128.0 * Math.ulp(1.0)
                * Math.max(Double.MIN_NORMAL,
                    Math.max(restrictedRss, sumSquares(target)));
            boolean perfectFirstStage = improvement > residualTolerance
                && unrestrictedRss <= residualTolerance;
            double partialRSquared = restrictedRss == 0.0
                ? Double.NaN : Math.min(1.0, improvement / restrictedRss);
            double classicalF = perfectFirstStage
                ? Double.POSITIVE_INFINITY
                : (improvement / excludedColumns)
                    / (unrestrictedRss / (observations - instruments));

            double wald;
            double effectiveF;
            double pValue;
            double denominatorDegreesOfFreedom;
            StatisticDistribution distribution;
            if (estimator == InstrumentalVariableCovariance.HOMOSKEDASTIC) {
                wald = excludedColumns * classicalF;
                effectiveF = classicalF;
                denominatorDegreesOfFreedom = observations - instruments;
                pValue = F.cumulative(classicalF, excludedColumns,
                    denominatorDegreesOfFreedom, false, false);
                distribution = StatisticDistribution.F;
            } else if (perfectFirstStage) {
                wald = Double.POSITIVE_INFINITY;
                effectiveF = Double.POSITIVE_INFINITY;
                pValue = 0.0;
                if (estimator.requiresClusters()) {
                    denominatorDegreesOfFreedom = clusterCount - instruments;
                    distribution = StatisticDistribution.F;
                } else {
                    denominatorDegreesOfFreedom = Double.POSITIVE_INFINITY;
                    distribution = StatisticDistribution.CHI_SQUARE;
                }
            } else {
                double[] firstStageCovariance = covariance(
                    firstStages[column].unscaledCovariance(), instrumentDesign,
                    unrestrictedResiduals, observations, instruments, clusters,
                    clusterCount, estimator);
                wald = excludedWald(firstStages[column].coefficients(),
                    firstStageCovariance, exogenousColumns, excludedColumns,
                    instruments, backend);
                effectiveF = wald / excludedColumns;
                if (estimator.requiresClusters()) {
                    denominatorDegreesOfFreedom = clusterCount - instruments;
                    pValue = F.cumulative(effectiveF, excludedColumns,
                        denominatorDegreesOfFreedom, false, false);
                    distribution = StatisticDistribution.F;
                } else {
                    denominatorDegreesOfFreedom = Double.POSITIVE_INFINITY;
                    pValue = Gamma.cumulative(wald, excludedColumns / 2.0,
                        2.0, false, false);
                    distribution = StatisticDistribution.CHI_SQUARE;
                }
            }
            result.add(new InstrumentStrengthDiagnostic(column,
                partialRSquared, classicalF, wald, effectiveF,
                excludedColumns, denominatorDegreesOfFreedom, pValue,
                distribution));
        }
        return List.copyOf(result);
    }

    private static AssociationStatistics association(
            double[] coefficients,
            double[] standardErrors,
            double degreesOfFreedom,
            InstrumentalVariableCovariance estimator) {
        if (Double.isInfinite(degreesOfFreedom)) {
            return AssociationStatistics.normal(coefficients, standardErrors);
        }
        DegreesOfFreedomMethod method = estimator.requiresClusters()
            ? DegreesOfFreedomMethod.CLUSTER
            : DegreesOfFreedomMethod.RESIDUAL;
        return AssociationStatistics.studentT(coefficients, standardErrors,
            degreesOfFreedom, method);
    }

    private static double[][] intervals(
            double[] coefficients,
            double[] standardErrors,
            double degreesOfFreedom,
            double confidenceLevel) {
        double probability = (1.0 + confidenceLevel) / 2.0;
        double critical = Double.isFinite(degreesOfFreedom)
            ? T.quantile(probability, degreesOfFreedom, true, false)
            : Normal.quantile(probability, 0.0, 1.0, true, false);
        double[] lower = new double[coefficients.length];
        double[] upper = new double[coefficients.length];
        for (int index = 0; index < coefficients.length; index++) {
            double margin = critical * standardErrors[index];
            lower[index] = coefficients[index] - margin;
            upper[index] = coefficients[index] + margin;
        }
        return new double[][] {lower, upper};
    }

    private static double[] covariance(
            double[] bread,
            double[] scoreDesign,
            double[] residuals,
            int observations,
            int parameters,
            int[] clusters,
            int clusterCount,
            InstrumentalVariableCovariance estimator) {
        if (estimator == InstrumentalVariableCovariance.HOMOSKEDASTIC) {
            double[] result = bread.clone();
            scale(result, sumSquares(residuals) / (observations - parameters));
            return result;
        }
        double[] meat = new double[parameters * parameters];
        if (estimator.requiresClusters()) {
            Map<Integer, double[]> scores = new LinkedHashMap<>();
            for (int row = 0; row < observations; row++) {
                double[] score = scores.computeIfAbsent(clusters[row],
                    ignored -> new double[parameters]);
                for (int column = 0; column < parameters; column++) {
                    score[column] += scoreDesign[row * parameters + column]
                        * residuals[row];
                }
            }
            for (double[] score : scores.values()) addOuter(meat, score,
                parameters);
            if (estimator == InstrumentalVariableCovariance.CLUSTER_CR1) {
                scale(meat, (double) clusterCount / (clusterCount - 1.0)
                    * (observations - 1.0) / (observations - parameters));
            }
        } else {
            double[] score = new double[parameters];
            for (int row = 0; row < observations; row++) {
                for (int column = 0; column < parameters; column++) {
                    score[column] = scoreDesign[row * parameters + column]
                        * residuals[row];
                }
                addOuter(meat, score, parameters);
            }
            if (estimator == InstrumentalVariableCovariance.HC1) {
                scale(meat, (double) observations
                    / (observations - parameters));
            }
        }
        return sandwich(bread, meat, parameters);
    }

    private static double[] sandwich(
            double[] bread, double[] meat, int parameters) {
        double[] temporary = new double[parameters * parameters];
        double[] result = new double[parameters * parameters];
        for (int row = 0; row < parameters; row++) {
            for (int shared = 0; shared < parameters; shared++) {
                double left = bread[row * parameters + shared];
                for (int column = 0; column < parameters; column++) {
                    temporary[row * parameters + column] += left
                        * meat[shared * parameters + column];
                }
            }
        }
        for (int row = 0; row < parameters; row++) {
            for (int shared = 0; shared < parameters; shared++) {
                double left = temporary[row * parameters + shared];
                for (int column = 0; column < parameters; column++) {
                    result[row * parameters + column] += left
                        * bread[shared * parameters + column];
                }
            }
        }
        for (int row = 0; row < parameters; row++) {
            for (int column = 0; column < row; column++) {
                double value = 0.5 * (result[row * parameters + column]
                    + result[column * parameters + row]);
                result[row * parameters + column] = value;
                result[column * parameters + row] = value;
            }
        }
        return result;
    }

    private static double excludedWald(
            double[] coefficients,
            double[] covariance,
            int offset,
            int count,
            int dimension,
            ComputeBackend backend) {
        double[] correlation = new double[count * count];
        double[] standardized = new double[count];
        for (int row = 0; row < count; row++) {
            double variance = covariance[(offset + row) * dimension
                + offset + row];
            if (!(variance > 0.0) || !Double.isFinite(variance)) {
                throw new IllegalArgumentException(
                    "excluded-instrument Wald covariance is singular");
            }
            double standardError = Math.sqrt(variance);
            standardized[row] = coefficients[offset + row] / standardError;
            for (int column = 0; column < count; column++) {
                double otherVariance = covariance[(offset + column)
                    * dimension + offset + column];
                if (!(otherVariance > 0.0)
                        || !Double.isFinite(otherVariance)) {
                    throw new IllegalArgumentException(
                        "excluded-instrument Wald covariance is singular");
                }
                correlation[row * count + column] = covariance[
                    (offset + row) * dimension + offset + column]
                    / standardError / Math.sqrt(otherVariance);
            }
        }
        try {
            double[] solved = backend.dpotrf(correlation, count)
                .solve(standardized);
            return Math.max(0.0, backend.ddot(count, standardized, 0, 1,
                solved, 0, 1));
        } catch (IllegalArgumentException | IllegalStateException exception) {
            throw new IllegalArgumentException(
                "excluded-instrument Wald covariance is rank deficient",
                exception);
        }
    }

    private static double identificationConditionNumber(
            double[] structuralDesign,
            double[] projectedDesign,
            int observations,
            int parameters,
            ComputeBackend backend) {
        double[] standardized = projectedDesign.clone();
        for (int column = 0; column < parameters; column++) {
            double scale = 0.0;
            for (int row = 0; row < observations; row++) {
                scale = Math.max(scale,
                    Math.abs(structuralDesign[row * parameters + column]));
            }
            if (!(scale > 0.0)) return Double.POSITIVE_INFINITY;
            double sumSquares = 0.0;
            for (int row = 0; row < observations; row++) {
                double value = structuralDesign[row * parameters + column]
                    / scale;
                sumSquares += value * value;
            }
            double norm = scale * Math.sqrt(sumSquares);
            for (int row = 0; row < observations; row++) {
                standardized[row * parameters + column]
                    = projectedDesign[row * parameters + column] / norm;
            }
        }
        SingularValueDecomposition decomposition = backend.dgesvd(
            standardized, observations, parameters);
        double[] singularValues = decomposition.singularValues();
        double smallest = singularValues[singularValues.length - 1];
        double tolerance = Math.max(observations, parameters) * Math.ulp(1.0)
            * singularValues[0];
        if (!(smallest > tolerance)) {
            throw new IllegalArgumentException(
                "projected structural design is rank deficient; endogenous regressors "
                + "are not identified by the supplied instruments");
        }
        return smallest > 0.0 ? singularValues[0] / smallest
            : Double.POSITIVE_INFINITY;
    }

    private static void requireFullRank(
            double[] design,
            int observations,
            int parameters,
            String name,
            ComputeBackend backend) {
        try {
            LeastSquaresSolver.solve(design, new double[observations],
                observations, parameters, false, backend);
        } catch (IllegalArgumentException exception) {
            throw new IllegalArgumentException(name
                + " is rank deficient", exception);
        }
    }

    private static int validateClusters(
            int[] clusters,
            int observations,
            int structuralParameters,
            int firstStageParameters) {
        if (clusters.length != observations) {
            throw new IllegalArgumentException(
                "cluster identifiers must match the response length");
        }
        int count = (int) Arrays.stream(clusters).distinct().count();
        int required = Math.max(structuralParameters, firstStageParameters);
        if (count <= required) {
            throw new IllegalArgumentException(
                "cluster inference requires more independent clusters than "
                + "both structural and first-stage parameters");
        }
        return count;
    }

    private static MatrixData matrix(
            double[][] matrix,
            int rows,
            String name,
            boolean allowZeroColumns) {
        if (matrix == null || matrix.length != rows || matrix[0] == null) {
            throw new IllegalArgumentException(name
                + " must have one rectangular row per response");
        }
        int columns = matrix[0].length;
        if ((!allowZeroColumns && columns < 1)
                || (long) rows * columns > Integer.MAX_VALUE) {
            throw new IllegalArgumentException(name + " has invalid dimensions");
        }
        double[] result = new double[rows * columns];
        for (int row = 0; row < rows; row++) {
            if (matrix[row] == null || matrix[row].length != columns) {
                throw new IllegalArgumentException(name
                    + " must be rectangular");
            }
            for (int column = 0; column < columns; column++) {
                double value = matrix[row][column];
                if (!Double.isFinite(value)) {
                    throw new IllegalArgumentException(name
                        + " must contain only finite values");
                }
                result[row * columns + column] = value;
            }
        }
        return new MatrixData(columns, result);
    }

    private static double[] concatenate(
            double[] left,
            int leftColumns,
            double[] right,
            int rightColumns,
            int rows) {
        int columns = leftColumns + rightColumns;
        double[] result = new double[rows * columns];
        for (int row = 0; row < rows; row++) {
            System.arraycopy(left, row * leftColumns, result,
                row * columns, leftColumns);
            System.arraycopy(right, row * rightColumns, result,
                row * columns + leftColumns, rightColumns);
        }
        return result;
    }

    private static double[] column(
            double[] matrix,
            int rows,
            int columns,
            int selected) {
        double[] result = new double[rows];
        for (int row = 0; row < rows; row++) {
            result[row] = matrix[row * columns + selected];
        }
        return result;
    }

    private static double[] diagonalStandardErrors(
            double[] covariance, int parameters) {
        double[] result = new double[parameters];
        for (int index = 0; index < parameters; index++) {
            double variance = covariance[index * parameters + index];
            if (variance < -1e-12 || !Double.isFinite(variance)) {
                throw new IllegalArgumentException(
                    "IV covariance has an invalid diagonal");
            }
            result[index] = Math.sqrt(Math.max(0.0, variance));
        }
        return result;
    }

    private static void addOuter(
            double[] destination, double[] vector, int dimension) {
        for (int row = 0; row < dimension; row++) {
            for (int column = 0; column < dimension; column++) {
                destination[row * dimension + column] += vector[row]
                    * vector[column];
            }
        }
    }

    private static double sumSquares(double[] values) {
        double result = 0.0;
        for (double value : values) result += value * value;
        return Math.max(0.0, result);
    }

    private static void scale(double[] values, double multiplier) {
        for (int index = 0; index < values.length; index++) {
            values[index] *= multiplier;
        }
    }

    private record MatrixData(int columns, double[] values) { }
}
