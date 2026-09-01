/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.List;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.internal.LeastSquaresSolver;

/** Package-private weighted meta-analytic linear algebra. */
final class MetaMath {
    private static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);
    private MetaMath() { }

    static Data data(List<MetaStudy> studies) {
        if (studies == null || studies.size() < 2)
            throw new IllegalArgumentException("at least two studies are required");
        double[] effects = new double[studies.size()];
        double[] variances = new double[studies.size()];
        for (int index = 0; index < studies.size(); index++) {
            MetaStudy study = studies.get(index);
            if (study == null)
                throw new IllegalArgumentException("studies must not contain null");
            effects[index] = study.effectSize();
            variances[index] = study.variance();
        }
        return new Data(effects, variances);
    }

    static Fit fit(Data data, double[] design, int columns,
                   double tauSquared, ComputeBackend backend) {
        int rows = data.effects().length;
        if (columns < 1 || rows <= columns || design == null
                || design.length != rows * columns)
            throw new IllegalArgumentException(
                "meta-regression needs more studies than full-rank coefficients");
        double[] weights = new double[rows];
        double[] weightedDesign = new double[design.length];
        double[] weightedResponse = new double[rows];
        double sumLogVariance = 0.0;
        for (int row = 0; row < rows; row++) {
            double totalVariance = data.variances()[row] + tauSquared;
            double rootWeight = 1.0 / Math.sqrt(totalVariance);
            weights[row] = rootWeight * rootWeight;
            weightedResponse[row] = data.effects()[row] * rootWeight;
            sumLogVariance += Math.log(totalVariance);
            for (int column = 0; column < columns; column++)
                weightedDesign[row * columns + column] =
                    design[row * columns + column] * rootWeight;
        }
        LeastSquaresSolver.Solution solved = LeastSquaresSolver.solve(
            weightedDesign, weightedResponse, rows, columns, false, backend);
        double[] beta = solved.coefficients();
        double qe = 0.0;
        for (int row = 0; row < rows; row++) {
            double fitted = 0.0;
            for (int column = 0; column < columns; column++)
                fitted += design[row * columns + column] * beta[column];
            double residual = data.effects()[row] - fitted;
            qe += weights[row] * residual * residual;
        }
        double[] information = crossProduct(design, weights, rows, columns, 1);
        CholeskyFactor factor = backend.dpotrf(information, columns);
        double restrictedObjective = sumLogVariance + factor.logDeterminant()
            + qe + (rows - columns) * LOG_TWO_PI;
        return new Fit(beta, solved.unscaledCovariance(), weights, qe,
            restrictedObjective);
    }

    static double estimateTauSquared(
            Data data, double[] design, int columns,
            MetaAnalysisOptions options, ComputeBackend backend) {
        if (options.method() == MetaAnalysisMethod.FIXED_EFFECT) return 0.0;
        return switch (options.tauSquaredEstimator()) {
            case DERSIMONIAN_LAIRD -> dersimonianLaird(
                data, design, columns, backend);
            case PAULE_MANDEL -> pauleMandel(
                data, design, columns, options, backend);
            case REML -> reml(data, design, columns, options, backend);
        };
    }

    private static double dersimonianLaird(
            Data data, double[] design, int columns, ComputeBackend backend) {
        Fit fixed = fit(data, design, columns, 0.0, backend);
        int rows = data.effects().length;
        double[] second = crossProduct(
            design, fixed.weights(), rows, columns, 2);
        double trace = 0.0;
        for (int row = 0; row < columns; row++)
            for (int column = 0; column < columns; column++)
                trace += fixed.covariance()[row * columns + column]
                    * second[column * columns + row];
        double sumWeights = 0.0;
        for (double weight : fixed.weights()) sumWeights += weight;
        double denominator = sumWeights - trace;
        return denominator > 0.0
            ? Math.max(0.0, (fixed.qe() - (rows - columns)) / denominator)
            : 0.0;
    }

    private static double pauleMandel(
            Data data, double[] design, int columns,
            MetaAnalysisOptions options, ComputeBackend backend) {
        double target = data.effects().length - columns;
        if (fit(data, design, columns, 0.0, backend).qe() <= target) return 0.0;
        double upper = startingUpper(data);
        while (fit(data, design, columns, upper, backend).qe() > target
                && upper < 1e12) upper *= 4.0;
        double lower = 0.0;
        for (int iteration = 0; iteration < options.maximumIterations(); iteration++) {
            double middle = 0.5 * (lower + upper);
            if (fit(data, design, columns, middle, backend).qe() > target)
                lower = middle;
            else upper = middle;
            if (upper - lower <= options.tolerance() * Math.max(1.0, upper))
                break;
        }
        return 0.5 * (lower + upper);
    }

    private static double reml(
            Data data, double[] design, int columns,
            MetaAnalysisOptions options, ComputeBackend backend) {
        double upper = startingUpper(data);
        double previous = fit(data, design, columns, 0.0, backend)
            .restrictedObjective();
        double atUpper = fit(data, design, columns, upper, backend)
            .restrictedObjective();
        while (atUpper < previous && upper < 1e12) {
            previous = atUpper;
            upper *= 4.0;
            atUpper = fit(data, design, columns, upper, backend)
                .restrictedObjective();
        }
        double left = 0.0;
        double right = upper;
        double ratio = (Math.sqrt(5.0) - 1.0) / 2.0;
        double first = right - ratio * (right - left);
        double second = left + ratio * (right - left);
        double firstValue = fit(data, design, columns, first, backend)
            .restrictedObjective();
        double secondValue = fit(data, design, columns, second, backend)
            .restrictedObjective();
        for (int iteration = 0; iteration < options.maximumIterations(); iteration++) {
            if (firstValue < secondValue) {
                right = second; second = first; secondValue = firstValue;
                first = right - ratio * (right - left);
                firstValue = fit(data, design, columns, first, backend)
                    .restrictedObjective();
            } else {
                left = first; first = second; firstValue = secondValue;
                second = left + ratio * (right - left);
                secondValue = fit(data, design, columns, second, backend)
                    .restrictedObjective();
            }
            if (right - left <= options.tolerance() * Math.max(1.0, right))
                break;
        }
        double candidate = 0.5 * (left + right);
        double atZero = fit(data, design, columns, 0.0, backend)
            .restrictedObjective();
        return atZero <= fit(data, design, columns, candidate, backend)
            .restrictedObjective() ? 0.0 : candidate;
    }

    private static double startingUpper(Data data) {
        double mean = 0.0;
        for (double value : data.effects()) mean += value;
        mean /= data.effects().length;
        double variance = 0.0;
        for (double value : data.effects()) variance += (value - mean) * (value - mean);
        return Math.max(1e-8, variance / Math.max(1, data.effects().length - 1));
    }

    private static double[] crossProduct(
            double[] design, double[] weights, int rows, int columns,
            int weightPower) {
        double[] result = new double[columns * columns];
        for (int row = 0; row < rows; row++) {
            double weight = weightPower == 1
                ? weights[row] : weights[row] * weights[row];
            for (int left = 0; left < columns; left++)
                for (int right = 0; right < columns; right++)
                    result[left * columns + right] += weight
                        * design[row * columns + left]
                        * design[row * columns + right];
        }
        return result;
    }

    record Data(double[] effects, double[] variances) { }
    record Fit(double[] beta, double[] covariance, double[] weights,
               double qe, double restrictedObjective) { }
}
