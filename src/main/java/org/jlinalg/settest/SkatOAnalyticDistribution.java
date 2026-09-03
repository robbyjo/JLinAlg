/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import java.util.List;
import jdistlib.ChiSquare;
import jdistlib.accelerator.ComputeBackend;

/** Moment-matched one-dimensional SKAT-O minimum-p calibration. */
final class SkatOAnalyticDistribution {
    private static final double CHI_SQUARE_UPPER = 40.0;
    private static final double INTEGRATION_TOLERANCE = 1e-12;
    private static final int MAXIMUM_DEPTH = 20;
    private static final double HALF_NORMAL_SCALE = Math.sqrt(2.0 / Math.PI);

    private SkatOAnalyticDistribution() { }

    static double adjustedPValue(
            double minimumP, double[] rho,
            List<SkatOResult.Component> components,
            SetTestScoreState base, ComputeBackend backend) {
        if (!(minimumP > 0.0 && minimumP <= 1.0)
                || rho.length != components.size())
            throw new IllegalArgumentException(
                "SKAT-O components and minimum p-value are invalid");
        double[] qValues = componentQuantiles(
            minimumP, components);
        if (qValues == null)
            return bonferroni(minimumP, rho.length);
        double[] information = base.informationView();
        int variants = base.variants();
        double[] rowSums = new double[variants];
        double informationSum = 0.0;
        for (int row = 0; row < variants; row++) {
            for (int column = 0; column < variants; column++)
                rowSums[row] += information[row * variants + column];
            informationSum += rowSums[row];
        }
        if (!(informationSum > 0.0) || !Double.isFinite(informationSum))
            return bonferroni(minimumP, rho.length);

        double[] rankOne = new double[information.length];
        double[] orthogonal = new double[information.length];
        double rowSumSquares = 0.0;
        for (int row = 0; row < variants; row++) {
            rowSumSquares += rowSums[row] * rowSums[row];
            for (int column = 0; column < variants; column++) {
                int index = row * variants + column;
                rankOne[index] = rowSums[row] * rowSums[column]
                    / informationSum;
                orthogonal[index] = information[index] - rankOne[index];
            }
        }
        double[] eigenvalues = SetTests.eigenvalues(
            orthogonal, variants, backend);
        if (eigenvalues.length == 0)
            return bonferroni(minimumP, rho.length);
        Moments moments = moments(eigenvalues);
        double cross = 0.0;
        for (int index = 0; index < orthogonal.length; index++)
            cross += rankOne[index] * orthogonal[index];
        double variance = 2.0 * moments.sumSquares + 4.0 * cross;
        if (!(moments.degrees > 0.0) || !(variance > 0.0)
                || !Double.isFinite(variance))
            return bonferroni(minimumP, rho.length);

        double orthogonalBurden = rowSumSquares / informationSum;
        double[] tau = new double[rho.length];
        double upper = Math.sqrt(CHI_SQUARE_UPPER);
        for (int index = 0; index < rho.length; index++) {
            tau[index] = rho[index] * informationSum
                + (1.0 - rho[index]) * orthogonalBurden;
            if (rho[index] == 1.0 && tau[index] > 0.0
                    && qValues[index] > 0.0)
                upper = Math.min(upper,
                    Math.sqrt(qValues[index] / tau[index]));
        }
        Integrand integrand = new Integrand(
            rho, qValues, tau, moments.sum, variance, moments.degrees);
        double middle = 0.5 * upper;
        double left = integrand.value(0.0);
        double center = integrand.value(middle);
        double right = integrand.value(upper);
        double whole = simpson(0.0, upper, left, center, right);
        double cumulative = adaptiveSimpson(integrand, 0.0, upper,
            left, center, right, whole, INTEGRATION_TOLERANCE,
            MAXIMUM_DEPTH);
        double adjusted = Math.min(
            1.0 - Math.max(0.0, Math.min(1.0, cumulative)),
            bonferroni(minimumP, rho.length));
        return Math.max(Double.MIN_VALUE, Math.min(1.0, adjusted));
    }

    private static double[] componentQuantiles(
            double minimumP, List<SkatOResult.Component> components) {
        double[] result = new double[components.size()];
        for (int index = 0; index < result.length; index++) {
            double[] eigenvalues = components.get(index)
                .result().eigenvalues();
            Moments moments = moments(eigenvalues);
            if (!(moments.degrees > 0.0)
                    || !Double.isFinite(moments.degrees))
                return null;
            double quantile = ChiSquare.quantile(
                minimumP, moments.degrees, false, false);
            result[index] = (quantile - moments.degrees)
                / Math.sqrt(2.0 * moments.degrees)
                * Math.sqrt(2.0 * moments.sumSquares) + moments.sum;
            if (!Double.isFinite(result[index])) return null;
        }
        return result;
    }

    private static Moments moments(double[] eigenvalues) {
        double sum = 0.0;
        double sumSquares = 0.0;
        double sumFourth = 0.0;
        for (double eigenvalue : eigenvalues) {
            double square = eigenvalue * eigenvalue;
            sum += eigenvalue;
            sumSquares += square;
            sumFourth += square * square;
        }
        double degrees = sumFourth > 0.0
            ? sumSquares * sumSquares / sumFourth : Double.NaN;
        return new Moments(sum, sumSquares, degrees);
    }

    private static double adaptiveSimpson(
            Integrand function, double left, double right,
            double leftValue, double centerValue, double rightValue,
            double whole, double tolerance, int depth) {
        double center = 0.5 * (left + right);
        double leftCenter = 0.5 * (left + center);
        double rightCenter = 0.5 * (center + right);
        double leftCenterValue = function.value(leftCenter);
        double rightCenterValue = function.value(rightCenter);
        double leftIntegral = simpson(left, center,
            leftValue, leftCenterValue, centerValue);
        double rightIntegral = simpson(center, right,
            centerValue, rightCenterValue, rightValue);
        double refined = leftIntegral + rightIntegral;
        if (depth == 0 || Math.abs(refined - whole) <= 15.0 * tolerance)
            return refined + (refined - whole) / 15.0;
        return adaptiveSimpson(function, left, center,
                leftValue, leftCenterValue, centerValue,
                leftIntegral, tolerance * 0.5, depth - 1)
            + adaptiveSimpson(function, center, right,
                centerValue, rightCenterValue, rightValue,
                rightIntegral, tolerance * 0.5, depth - 1);
    }

    private static double simpson(
            double left, double right,
            double leftValue, double centerValue, double rightValue) {
        return (right - left)
            * (leftValue + 4.0 * centerValue + rightValue) / 6.0;
    }

    private static double bonferroni(double minimumP, int tests) {
        return Math.min(1.0, minimumP * tests);
    }

    private record Moments(double sum, double sumSquares, double degrees) { }

    private record Integrand(
            double[] rho, double[] qValues, double[] tau,
            double mean, double variance, double degrees) {
        double value(double rootChiSquare) {
            double chiSquare = rootChiSquare * rootChiSquare;
            double threshold = Double.POSITIVE_INFINITY;
            for (int index = 0; index < rho.length; index++) {
                double numerator = qValues[index] - tau[index] * chiSquare;
                if (rho[index] == 1.0) {
                    if (numerator < 0.0) return 0.0;
                    continue;
                }
                threshold = Math.min(
                    threshold, numerator / (1.0 - rho[index]));
            }
            double transformed = (threshold - mean) / Math.sqrt(variance)
                * Math.sqrt(2.0 * degrees) + degrees;
            double probability = transformed <= 0.0 ? 0.0
                : ChiSquare.cumulative(transformed, degrees, true, false);
            return probability * HALF_NORMAL_SCALE
                * Math.exp(-0.5 * chiSquare);
        }
    }
}
