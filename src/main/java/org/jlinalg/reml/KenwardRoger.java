/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.reml;

import java.util.List;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.internal.MatrixOps;

/** Internal Kenward-Roger calculations for a linear covariance model. */
final class KenwardRoger {
    private static final double MINIMUM_DENOMINATOR = 1e-14;

    private KenwardRoger() { }

    static Result calculate(
            double[] fixedCovariance,
            double[] inverseCovariance,
            double[] inverseCovarianceFixed,
            double[] logVarianceCovariance,
            double[] variances,
            List<VarianceComponent> components,
            int observations,
            int coefficients,
            ComputeBackend backend) {
        int count = components.size();
        double[] rawVarianceCovariance = new double[count * count];
        for (int row = 0; row < count; row++) {
            for (int column = 0; column < count; column++) {
                rawVarianceCovariance[row * count + column] = variances[row]
                    * logVarianceCovariance[row * count + column]
                    * variances[column];
            }
        }

        double[][] covarianceDerivativeFixed = new double[count][];
        double[][] precisionDerivative = new double[count][];
        double[][] inversePrecisionDerivative = new double[count][];
        for (int component = 0; component < count; component++) {
            covarianceDerivativeFixed[component] = MatrixOps.multiply(
                backend, components.get(component).covarianceView(),
                observations, observations, inverseCovarianceFixed,
                coefficients);
            precisionDerivative[component] = MatrixOps.transposeMultiply(
                backend, covarianceDerivativeFixed[component], observations,
                coefficients, inverseCovarianceFixed, coefficients);
            for (int index = 0;
                    index < precisionDerivative[component].length; index++) {
                precisionDerivative[component][index] =
                    -precisionDerivative[component][index];
            }
            symmetrize(precisionDerivative[component], coefficients);
            inversePrecisionDerivative[component] = MatrixOps.multiply(
                backend, inverseCovariance, observations, observations,
                covarianceDerivativeFixed[component], coefficients);
        }

        double[] adjustmentCore = new double[coefficients * coefficients];
        for (int left = 0; left < count; left++) {
            double[] leftTimesCovariance = MatrixOps.multiply(
                backend, precisionDerivative[left], coefficients, coefficients,
                fixedCovariance, coefficients);
            for (int right = 0; right < count; right++) {
                double weight = rawVarianceCovariance[left * count + right];
                if (weight == 0.0) {
                    continue;
                }
                double[] secondDerivative = MatrixOps.transposeMultiply(
                    backend, covarianceDerivativeFixed[left], observations,
                    coefficients, inversePrecisionDerivative[right],
                    coefficients);
                double[] firstDerivativeProduct = MatrixOps.multiply(
                    backend, leftTimesCovariance, coefficients, coefficients,
                    precisionDerivative[right], coefficients);
                for (int index = 0; index < adjustmentCore.length; index++) {
                    adjustmentCore[index] += weight
                        * (secondDerivative[index]
                            - firstDerivativeProduct[index]);
                }
            }
        }
        symmetrize(adjustmentCore, coefficients);
        double[] leftAdjustment = MatrixOps.multiply(
            backend, fixedCovariance, coefficients, coefficients,
            adjustmentCore, coefficients);
        double[] adjustment = MatrixOps.multiply(
            backend, leftAdjustment, coefficients, coefficients,
            fixedCovariance, coefficients);
        double[] adjustedCovariance = fixedCovariance.clone();
        for (int index = 0; index < adjustedCovariance.length; index++) {
            adjustedCovariance[index] += 2.0 * adjustment[index];
        }
        symmetrize(adjustedCovariance, coefficients);

        double[] degreesOfFreedom = coefficientDegreesOfFreedom(
            fixedCovariance, precisionDerivative, rawVarianceCovariance,
            coefficients, count);
        return new Result(adjustedCovariance, degreesOfFreedom);
    }

    private static double[] coefficientDegreesOfFreedom(
            double[] fixedCovariance,
            double[][] precisionDerivative,
            double[] rawVarianceCovariance,
            int coefficients,
            int componentCount) {
        double[] result = new double[coefficients];
        double[] relativeDerivatives = new double[componentCount];
        for (int coefficient = 0;
                coefficient < coefficients; coefficient++) {
            double variance = fixedCovariance[
                coefficient * coefficients + coefficient];
            if (!(variance > 0.0)) {
                result[coefficient] = Double.NaN;
                continue;
            }
            for (int component = 0;
                    component < componentCount; component++) {
                double value = 0.0;
                double[] derivative = precisionDerivative[component];
                for (int row = 0; row < coefficients; row++) {
                    double left = fixedCovariance[
                        coefficient * coefficients + row] / variance;
                    for (int column = 0; column < coefficients; column++) {
                        value += left
                            * derivative[row * coefficients + column]
                            * fixedCovariance[
                                column * coefficients + coefficient];
                    }
                }
                relativeDerivatives[component] = value;
            }
            double moment = 0.0;
            for (int left = 0; left < componentCount; left++) {
                for (int right = 0; right < componentCount; right++) {
                    moment += relativeDerivatives[left]
                        * rawVarianceCovariance[left * componentCount + right]
                        * relativeDerivatives[right];
                }
            }
            result[coefficient] = moment > MINIMUM_DENOMINATOR
                    && Double.isFinite(moment)
                ? 2.0 / moment : Double.NaN;
        }
        return result;
    }

    private static void symmetrize(double[] matrix, int dimension) {
        for (int row = 0; row < dimension; row++) {
            for (int column = row + 1; column < dimension; column++) {
                double value = 0.5 * (matrix[row * dimension + column]
                    + matrix[column * dimension + row]);
                matrix[row * dimension + column] = value;
                matrix[column * dimension + row] = value;
            }
        }
    }

    record Result(double[] adjustedCovariance, double[] degreesOfFreedom) { }
}
