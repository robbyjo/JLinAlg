/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import java.util.ArrayList;
import java.util.List;
import jdistlib.ChiSquare;
import jdistlib.Normal;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTranspose;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.model.MissingDataPolicy;

/** Observed-variable covariance-structure maximum-likelihood SEM. */
public final class Sem {
    private Sem() { }

    public static SemFitResult fit(double[][] data, SemModel model) {
        return fit(data, model, SemOptions.defaults(), BackendPolicy.PREFERRED);
    }

    public static SemFitResult fit(
            double[][] data, SemModel model, SemOptions options,
            BackendPolicy backendPolicy) {
        if (data == null || data.length < 2 || model == null
                || options == null || backendPolicy == null) {
            throw new IllegalArgumentException("data, model, options, and backend are required");
        }
        int variables = model.variables().size();
        int complete = 0;
        double[] means = new double[variables];
        for (int row = 0; row < data.length; row++) {
            if (data[row] == null || data[row].length != variables)
                throw new IllegalArgumentException("SEM data must be rectangular and match variables");
            boolean finite = true;
            for (double value : data[row]) finite &= Double.isFinite(value);
            if (finite) {
                complete++;
                for (int variable = 0; variable < variables; variable++)
                    means[variable] += data[row][variable];
            }
            else if (options.missingDataPolicy() == MissingDataPolicy.ERROR)
                throw new IllegalArgumentException("non-finite SEM row: " + row);
        }
        if (complete <= variables)
            throw new IllegalArgumentException("too few complete rows for SEM covariance");
        for (int variable = 0; variable < variables; variable++) means[variable] /= complete;
        double[] covariance = new double[variables * variables];
        for (double[] row : data) {
            boolean finite = true;
            for (double value : row) finite &= Double.isFinite(value);
            if (!finite) continue;
            for (int first = 0; first < variables; first++) {
                double centeredFirst = row[first] - means[first];
                for (int second = 0; second <= first; second++) {
                    double value = centeredFirst * (row[second] - means[second]);
                    covariance[first * variables + second] += value;
                    if (first != second) covariance[second * variables + first] += value;
                }
            }
        }
        for (int index = 0; index < covariance.length; index++) covariance[index] /= complete;
        return fitCovariance(covariance, complete, model, options, backendPolicy);
    }

    public static SemFitResult fitCovariance(
            double[] sampleCovariance, int observations, SemModel model,
            SemOptions options, BackendPolicy backendPolicy) {
        int variables = model.variables().size();
        if (sampleCovariance == null || sampleCovariance.length != variables * variables
                || observations <= variables || options == null || backendPolicy == null) {
            throw new IllegalArgumentException("sample covariance dimensions are invalid");
        }
        int moments = variables * (variables + 1) / 2;
        int degreesOfFreedom = moments - model.freeParameterCount();
        if (degreesOfFreedom < 0)
            throw new IllegalArgumentException("SEM has more free parameters than covariance moments");
        double[] initial = initial(model);
        double[] lower = new double[initial.length];
        double[] upper = new double[initial.length];
        for (int index = 0; index < initial.length; index++) {
            lower[index] = model.varianceParameter(index) ? -20.0 : -20.0;
            upper[index] = 20.0;
        }
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend();
            CholeskyFactor sampleFactor;
            try {
                sampleFactor = backend.dpotrf(sampleCovariance, variables);
            } catch (IllegalArgumentException | IllegalStateException exception) {
                throw new IllegalArgumentException("sample covariance must be positive definite", exception);
            }
            Objective objective = point -> evaluate(
                point, model, sampleCovariance, observations, backend, true);
            Optimum optimum = optimize(initial, lower, upper, objective,
                options.maximumEvaluations(), options.tolerance(), observations);
            Evaluation fitted = evaluate(optimum.parameters(), model,
                sampleCovariance, observations, backend, false);
            double[] rawCovariance = parameterCovariance(
                optimum.parameters(), model, fitted, observations, backend);
            List<SemParameterEstimate> estimates = parameterEstimates(
                optimum.parameters(), rawCovariance, model);
            double logDetSample = sampleFactor.logDeterminant();
            double discrepancy = fitted.logDeterminant() + fitted.trace()
                - logDetSample - variables;
            double chiSquare = Math.max(0.0, observations * discrepancy);
            double pValue = degreesOfFreedom > 0
                ? ChiSquare.cumulative(chiSquare, degreesOfFreedom, false, false)
                : Double.NaN;
            double baselineDiscrepancy = baselineDiscrepancy(
                sampleCovariance, variables, logDetSample);
            double baselineChi = observations * baselineDiscrepancy;
            int baselineDf = variables * (variables - 1) / 2;
            double modelExcess = Math.max(0.0, chiSquare - degreesOfFreedom);
            double baselineExcess = Math.max(modelExcess,
                baselineChi - baselineDf);
            double cfi = baselineExcess > 0.0 ? 1.0 - modelExcess / baselineExcess : 1.0;
            double tli = degreesOfFreedom > 0 && baselineDf > 0 && baselineChi > 0.0
                ? 1.0 - (chiSquare / degreesOfFreedom - 1.0)
                    / (baselineChi / baselineDf - 1.0) : Double.NaN;
            double rmsea = degreesOfFreedom > 0
                ? Math.sqrt(Math.max(0.0,
                    (chiSquare - degreesOfFreedom)
                        / (degreesOfFreedom * observations))) : 0.0;
            double srmr = srmr(sampleCovariance, fitted.covariance(), variables);
            int free = model.freeParameterCount();
            double logLikelihood = -fitted.negativeLogLikelihood();
            return new SemFitResult(estimates, fitted.covariance(), logLikelihood,
                chiSquare, degreesOfFreedom, pValue, cfi, tli, rmsea, srmr,
                -2.0 * logLikelihood + 2.0 * free,
                -2.0 * logLikelihood + Math.log(observations) * free,
                observations, optimum.evaluations(), optimum.converged(),
                context.provenance());
        }
    }

    private static Evaluation evaluate(
            double[] parameters, SemModel model, double[] sample,
            int observations, ComputeBackend backend, boolean withGradient) {
        int size = model.variables().size();
        double[] paths = new double[size * size];
        double[] residual = new double[size * size];
        for (SemModel.Element element : model.elements()) {
            double value = element.fixed() ? element.start()
                : parameters[model.freeIndex(element.label())];
            if (!element.fixed() && element.kind() == SemModel.Kind.VARIANCE)
                value = Math.exp(value);
            if (element.kind() == SemModel.Kind.REGRESSION) {
                paths[element.first() * size + element.second()] = value;
            } else if (element.kind() == SemModel.Kind.VARIANCE) {
                residual[element.first() * size + element.first()] = value;
            } else {
                residual[element.first() * size + element.second()] = value;
                residual[element.second() * size + element.first()] = value;
            }
        }
        double[] identityMinusPaths = MatrixOps.identity(size);
        for (int index = 0; index < paths.length; index++) identityMinusPaths[index] -= paths[index];
        double[] inverse;
        try {
            inverse = inverse(identityMinusPaths, size);
        } catch (IllegalArgumentException exception) {
            return Evaluation.invalid(size);
        }
        double[] temporary = MatrixOps.multiply(
            backend, inverse, size, size, residual, size);
        double[] implied = new double[size * size];
        backend.dgemm(MatrixTranspose.NONE, MatrixTranspose.TRANSPOSE,
            size, size, size, 1.0, temporary, inverse, 0.0, implied);
        try {
            CholeskyFactor factor = backend.dpotrf(implied, size);
            double[] solved = factor.solve(sample, size);
            double trace = 0.0;
            for (int index = 0; index < size; index++) trace += solved[index * size + index];
            double logDet = factor.logDeterminant();
            double nll = 0.5 * observations
                * (size * Math.log(2.0 * Math.PI) + logDet + trace);
            if (!withGradient)
                return new Evaluation(implied, inverse, nll, logDet, trace, null);
            double[] precision = factor.solve(MatrixOps.identity(size), size);
            double[] precisionSamplePrecision = MatrixOps.multiply(
                backend, solved, size, size, precision, size);
            double[] score = new double[size * size];
            for (int index = 0; index < score.length; index++)
                score[index] = 0.5 * (precision[index] - precisionSamplePrecision[index]);
            double[] scoreTimesInverse = MatrixOps.multiply(
                backend, score, size, size, inverse, size);
            double[] residualScore = MatrixOps.transposeMultiply(
                backend, inverse, size, size, scoreTimesInverse, size);
            double[] covarianceScore = MatrixOps.multiply(
                backend, implied, size, size, scoreTimesInverse, size);
            double[] gradient = new double[model.freeParameterCount()];
            for (SemModel.Element element : model.elements()) {
                if (element.fixed()) continue;
                int parameter = model.freeIndex(element.label());
                if (element.kind() == SemModel.Kind.REGRESSION) {
                    gradient[parameter] += 2.0 * covarianceScore[
                        element.second() * size + element.first()];
                } else if (element.kind() == SemModel.Kind.VARIANCE) {
                    gradient[parameter] += residualScore[
                        element.first() * size + element.first()]
                        * Math.exp(parameters[parameter]);
                } else {
                    gradient[parameter] += 2.0 * residualScore[
                        element.first() * size + element.second()];
                }
            }
            return new Evaluation(implied, inverse, nll, logDet, trace, gradient);
        } catch (IllegalArgumentException | IllegalStateException exception) {
            return Evaluation.invalid(size);
        }
    }

    private static double[] initial(SemModel model) {
        double[] result = new double[model.freeParameterCount()];
        boolean[] set = new boolean[result.length];
        for (SemModel.Element element : model.elements()) {
            if (element.fixed()) continue;
            int index = model.freeIndex(element.label());
            if (!set[index]) {
                result[index] = element.kind() == SemModel.Kind.VARIANCE
                    ? Math.log(element.start()) : element.start();
                set[index] = true;
            }
        }
        return result;
    }

    private static List<SemParameterEstimate> parameterEstimates(
            double[] raw, double[] covariance, SemModel model) {
        List<SemParameterEstimate> result = new ArrayList<>(raw.length);
        for (int index = 0; index < raw.length; index++) {
            double estimate = model.varianceParameter(index) ? Math.exp(raw[index]) : raw[index];
            double derivative = model.varianceParameter(index) ? estimate : 1.0;
            double se = derivative * Math.sqrt(Math.max(0.0,
                covariance[index * raw.length + index]));
            double z = se == 0.0 ? Math.copySign(Double.POSITIVE_INFINITY, estimate) : estimate / se;
            double p = Math.min(1.0, 2.0 * Normal.cumulative(
                Math.abs(z), 0.0, 1.0, false, false));
            result.add(new SemParameterEstimate(
                model.freeParameterLabels().get(index), estimate, se, z, p));
        }
        return List.copyOf(result);
    }

    private static Optimum optimize(double[] initial, double[] lower, double[] upper,
            Objective objective, int maximum, double tolerance, int objectiveScale) {
        Evaluation current = objective.evaluate(initial);
        if (initial.length == 0)
            return new Optimum(initial, current.negativeLogLikelihood(), 1, true);
        double[] point = initial.clone();
        int size = point.length;
        double value = current.negativeLogLikelihood() / objectiveScale;
        double[] gradient = current.gradient();
        double[] inverseHessian = MatrixOps.identity(size);
        int evaluations = 1;
        boolean converged = gradientNorm(gradient)
            <= tolerance * (1.0 + Math.abs(value));
        while (evaluations < maximum && !converged) {
            double[] direction = new double[size];
            for (int row = 0; row < size; row++)
                for (int column = 0; column < size; column++)
                    direction[row] -= inverseHessian[row * size + column] * gradient[column];
            double directionalDerivative = dot(gradient, direction);
            if (!(directionalDerivative < 0.0)) {
                for (int index = 0; index < size; index++) direction[index] = -gradient[index];
                directionalDerivative = -dot(gradient, gradient);
                inverseHessian = MatrixOps.identity(size);
            }
            double step = 1.0;
            Evaluation candidate = null;
            double[] candidatePoint = new double[size];
            while (evaluations < maximum && step >= 1e-12) {
                boolean moved = false;
                for (int index = 0; index < size; index++) {
                    candidatePoint[index] = Math.max(lower[index],
                        Math.min(upper[index], point[index] + step * direction[index]));
                    moved |= candidatePoint[index] != point[index];
                }
                if (!moved) break;
                candidate = objective.evaluate(candidatePoint);
                evaluations++;
                double candidateValue = candidate.negativeLogLikelihood() / objectiveScale;
                if (Double.isFinite(candidateValue)
                        && candidateValue
                            <= value + 1e-4 * step * directionalDerivative) break;
                candidate = null;
                step *= 0.5;
            }
            if (candidate == null) break;
            double[] nextGradient = candidate.gradient();
            double[] displacement = new double[size];
            double[] gradientChange = new double[size];
            for (int index = 0; index < size; index++) {
                displacement[index] = candidatePoint[index] - point[index];
                gradientChange[index] = nextGradient[index] - gradient[index];
            }
            double curvature = dot(displacement, gradientChange);
            if (curvature > 1e-12 * Math.sqrt(
                    dot(displacement, displacement) * dot(gradientChange, gradientChange))) {
                double[] hessianTimesChange = new double[size];
                for (int row = 0; row < size; row++)
                    for (int column = 0; column < size; column++)
                        hessianTimesChange[row] += inverseHessian[row * size + column]
                            * gradientChange[column];
                double changeQuadratic = dot(gradientChange, hessianTimesChange);
                double coefficient = (curvature + changeQuadratic)
                    / (curvature * curvature);
                for (int row = 0; row < size; row++) {
                    for (int column = 0; column < size; column++) {
                        inverseHessian[row * size + column] += coefficient
                            * displacement[row] * displacement[column]
                            - (hessianTimesChange[row] * displacement[column]
                                + displacement[row] * hessianTimesChange[column]) / curvature;
                    }
                }
            } else inverseHessian = MatrixOps.identity(size);
            double relativeStep = 0.0;
            for (int index = 0; index < size; index++)
                relativeStep = Math.max(relativeStep, Math.abs(displacement[index])
                    / (1.0 + Math.abs(candidatePoint[index])));
            point = candidatePoint.clone();
            value = candidate.negativeLogLikelihood() / objectiveScale;
            gradient = nextGradient;
            converged = gradientNorm(gradient)
                    <= tolerance * (1.0 + Math.abs(value))
                || relativeStep <= tolerance;
        }
        return new Optimum(point, value * objectiveScale, evaluations, converged);
    }

    private static double[] parameterCovariance(double[] point, SemModel model,
            Evaluation fitted, int observations, ComputeBackend backend) {
        int parameters = point.length;
        if (parameters == 0) return new double[0];
        int variables = model.variables().size();
        double[][] derivatives = new double[parameters][variables * variables];
        for (SemModel.Element element : model.elements()) {
            if (element.fixed()) continue;
            int parameter = model.freeIndex(element.label());
            double[] derivative = derivatives[parameter];
            for (int row = 0; row < variables; row++) {
                for (int column = 0; column < variables; column++) {
                    double value;
                    if (element.kind() == SemModel.Kind.REGRESSION) {
                        value = fitted.inversePath()[row * variables + element.first()]
                            * fitted.covariance()[element.second() * variables + column]
                            + fitted.covariance()[row * variables + element.second()]
                            * fitted.inversePath()[column * variables + element.first()];
                    } else if (element.kind() == SemModel.Kind.VARIANCE) {
                        value = Math.exp(point[parameter])
                            * fitted.inversePath()[row * variables + element.first()]
                            * fitted.inversePath()[column * variables + element.first()];
                    } else {
                        value = fitted.inversePath()[row * variables + element.first()]
                            * fitted.inversePath()[column * variables + element.second()]
                            + fitted.inversePath()[row * variables + element.second()]
                            * fitted.inversePath()[column * variables + element.first()];
                    }
                    derivative[row * variables + column] += value;
                }
            }
        }
        CholeskyFactor factor = backend.dpotrf(fitted.covariance(), variables);
        double[][] precisionDerivatives = new double[parameters][];
        for (int parameter = 0; parameter < parameters; parameter++)
            precisionDerivatives[parameter] = factor.solve(derivatives[parameter], variables);
        double[] information = new double[parameters * parameters];
        for (int first = 0; first < parameters; first++) {
            for (int second = 0; second <= first; second++) {
                double trace = 0.0;
                for (int row = 0; row < variables; row++)
                    for (int column = 0; column < variables; column++)
                        trace += precisionDerivatives[first][row * variables + column]
                            * precisionDerivatives[second][column * variables + row];
                double informationValue = 0.5 * observations * trace;
                information[first * parameters + second] = informationValue;
                information[second * parameters + first] = informationValue;
            }
        }
        double ridge = 1e-8;
        for (int attempt = 0; attempt < 12; attempt++) {
            double[] value = information.clone();
            for (int index = 0; index < parameters; index++)
                value[index * parameters + index] += ridge;
            try {
                return backend.dpotrf(value, parameters)
                    .solve(MatrixOps.identity(parameters), parameters);
            } catch (IllegalArgumentException | IllegalStateException ignored) { ridge *= 10.0; }
        }
        double[] result = new double[parameters * parameters];
        java.util.Arrays.fill(result, Double.NaN);
        return result;
    }

    private static double dot(double[] first, double[] second) {
        double result = 0.0;
        for (int index = 0; index < first.length; index++)
            result += first[index] * second[index];
        return result;
    }

    private static double gradientNorm(double[] gradient) {
        double result = 0.0;
        for (double value : gradient) result = Math.max(result, Math.abs(value));
        return result;
    }

    private static double baselineDiscrepancy(double[] sample, int size, double logDet) {
        double sumLog = 0.0;
        for (int index = 0; index < size; index++) sumLog += Math.log(sample[index * size + index]);
        return sumLog - logDet;
    }

    private static double srmr(double[] sample, double[] implied, int size) {
        double sum = 0.0;
        int count = 0;
        for (int first = 0; first < size; first++) {
            for (int second = 0; second <= first; second++) {
                double observed = sample[first * size + second]
                    / Math.sqrt(sample[first * size + first] * sample[second * size + second]);
                double fitted = implied[first * size + second]
                    / Math.sqrt(sample[first * size + first] * sample[second * size + second]);
                double difference = observed - fitted;
                sum += difference * difference;
                count++;
            }
        }
        return Math.sqrt(sum / count);
    }

    private static double[] inverse(double[] matrix, int size) {
        double[] work = matrix.clone();
        double[] result = MatrixOps.identity(size);
        for (int column = 0; column < size; column++) {
            int pivot = column;
            for (int row = column + 1; row < size; row++)
                if (Math.abs(work[row * size + column]) > Math.abs(work[pivot * size + column])) pivot = row;
            if (Math.abs(work[pivot * size + column]) < 1e-12)
                throw new IllegalArgumentException("singular path system");
            if (pivot != column) {
                for (int index = 0; index < size; index++) {
                    double temporary = work[column * size + index];
                    work[column * size + index] = work[pivot * size + index];
                    work[pivot * size + index] = temporary;
                    temporary = result[column * size + index];
                    result[column * size + index] = result[pivot * size + index];
                    result[pivot * size + index] = temporary;
                }
            }
            double diagonal = work[column * size + column];
            for (int index = 0; index < size; index++) {
                work[column * size + index] /= diagonal;
                result[column * size + index] /= diagonal;
            }
            for (int row = 0; row < size; row++) {
                if (row == column) continue;
                double scale = work[row * size + column];
                for (int index = 0; index < size; index++) {
                    work[row * size + index] -= scale * work[column * size + index];
                    result[row * size + index] -= scale * result[column * size + index];
                }
            }
        }
        return result;
    }

    @FunctionalInterface
    private interface Objective { Evaluation evaluate(double[] point); }

    private record Evaluation(double[] covariance, double[] inversePath,
                              double negativeLogLikelihood,
                              double logDeterminant, double trace,
                              double[] gradient) {
        static Evaluation invalid(int size) {
            return new Evaluation(new double[size * size], new double[size * size],
                Double.MAX_VALUE / 8.0, Double.NaN, Double.NaN, null);
        }
    }
    private record Optimum(double[] parameters, double objective,
                           int evaluations, boolean converged) { }
}
