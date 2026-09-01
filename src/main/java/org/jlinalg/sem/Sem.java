/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import jdistlib.ChiSquare;
import jdistlib.Normal;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTranspose;
import jdistlib.math.MultivariableFunction;
import jdistlib.math.opt.Bobyqa;
import jdistlib.math.opt.OptimizationResult;
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
        List<double[]> complete = new ArrayList<>();
        for (int row = 0; row < data.length; row++) {
            if (data[row] == null || data[row].length != variables)
                throw new IllegalArgumentException("SEM data must be rectangular and match variables");
            boolean finite = true;
            for (double value : data[row]) finite &= Double.isFinite(value);
            if (finite) complete.add(data[row]);
            else if (options.missingDataPolicy() == MissingDataPolicy.ERROR)
                throw new IllegalArgumentException("non-finite SEM row: " + row);
        }
        if (complete.size() <= variables)
            throw new IllegalArgumentException("too few complete rows for SEM covariance");
        double[] means = new double[variables];
        for (double[] row : complete)
            for (int variable = 0; variable < variables; variable++) means[variable] += row[variable];
        for (int variable = 0; variable < variables; variable++) means[variable] /= complete.size();
        double[] covariance = new double[variables * variables];
        for (double[] row : complete) {
            for (int first = 0; first < variables; first++) {
                for (int second = 0; second < variables; second++) {
                    covariance[first * variables + second] +=
                        (row[first] - means[first]) * (row[second] - means[second]);
                }
            }
        }
        for (int index = 0; index < covariance.length; index++) covariance[index] /= complete.size();
        return fitCovariance(covariance, complete.size(), model, options, backendPolicy);
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
            MultivariableFunction objective = point -> evaluate(
                point, model, sampleCovariance, observations, backend).negativeLogLikelihood();
            Optimum optimum = optimize(initial, lower, upper, objective,
                options.maximumEvaluations(), options.tolerance());
            Evaluation fitted = evaluate(optimum.parameters(), model,
                sampleCovariance, observations, backend);
            double[] rawCovariance = inverseHessian(
                optimum.parameters(), objective, backend);
            List<SemParameterEstimate> estimates = parameterEstimates(
                optimum.parameters(), rawCovariance, model);
            double logDetSample = sampleFactor.logDeterminant();
            double discrepancy = fitted.logDeterminant() + fitted.trace()
                - logDetSample - variables;
            double chiSquare = Math.max(0.0, (observations - 1.0) * discrepancy);
            double pValue = degreesOfFreedom > 0
                ? ChiSquare.cumulative(chiSquare, degreesOfFreedom, false, false)
                : Double.NaN;
            double baselineDiscrepancy = baselineDiscrepancy(
                sampleCovariance, variables, logDetSample);
            double baselineChi = (observations - 1.0) * baselineDiscrepancy;
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
                        / (degreesOfFreedom * (observations - 1.0)))) : 0.0;
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
            int observations, ComputeBackend backend) {
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
            return new Evaluation(implied, nll, logDet, trace);
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
            MultivariableFunction objective, int maximum, double tolerance) {
        if (initial.length == 0) return new Optimum(initial, objective.eval(initial), 1, true);
        try {
            OptimizationResult result = Bobyqa.bobyqa(initial, lower, upper,
                objective, 2 * initial.length + 1, 0.5, tolerance, maximum, true);
            if (result.mX != null && Double.isFinite(result.mF))
                return new Optimum(result.mX, result.mF, result.numFunctionCalls, true);
        } catch (RuntimeException ignored) { }
        double[] point = initial.clone();
        double value = objective.eval(point);
        int evaluations = 1;
        double step = 0.5;
        while (evaluations < maximum && step > tolerance) {
            boolean improved = false;
            for (int dimension = 0; dimension < point.length && evaluations < maximum; dimension++) {
                double original = point[dimension];
                for (int direction : new int[] {-1, 1}) {
                    point[dimension] = Math.max(lower[dimension],
                        Math.min(upper[dimension], original + direction * step));
                    double candidate = objective.eval(point);
                    evaluations++;
                    if (candidate < value) {
                        value = candidate; original = point[dimension]; improved = true;
                    } else point[dimension] = original;
                }
            }
            if (!improved) step *= 0.5;
        }
        return new Optimum(point, value, evaluations, step <= Math.sqrt(tolerance));
    }

    private static double[] inverseHessian(
            double[] point, MultivariableFunction objective, ComputeBackend backend) {
        int size = point.length;
        if (size == 0) return new double[0];
        double[] hessian = new double[size * size];
        double center = objective.eval(point);
        for (int first = 0; first < size; first++) {
            double h1 = 1e-4 * (1.0 + Math.abs(point[first]));
            double[] plus = point.clone(); plus[first] += h1;
            double[] minus = point.clone(); minus[first] -= h1;
            hessian[first * size + first] =
                (objective.eval(plus) - 2.0 * center + objective.eval(minus)) / (h1 * h1);
            for (int second = 0; second < first; second++) {
                double h2 = 1e-4 * (1.0 + Math.abs(point[second]));
                double[] pp = point.clone(); pp[first] += h1; pp[second] += h2;
                double[] pm = point.clone(); pm[first] += h1; pm[second] -= h2;
                double[] mp = point.clone(); mp[first] -= h1; mp[second] += h2;
                double[] mm = point.clone(); mm[first] -= h1; mm[second] -= h2;
                double value = (objective.eval(pp) - objective.eval(pm)
                    - objective.eval(mp) + objective.eval(mm)) / (4.0 * h1 * h2);
                hessian[first * size + second] = value;
                hessian[second * size + first] = value;
            }
        }
        double ridge = 1e-8;
        for (int attempt = 0; attempt < 12; attempt++) {
            double[] value = hessian.clone();
            for (int index = 0; index < size; index++) value[index * size + index] += ridge;
            try {
                return backend.dpotrf(value, size).solve(MatrixOps.identity(size), size);
            } catch (IllegalArgumentException | IllegalStateException ignored) { ridge *= 10.0; }
        }
        double[] result = new double[size * size];
        Arrays.fill(result, Double.NaN);
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
                    / Math.sqrt(implied[first * size + first] * implied[second * size + second]);
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

    private record Evaluation(double[] covariance, double negativeLogLikelihood,
                              double logDeterminant, double trace) {
        static Evaluation invalid(int size) {
            return new Evaluation(new double[size * size],
                Double.MAX_VALUE / 8.0, Double.NaN, Double.NaN);
        }
    }
    private record Optimum(double[] parameters, double objective,
                           int evaluations, boolean converged) { }
}
