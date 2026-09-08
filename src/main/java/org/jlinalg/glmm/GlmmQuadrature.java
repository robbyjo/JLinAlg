/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import jdistlib.math.MathFunctions;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glm.GlmFamily;

/** Maximum marginal likelihood for independent Gaussian random intercepts.
 * Supports built-in binomial/logit and Poisson/log families, dispersion one.
 * Integrates around each group's posterior mode and observed curvature and
 * checks successive Hermite rules. This is numerical, not exact, integration;
 * correlated pedigree effects and random slopes are not supported. */
public final class GlmmQuadrature {
    private GlmmQuadrature() { }

    public static GlmmQuadratureResult fit(double[] y, double[][] x,
            List<String> groups, GlmFamily family) {
        return fit(y, x, groups, family, GlmmQuadratureOptions.defaults());
    }
    public static GlmmQuadratureResult fit(double[] y, double[][] x,
            List<String> groups, GlmFamily family, GlmmQuadratureOptions options) {
        return fit(y, x, groups, family, null, null, options);
    }

    /** Binomial responses are proportions; trials are positive integer counts.
     * For Poisson, trials must be null or all one. Null offsets mean zero. */
    public static GlmmQuadratureResult fit(double[] y, double[][] x,
            List<String> groups, GlmFamily family, double[] trials, double[] offsets,
            GlmmQuadratureOptions options) {
        Data data = new Data(y, x, groups, family, trials, offsets, options);
        if (data.rows.size() < 2) throw new IllegalArgumentException("fitting needs at least two independent groups");
        if (data.p >= y.length) throw new IllegalArgumentException("fixed design must have fewer columns than observations");
        data.checkRank();
        QuadratureOptimizer.Objective objective = new QuadratureOptimizer.Objective() {
            @Override public double value(double[] parameters) {
                data.evaluations++;
                GlmmQuadratureEvaluation value = data.evaluate(data.physical(parameters));
                return value.converged() ? -value.logLikelihood() : Double.POSITIVE_INFINITY;
            }
            @Override public double[] exactGradient(double[] parameters) {
                return parameters[data.p] == 0 ? data.boundaryGradient(parameters) : null;
            }
        };
        // Fit the actual zero-variance model as well as two interior starts.
        QuadratureOptimizer.Fit boundary = QuadratureOptimizer.minimize(objective,
            new double[data.p + 1], data.p, true, options);
        QuadratureOptimizer.Fit best = boundary;
        for (double coordinate : new double[] {Math.log(2), Math.log1p(2 / data.varianceScale)}) {
            double[] start = boundary.parameters().clone(); start[data.p] = coordinate;
            QuadratureOptimizer.Fit candidate = QuadratureOptimizer.minimize(objective,
                start, data.p, false, options);
            if (candidate.value() < best.value()) best = candidate;
        }
        double[] point = best.parameters().clone(), parameters = data.physical(point);
        double[] fullGradient = QuadratureOptimizer.gradient(objective, point, data.p, false);
        double gradientNorm = QuadratureOptimizer.projectedNorm(fullGradient, point, data.p, false);
        GlmmQuadratureEvaluation integral = data.evaluate(parameters);
        boolean boundaryVariance = parameters[data.p] == 0;
        boolean stationary = Double.isFinite(gradientNorm) && gradientNorm <= options.gradientTolerance();
        double[][] covariance = QuadratureOptimizer.covariance(objective, point,
            boundaryVariance ? data.p : data.p + 1);
        boolean informationValid = covariance != null;
        // A vanishing score at infinity (e.g. complete separation) is not a
        // finite optimum: its Newton correction remains appreciable.
        if (informationValid) for (int i = 0; i < covariance.length; i++) {
            double correction = 0;
            for (int j = 0; j < covariance.length; j++) correction += covariance[i][j] * fullGradient[j];
            stationary &= Math.abs(correction) < 1e-3 * (1 + Math.abs(point[i]));
        }
        boolean converged = stationary && integral.converged() && informationValid;
        double[][] joint = new double[data.p + 1][data.p + 1];
        for (double[] row : joint) Arrays.fill(row, Double.NaN);
        if (informationValid) for (int i = 0; i < covariance.length; i++)
            for (int j = 0; j < covariance.length; j++)
                joint[i][j] = covariance[i][j]
                    * (i < data.p ? data.coefficientMultipliers[i] / data.scales[i] : parameters[data.p] + data.varianceScale)
                    * (j < data.p ? data.coefficientMultipliers[j] / data.scales[j] : parameters[data.p] + data.varianceScale);
        if (!converged) for (double[] row : joint) Arrays.fill(row, Double.NaN);
        double[][] betaCovariance = new double[data.p][data.p];
        for (int i = 0; i < data.p; i++) System.arraycopy(joint[i], 0, betaCovariance[i], 0, data.p);
        double[] beta = Arrays.copyOf(parameters, data.p);
        for (int i = 0; i < beta.length; i++) beta[i] /= data.scales[i];
        String status = !integral.converged() ? "quadrature tolerance not met"
            : !stationary ? "marginal likelihood score tolerance not met"
            : !informationValid ? "observed information is not positive definite"
            : boundaryVariance ? "converged at zero variance; fixed-effect covariance conditional on variance zero"
            : "converged";
        return new GlmmQuadratureResult(family.name(), beta, Math.sqrt(parameters[data.p]),
            data.marginalMeans(parameters), integral.logLikelihood(), integral.nodes(), data.evaluations, converged,
            betaCovariance, joint, gradientNorm, integral.estimatedError(), integral.converged(),
            converged && !boundaryVariance, boundaryVariance, status);
    }

    /** Evaluate the ordered-observation likelihood at supplied parameters. */
    public static GlmmQuadratureEvaluation evaluate(double[] y, double[][] x,
            List<String> groups, GlmFamily family, double[] beta, double sd,
            GlmmQuadratureOptions options) {
        return evaluate(y, x, groups, family, null, null, beta, sd, options);
    }
    public static GlmmQuadratureEvaluation evaluate(double[] y, double[][] x,
            List<String> groups, GlmFamily family, double[] trials, double[] offsets,
            double[] beta, double sd, GlmmQuadratureOptions options) {
        Data data = new Data(y, x, groups, family, trials, offsets, options);
        if (beta == null || beta.length != data.p || !Double.isFinite(sd) || sd < 0
                || !Double.isFinite(sd * sd)) throw new IllegalArgumentException("invalid quadrature parameters");
        double[] parameters = new double[data.p + 1];
        for (int j = 0; j < data.p; j++) {
            if (!Double.isFinite(beta[j])) throw new IllegalArgumentException("nonfinite coefficient");
            parameters[j] = beta[j] * data.scales[j];
        }
        parameters[data.p] = sd * sd;
        return data.evaluate(parameters);
    }

    private static final class Data {
        final int p;
        final boolean binomial;
        final double[][] x;
        final double[] successes, trials, offsets, constants, scales;
        final double[] coefficientMultipliers;
        final double varianceScale;
        final List<int[]> rows = new ArrayList<>();
        final GlmmQuadratureOptions options;
        int evaluations;

        Data(double[] y, double[][] design, List<String> groups, GlmFamily family,
                double[] trialCounts, double[] offset, GlmmQuadratureOptions controls) {
            if (family != GlmFamilies.binomial() && family != GlmFamilies.poisson())
                throw new IllegalArgumentException("quadrature supports only built-in binomial(logit) and Poisson(log), dispersion one");
            if (y == null || design == null || groups == null || controls == null || y.length == 0
                    || y.length != design.length || y.length != groups.size() || design[0] == null
                    || design[0].length == 0 || (trialCounts != null && trialCounts.length != y.length)
                    || (offset != null && offset.length != y.length))
                throw new IllegalArgumentException("invalid quadrature data dimensions");
            options = controls; binomial = family == GlmFamilies.binomial(); p = design[0].length;
            x = new double[y.length][p]; successes = new double[y.length]; trials = new double[y.length];
            offsets = new double[y.length]; constants = new double[y.length]; scales = new double[p];
            LinkedHashMap<String, List<Integer>> map = new LinkedHashMap<>();
            for (int i = 0; i < y.length; i++) {
                if (design[i] == null || design[i].length != p || groups.get(i) == null || groups.get(i).isBlank())
                    throw new IllegalArgumentException("invalid quadrature row");
                trials[i] = trialCounts == null ? 1 : trialCounts[i];
                offsets[i] = offset == null ? 0 : offset[i];
                if (!Double.isFinite(y[i]) || !Double.isFinite(trials[i]) || !Double.isFinite(offsets[i])
                        || trials[i] < 1 || trials[i] != Math.rint(trials[i])
                        || (!binomial && trials[i] != 1)) throw new IllegalArgumentException("invalid trials, offset or response");
                double s = binomial ? y[i] * trials[i] : y[i];
                // A count divided by its trial count can multiply back one ULP
                // away from the integer. An absolute tolerance rejects valid
                // large denominators; genuinely fractional counts still fail.
                double integerTolerance = binomial && trials[i] > 1 ? 2 * Math.ulp(s) : 0;
                if (!Double.isFinite(s) || s < 0 || Math.abs(s - Math.rint(s)) > integerTolerance || (binomial && s > trials[i]))
                    throw new IllegalArgumentException("responses must encode integer counts (Bernoulli without trials)");
                successes[i] = Math.rint(s);
                constants[i] = -MathFunctions.lgammafn(successes[i] + 1);
                if (binomial) constants[i] += MathFunctions.lgammafn(trials[i] + 1)
                    - MathFunctions.lgammafn(trials[i] - successes[i] + 1);
                if (binomial && trials[i] > 64 && successes[i] > 0 && successes[i] < trials[i]) {
                    double failures = trials[i] - successes[i];
                    constants[i] = MathFunctions.stirlerr(trials[i]) - MathFunctions.stirlerr(successes[i])
                        - MathFunctions.stirlerr(failures)
                        - .5 * (Math.log(2 * Math.PI) + Math.log(successes[i]) + Math.log(failures / trials[i]));
                } else if (!binomial && successes[i] > 64)
                    constants[i] = -MathFunctions.stirlerr(successes[i]) - .5 * (Math.log(2 * Math.PI) + Math.log(successes[i]));
                for (int j = 0; j < p; j++) {
                    if (!Double.isFinite(design[i][j])) throw new IllegalArgumentException("nonfinite fixed design");
                    x[i][j] = design[i][j]; scales[j] = Math.max(scales[j], Math.abs(x[i][j]));
                }
                map.computeIfAbsent(groups.get(i), ignored -> new ArrayList<>()).add(i);
            }
            for (int j = 0; j < p; j++) {
                if (scales[j] == 0) scales[j] = 1;
                for (double[] row : x) row[j] /= scales[j];
            }
            for (List<Integer> group : map.values()) rows.add(group.stream().mapToInt(Integer::intValue).toArray());
            coefficientMultipliers = new double[p];
            double maxGroupInformation = 1;
            for (int[] group : rows) {
                double information = 0;
                for (int row : group) {
                    double weight = binomial ? trials[row] / 4 : Math.max(1, successes[row]);
                    information += weight;
                    for (int j = 0; j < p; j++) coefficientMultipliers[j] += weight * x[row][j] * x[row][j];
                }
                maxGroupInformation = Math.max(maxGroupInformation, information);
            }
            varianceScale = 1 / maxGroupInformation;
            for (int j = 0; j < p; j++) coefficientMultipliers[j] = 1 / Math.sqrt(Math.max(1, coefficientMultipliers[j]));
        }
        // Optimize t=log1p(variance/v0), with v0 a within-group information
        // scale. This resolves tiny concentrated optima without imposing a
        // positive lower bound or using an absolute variance derivative step.
        double[] physical(double[] point) {
            double[] result = point.clone();
            for (int j = 0; j < p; j++) result[j] *= coefficientMultipliers[j];
            result[p] = varianceScale * Math.expm1(point[p]);
            return result;
        }
        double[] boundaryGradient(double[] point) {
            double[] eta = eta(physical(point)), gradient = new double[p + 1];
            for (int[] group : rows) {
                double score = 0, curvature = 0;
                for (int row : group) {
                    double residual, weight;
                    if (binomial) {
                        double probability = logistic(eta[row]), failure = logistic(-eta[row]);
                        residual = successes[row] * failure - (trials[row] - successes[row]) * probability;
                        weight = trials[row] * probability * failure;
                    } else {
                        double mean = Math.exp(eta[row]); residual = successes[row] - mean; weight = mean;
                    }
                    score += residual; curvature += weight;
                    for (int j = 0; j < p; j++) gradient[j] -= x[row][j] * residual * coefficientMultipliers[j];
                }
                // Gaussian heat equation at variance zero: d log L/dv=(l''+l'^2)/2.
                gradient[p] += .5 * (curvature - score * score) * varianceScale;
            }
            return gradient;
        }
        void checkRank() {
            double[][] gram = new double[p][p];
            for (double[] row : x) for (int i = 0; i < p; i++) for (int j = 0; j < p; j++)
                gram[i][j] += row[i] * row[j];
            if (QuadratureOptimizer.inversePositiveDefinite(gram) == null)
                throw new IllegalArgumentException("fixed design is rank deficient");
        }
        double[] eta(double[] parameters) {
            double[] eta = offsets.clone();
            for (int i = 0; i < x.length; i++) for (int j = 0; j < p; j++) eta[i] += x[i][j] * parameters[j];
            return eta;
        }
        GlmmQuadratureEvaluation evaluate(double[] parameters) {
            double variance = parameters[p];
            if (!Double.isFinite(variance) || variance < 0) return failed();
            double sd = Math.sqrt(variance), total = 0, error = 0;
            double[] eta = eta(parameters);
            for (double value : eta) if (!Double.isFinite(value)) return failed();
            int used = 0; boolean converged = true;
            for (int[] group : rows) {
                GlmmQuadratureEvaluation value = integrateGroup(group, eta, sd, options.quadratureTolerance() / rows.size());
                total += value.logLikelihood(); error += value.estimatedError(); used = Math.max(used, value.nodes());
                converged &= value.converged();
            }
            return new GlmmQuadratureEvaluation(total, used, error, converged && Double.isFinite(total));
        }
        private GlmmQuadratureEvaluation integrateGroup(int[] group, double[] eta, double sd, double tolerance) {
            if (sd == 0) {
                double value = logKernel(group, eta, 0, 0);
                return new GlmmQuadratureEvaluation(value, 0, 0, Double.isFinite(value));
            }
            double[] mode = mode(group, eta, sd);
            if (mode == null) return failed();
            int n = options.initialNodes(), stable = 0;
            double previous = integrate(group, eta, sd, mode, n), lastError = Double.POSITIVE_INFINITY;
            while (n < options.maximumNodes()) {
                n = Math.min(options.maximumNodes(), n + Math.max(4, n / 2));
                double value = integrate(group, eta, sd, mode, n);
                lastError = Math.abs(value - previous);
                stable = lastError <= tolerance ? stable + 1 : 0;
                previous = value;
                if (stable >= 2) break;
            }
            return new GlmmQuadratureEvaluation(previous, n, lastError, stable >= 2 && Double.isFinite(previous));
        }
        private GlmmQuadratureEvaluation failed() {
            return new GlmmQuadratureEvaluation(Double.NEGATIVE_INFINITY, 0, Double.POSITIVE_INFINITY, false);
        }
        // Standard-normal coordinates u=b/sd avoid degeneracy near variance zero.
        // Strict concavity guarantees a unique mode for both supported families.
        private double[] mode(int[] group, double[] eta, double sd) {
            double lo = -1, hi = 1;
            for (int i = 0; i < 60 && derivatives(group, eta, sd, lo)[0] < 0; i++) lo *= 2;
            for (int i = 0; i < 60 && derivatives(group, eta, sd, hi)[0] > 0; i++) hi *= 2;
            if (!(derivatives(group, eta, sd, lo)[0] >= 0) || !(derivatives(group, eta, sd, hi)[0] <= 0)) return null;
            double u = 0;
            for (int i = 0; i < 150; i++) {
                double[] k = derivatives(group, eta, sd, u);
                if (Double.isFinite(k[0]) && Double.isFinite(k[1])
                        && Math.abs(k[0]) / k[1] <= 2e-12 * (1 + Math.abs(u)))
                    return new double[] {u, 1 / Math.sqrt(k[1])};
                if (k[0] > 0) lo = u; else hi = u;
                double next = u + k[0] / k[1];
                // Newton moves only about one unit per iteration in the far
                // Poisson tail. Force bracket contraction instead of spending
                // the iteration budget hundreds of units away from the mode.
                double margin = .2 * (hi - lo);
                if (!Double.isFinite(next) || next <= lo + margin || next >= hi - margin) next = (lo + hi) / 2;
                u = next;
            }
            return null;
        }
        private double[] derivatives(int[] group, double[] eta, double sd, double u) {
            double score = -u, curvature = 1;
            for (int row : group) {
                double predictor = eta[row] + sd * u, s = successes[row];
                if (binomial) {
                    double probability = logistic(predictor), failure = logistic(-predictor);
                    score += sd * (s * failure - (trials[row] - s) * probability);
                    curvature += sd * sd * trials[row] * probability * failure;
                } else {
                    double mean = Math.exp(predictor);
                    if (sd != 0) { score += sd * (s - mean); curvature += sd * sd * mean; }
                }
            }
            return new double[] {score, curvature};
        }
        private double logKernel(int[] group, double[] eta, double sd, double u) {
            double value = -.5 * u * u, compensation = 0;
            for (int row : group) {
                double predictor = eta[row] + sd * u, s = successes[row];
                double contribution;
                if (binomial && trials[row] > 64 && s > 0 && s < trials[row]) {
                    double probability = logistic(predictor), failure = logistic(-predictor);
                    double expectedSuccesses = trials[row] * probability, expectedFailures = trials[row] * failure;
                    contribution = constants[row]
                        - deviancePart(s, expectedSuccesses, expectedSuccesses == 0 ? Math.log(trials[row]) - softplus(-predictor) : 0)
                        - deviancePart(trials[row] - s, expectedFailures, expectedFailures == 0 ? Math.log(trials[row]) - softplus(predictor) : 0);
                } else if (binomial) contribution = (trials[row] > 64 ? 0 : constants[row]) + (predictor >= 0
                        ? -(trials[row] - s) * predictor - trials[row] * Math.log1p(Math.exp(-predictor))
                        : s * predictor - trials[row] * Math.log1p(Math.exp(predictor)));
                else {
                    double mean = Math.exp(predictor);
                    contribution = s > 64 ? constants[row] - deviancePart(s, mean, predictor)
                        : constants[row] + (s == 0 ? 0 : s * predictor) - mean;
                }
                if (!Double.isFinite(contribution)) return Double.NEGATIVE_INFINITY;
                double corrected = contribution - compensation, next = value + corrected;
                compensation = (next - value) - corrected; value = next;
            }
            return value;
        }
        private double integrate(int[] group, double[] eta, double sd, double[] mode, int n) {
            AdaptiveHermiteRule rule = AdaptiveHermiteRule.of(n);
            double max = Double.NEGATIVE_INFINITY; double[] terms = new double[n];
            for (int i = 0; i < n; i++) {
                double node = rule.nodes[i];
                terms[i] = rule.logWeights[i] + logKernel(group, eta, sd,
                    mode[0] + Math.sqrt(2) * mode[1] * node) + node * node;
                max = Math.max(max, terms[i]);
            }
            if (!Double.isFinite(max)) return max;
            double sum = 0;
            for (double term : terms) sum += Math.exp(term - max);
            return Math.log(mode[1]) - .5 * Math.log(Math.PI) + max + Math.log(sum);
        }
        double[] marginalMeans(double[] parameters) {
            double[] eta = eta(parameters); double variance = parameters[p];
            if (!binomial) {
                for (int i = 0; i < eta.length; i++) eta[i] = Math.exp(eta[i] + variance / 2);
                return eta;
            }
            // E[p] is a one-success Bernoulli marginal likelihood. Adapt its
            // integrand too, rather than use a low-order rule around the prior.
            double[][] design = new double[eta.length][1];
            double[] one = new double[eta.length]; Arrays.fill(one, 1);
            List<String> labels = new ArrayList<>();
            for (int i = 0; i < eta.length; i++) labels.add(Integer.toString(i));
            Data means = new Data(one, design, labels, GlmFamilies.binomial(), null, eta, options);
            double[] result = new double[eta.length];
            for (int i = 0; i < result.length; i++) {
                int[] group = {i};
                if (variance == 0) { result[i] = logistic(eta[i]); continue; }
                GlmmQuadratureEvaluation mean = means.integrateGroup(group, eta, Math.sqrt(variance), options.quadratureTolerance());
                result[i] = mean.converged() ? Math.exp(mean.logLikelihood()) : Double.NaN;
            }
            return result;
        }
    }
    private static double logistic(double x) {
        return x >= 0 ? 1 / (1 + Math.exp(-x)) : Math.exp(x) / (1 + Math.exp(x));
    }
    /** x*log(x/mu)+mu-x, evaluated without cancellation around x=mu. */
    private static double deviancePart(double x, double mean, double logMean) {
        if (mean == 0) return x * (Math.log(x) - logMean) - x;
        return MathFunctions.bd0(x, mean);
    }
    private static double softplus(double x) { return Math.max(x, 0) + Math.log1p(Math.exp(-Math.abs(x))); }
}
