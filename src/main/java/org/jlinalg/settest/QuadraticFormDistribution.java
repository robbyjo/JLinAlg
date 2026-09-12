/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

import java.util.Arrays;
import jdistlib.ChiSquare;
import jdistlib.Normal;

/** Deterministic positive chi-square-mixture survival probabilities. */
final class QuadraticFormDistribution {
    private static final int MAXIMUM_SERIES_TERMS = 8192;
    private static final double RELATIVE_SERIES_TOLERANCE = 1e-12;

    private QuadraticFormDistribution() { }

    record Tail(double pValue, String method,
        double degreesOfFreedom, double scale) { }

    static Tail survival(double statistic, double[] suppliedEigenvalues) {
        if (!Double.isFinite(statistic) || statistic < 0
                || suppliedEigenvalues == null)
            throw new IllegalArgumentException(
                "quadratic-form statistic and eigenvalues are invalid");
        for (double value : suppliedEigenvalues)
            if (!Double.isFinite(value) || value < 0)
                throw new IllegalArgumentException("kernel eigenvalues must be finite and nonnegative");
        double maximum = Arrays.stream(suppliedEigenvalues).max().orElse(0);
        double[] eigenvalues = Arrays.stream(suppliedEigenvalues)
            .filter(value -> value > 0)
            .map(value -> value / maximum)
            .toArray();
        if (eigenvalues.length == 0)
            throw new IllegalArgumentException(
                "quadratic-form kernel has no positive eigenvalues");
        statistic /= maximum;
        double sum = Arrays.stream(eigenvalues).sum();
        double sumSquares = Arrays.stream(eigenvalues)
            .map(value -> value * value).sum();
        double degrees = sum * sum / sumSquares;
        double scale = sumSquares / sum;
        if (eigenvalues.length == 1)
            return new Tail(probability(ChiSquare.cumulative(
                statistic / eigenvalues[0], 1, false, false)),
                "exact-scaled-chi-square", 1, maximum * eigenvalues[0]);
        boolean equal = true;
        for (int index = 1; index < eigenvalues.length; index++)
            if (Math.abs(eigenvalues[index] - eigenvalues[0])
                    > 1e-12 * Math.max(eigenvalues[index], eigenvalues[0])) {
                equal = false;
                break;
            }
        if (equal)
            return new Tail(probability(ChiSquare.cumulative(
                statistic / eigenvalues[0], eigenvalues.length,
                false, false)), "exact-equal-eigenvalue-chi-square",
                eigenvalues.length, maximum * eigenvalues[0]);
        if (eigenvalues.length == 2) {
            // In two dimensions the Gaussian radius squared is chi-square(2)
            // independently of its uniform angle. This nonoscillatory integral
            // avoids Imhof truncation and moment matching for rank-two kernels.
            return new Tail(probability(angularProbability(statistic, eigenvalues, false)),
                "rank-two-angular-integral", degrees, maximum*scale);
        }
        try {
            return new Tail(probability(gammaSeries(statistic, eigenvalues, false)),
                "positive-gamma-series", degrees, maximum * scale);
        } catch (SeriesResolutionException unresolved) {
            return new Tail(probability(saddlepointSurvival(statistic, eigenvalues)),
                "lugannani-rice-saddlepoint", degrees, maximum * scale);
        }
    }

    static double critical(double[] eigenvalues, double survivalProbability) {
        if (!(survivalProbability > 0 && survivalProbability <= 1)
                || !Double.isFinite(survivalProbability))
            throw new IllegalArgumentException(
                "survival probability must be in (0,1]");
        // Invert the same tail used to compute component p-values. Substituting
        // a moment-matched quantile changes the simulated minimum-p event.
        survival(0, eigenvalues); // validate spectrum
        if (survivalProbability == 1) return 0;
        double maximum = Arrays.stream(eigenvalues).max().orElseThrow();
        double[] normalized = Arrays.stream(eigenvalues).filter(v -> v > 0).map(v -> v / maximum).toArray();
        double minimum = Arrays.stream(normalized).min().orElseThrow();
        double unitQuantile = ChiSquare.quantile(survivalProbability, normalized.length, false, false);
        if (minimum == 1) return maximum * unitQuantile;
        // min(lambda)*chi-square(rank) <= Q <= max(lambda)*chi-square(rank).
        // Positive rank, not array length, also gives a useful relative bracket.
        double left = minimum * unitQuantile;
        double right = unitQuantile;
        for (int iteration = 0; iteration < 1100; iteration++) {
            double middle = left + (right-left)/2;
            if (middle == left || middle == right || right-left <= 2e-13 * middle)
                return maximum * middle;
            // Near p=1 evaluate the lower tail directly, not 1 - survival.
            boolean above = survivalProbability > .5
                ? lowerProbability(middle, normalized) < 1-survivalProbability
                : survival(middle, normalized).pValue() > survivalProbability;
            if (above) left = middle;
            else right = middle;
        }
        throw new IllegalArgumentException("kernel quantile inversion did not converge");
    }

    private static double lowerProbability(double q, double[] eigenvalues) {
        if (eigenvalues.length == 2) return angularProbability(q, eigenvalues, true);
        try { return gammaSeries(q, eigenvalues, true); }
        catch (SeriesResolutionException unresolved) {
            return 1 - saddlepointSurvival(q, eigenvalues);
        }
    }

    private static double angularProbability(double q, double[] lambda, boolean lower) {
        if (q == 0) return lower ? 0 : 1;
        if (!Double.isFinite(q)) return lower ? 1 : 0;
        double left = 0, right = Math.PI / 2, middle = Math.PI / 4;
        double a = angular(left, q, lambda, lower), b = angular(middle, q, lambda, lower);
        double c = angular(right, q, lambda, lower);
        double integral = angularIntegral(left, right, a, b, c,
            (right-left)*(a+4*b+c)/6, q, lambda, 1e-12, 24, lower);
        return (lower ? -Math.expm1(-q/2) : Math.exp(-q/2)) * (2*integral/Math.PI);
    }

    private static double angular(double angle, double q, double[] lambda, boolean lower) {
        double sine = Math.sin(angle), cosine = Math.cos(angle);
        double variance = lambda[0]*cosine*cosine + lambda[1]*sine*sine;
        if (lower) return -Math.expm1(-q/(2*variance)) / -Math.expm1(-q/2);
        // The eigenvalues are normalized so max(lambda)=1. Remove the peak
        // exp(-q/2) before quadrature: an absolute tolerance on the unscaled
        // integrand is meaningless for small tail probabilities.
        return Math.exp(-q/2 * ((1-variance)/variance));
    }

    private static double angularIntegral(double left, double right, double a,
            double b, double c, double whole, double q, double[] lambda,
            double tolerance, int depth, boolean lower) {
        double middle = (left+right)/2;
        double d = angular((left+middle)/2,q,lambda,lower), e = angular((middle+right)/2,q,lambda,lower);
        double first = (middle-left)*(a+4*d+b)/6;
        double second = (right-middle)*(b+4*e+c)/6;
        double delta = first+second-whole;
        if (Math.abs(delta) <= 15*tolerance) return first+second+delta/15;
        if (depth == 0) throw new IllegalArgumentException("rank-two kernel integration did not converge");
        return angularIntegral(left,middle,a,d,b,first,q,lambda,tolerance/2,depth-1,lower)
            + angularIntegral(middle,right,b,e,c,second,q,lambda,tolerance/2,depth-1,lower);
    }

    /** Positive mixture of beta*chi-square(rank+2*k), beta=min(lambda).
     * The coefficient PGF is G(z)=prod_i[(beta/lambda_i)/(1-r_i*z)]^(1/2),
     * r_i=1-beta/lambda_i. For any 1<z<1/max(r), omitted coefficient mass
     * is <=G(z)/z^(K+1). Gamma tails are <=1, so this also bounds the omitted
     * probability, without cancellation or subtraction from total mass.
     * This is a truncation bound, not a claim of interval-arithmetic rounding.
     */
    private static double gammaSeries(double q, double[] lambda, boolean lower) {
        if (q == 0) return lower ? 0 : 1;
        if (!Double.isFinite(q)) return lower ? 1 : 0;
        double beta = Arrays.stream(lambda).min().orElseThrow();
        double bound = ChiSquare.cumulative(q, lambda.length, lower, false);
        if (!lower && bound == 0) return 0; // rigorous stochastic upper bound underflows
        if (!(beta > 0)) throw new IllegalArgumentException("kernel eigenvalue ratio exceeds numerical range");
        double[] ratios = new double[lambda.length], powers = new double[lambda.length];
        double logC = 0, maximumRatio = 0;
        for (int i = 0; i < lambda.length; i++) {
            ratios[i] = 1-beta/lambda[i];
            powers[i] = 1;
            maximumRatio = Math.max(maximumRatio, ratios[i]);
            logC += .5*Math.log(beta/lambda[i]);
        }
        if (maximumRatio == 0) return ChiSquare.cumulative(q/beta, lambda.length, lower, false);
        double z = 1 + .5*(1/maximumRatio-1), logZ = Math.log(z), logG = logC;
        for (double r : ratios) logG -= .5*Math.log1p(-r*z);
        if (!(logZ > 0) || !Double.isFinite(logG) || Math.exp(logC) == 0)
            throw new SeriesResolutionException();
        double[] weights = new double[MAXIMUM_SERIES_TERMS+1];
        double[] coefficients = new double[MAXIMUM_SERIES_TERMS+1];
        weights[0] = Math.exp(logC);
        double probability = 0;
        for (int k = 0; k <= MAXIMUM_SERIES_TERMS; k++) {
            if (k > 0) {
                double coefficient = 0;
                for (int i = 0; i < lambda.length; i++) {
                    powers[i] *= ratios[i];
                    coefficient += .5*powers[i];
                }
                coefficients[k] = coefficient;
                double sum = 0;
                for (int j = 1; j <= k; j++) sum += coefficients[j]*weights[k-j];
                weights[k] = sum/k;
            }
            probability += weights[k]*ChiSquare.cumulative(q/beta, lambda.length+2.0*k, lower, false);
            // Lower gamma CDF decreases with shape, tightening its remainder.
            double logRemainder = logG-(k+1)*logZ;
            if (lower) logRemainder += ChiSquare.cumulative(q/beta, lambda.length+2.0*(k+1), true, true);
            if (probability > 0 && logRemainder <= Math.log(probability)+Math.log(RELATIVE_SERIES_TOLERANCE))
                return Math.min(1, probability);
        }
        throw new SeriesResolutionException();
    }

    /** Lugannani-Rice tail for the weighted chi-square CGF. This is used only
     * when the explicitly bounded positive series cannot represent its leading
     * coefficient or exhausts its certified remainder budget. */
    private static double saddlepointSurvival(double q, double[] lambda) {
        if (q == 0) return 1;
        if (!Double.isFinite(q)) return 0;
        double mean = Arrays.stream(lambda).sum();
        double upper = .5 / Arrays.stream(lambda).max().orElseThrow();
        double left = q < mean ? -1 : 0;
        while (q < mean && cgfPrime(left, lambda) > q) left *= 2;
        double right = q < mean ? 0 : Math.nextDown(upper);
        double t = 0;
        for (int iteration = 0; q != mean && iteration < 160; iteration++) {
            t = left + (right-left)/2;
            // Compare K'(t)-K'(0) to q-mean without subtracting two
            // nearly equal CGF derivatives at the mean.
            double shift = 0;
            for (double value : lambda) shift += 2*t*value*value/(1-2*t*value);
            boolean below = Math.abs(t) <= .25*upper
                ? shift < q-mean : cgfPrime(t, lambda) < q;
            if (below) left = t; else right = t;
            if (right-left <= 2e-14 * Math.abs(t) || t == 0) break;
        }
        t = q == mean ? 0 : left + (right-left)/2;
        double deviance = 0, second = 0;
        for (double value : lambda) {
            double x = 2*t*value, denominator = 1-x;
            // At the saddle, 2*(t*K'(t)-K(t)) has this form.
            deviance += x/denominator + Math.log1p(-x);
            second += 2*value*value/(denominator*denominator);
        }
        double w, correction;
        if (Math.abs(t) <= .125 / Arrays.stream(lambda).max().orElseThrow()) {
            // A=2*(t*K'(t)-K(t))/t^2 and B=(K''(t)-A)/t.
            // Their convergent series avoid both cancellation in the deviance
            // and the singular subtraction 1/w-1/u. At t=0 this gives the
            // analytic LR limit kappa3/(6*kappa2^(3/2)), continuously.
            double a = 0, b = 0;
            for (double value : lambda) {
                double x = 2*t*value, power = 1;
                for (int j = 0; j < 40; j++) {
                    a += 4*value*value*(j+1.0)/(j+2)*power;
                    b += 4*value*value*value*(j+1.0)*(j+2)/(j+3)*power;
                    power *= x;
                }
            }
            double rootA = Math.sqrt(a), rootSecond = Math.sqrt(second);
            w = t*rootA;
            correction = b/(rootA*rootSecond*(rootA+rootSecond));
        } else {
            w = Math.copySign(Math.sqrt(deviance), t);
            correction = 1/w - 1/(t*Math.sqrt(second));
        }
        double upperNormal = Normal.cumulative(w, 0, 1, false, false);
        double density = Math.exp(-.5*w*w)/Math.sqrt(2*Math.PI);
        double result = upperNormal - density*correction;
        if (!Double.isFinite(result) || result < 0 || result > 1)
            throw new ArithmeticException("saddlepoint tail is outside the probability range");
        return result;
    }

    private static double cgfPrime(double t, double[] lambda) {
        double result = 0;
        for (double value : lambda) result += value/(1-2*t*value);
        return result;
    }

    private static final class SeriesResolutionException
            extends IllegalArgumentException {
        private static final long serialVersionUID = 1L;
    }

    private static double probability(double value) {
        return Math.max(Double.MIN_VALUE, Math.min(1, value));
    }
}
