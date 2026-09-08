/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glm;

import jdistlib.Binomial;
import jdistlib.Gamma;
import jdistlib.NegBinomial;
import jdistlib.Normal;
import jdistlib.Poisson;

/** Built-in canonical GLM families. */
public final class GlmFamilies {
    private static final GlmFamily GAUSSIAN = new GaussianIdentity();
    private static final GlmFamily BINOMIAL = new BinomialLogit();
    private static final GlmFamily POISSON = new PoissonLog();
    private static final GlmFamily GAMMA = new GammaLog();
    private static final GlmFamily INVERSE_GAUSSIAN = new InverseGaussianLog();
    private static final GlmFamily QUASI_BINOMIAL = new QuasiBinomialLogit();
    private static final GlmFamily QUASI_POISSON = new QuasiPoissonLog();

    private GlmFamilies() {
    }

    /** Gaussian family with identity link and estimated dispersion. */
    public static GlmFamily gaussian() { return GAUSSIAN; }

    /** Binomial proportions with logit link; prior weights represent trials. */
    public static GlmFamily binomial() { return BINOMIAL; }

    /** Poisson counts with log link. */
    public static GlmFamily poisson() { return POISSON; }

    /** Gamma responses with a log link and estimated dispersion. */
    public static GlmFamily gamma() { return GAMMA; }

    /** Inverse-Gaussian responses with a log link and estimated dispersion. */
    public static GlmFamily inverseGaussian() { return INVERSE_GAUSSIAN; }

    /** Negative-binomial counts with log link and caller-supplied size. */
    public static GlmFamily negativeBinomial(double size) {
        return new NegativeBinomialLog(size);
    }
    /** Tweedie quasi-likelihood with log link and caller-supplied variance power. */
    public static GlmFamily tweedie(double power) {
        return new TweedieLogFamily(power);
    }

    /** Quasi-binomial proportions with logit link and Pearson dispersion. */
    public static GlmFamily quasiBinomial() { return QUASI_BINOMIAL; }

    /** Quasi-Poisson nonnegative responses with estimated dispersion. */
    public static GlmFamily quasiPoisson() { return QUASI_POISSON; }

    private static final class GaussianIdentity implements GlmFamily {
        @Override public String name() { return "gaussian(identity)"; }
        @Override public void validateResponse(double response, double priorWeight) { }
        @Override public double initialMean(double response) { return response; }
        @Override public double link(double mean) { return mean; }
        @Override public double inverseLink(double predictor) { return predictor; }
        @Override public double meanDerivative(double predictor) { return 1.0; }
        @Override public double variance(double mean) { return 1.0; }
        @Override public double unitDeviance(double response, double mean) {
            double residual = response - mean;
            return residual * residual;
        }
        @Override public double logLikelihood(
                double response, double mean, double priorWeight, double dispersion) {
            double standardDeviation = Math.sqrt(dispersion / priorWeight);
            return Normal.density(response, mean, standardDeviation, true);
        }
        @Override public boolean fixedDispersion() { return false; }
    }

    private static class BinomialLogit implements GlmFamily {
        private static final double EPSILON = Double.MIN_NORMAL;
        private static final double MAXIMUM = Math.nextDown(1.0);

        @Override public String name() { return "binomial(logit)"; }
        @Override public void validateResponse(double response, double priorWeight) {
            if (response < 0.0 || response > 1.0) {
                throw new IllegalArgumentException(
                    "binomial responses must lie between zero and one");
            }
        }
        @Override public double initialMean(double response) {
            return clamp((response + 0.5) / 2.0, EPSILON, MAXIMUM);
        }
        @Override public double link(double mean) {
            double bounded = clamp(mean, EPSILON, MAXIMUM);
            return Math.log(bounded) - Math.log1p(-bounded);
        }
        @Override public double inverseLink(double predictor) {
            double result;
            if (predictor >= 0.0) {
                double exponential = Math.exp(-predictor);
                result = 1.0 / (1.0 + exponential);
            } else {
                double exponential = Math.exp(predictor);
                result = exponential / (1.0 + exponential);
            }
            return result;
        }
        @Override public double meanDerivative(double predictor) {
            double tail = Math.exp(-Math.abs(predictor));
            return tail / ((1.0 + tail) * (1.0 + tail));
        }
        @Override public double varianceAtPredictor(double predictor, double mean) {
            return meanDerivative(predictor);
        }
        @Override public double residualAtPredictor(double response, double predictor, double mean) {
            double tail = Math.exp(-Math.abs(predictor));
            double small = tail / (1.0 + tail), large = 1.0 / (1.0 + tail);
            return predictor >= 0.0 ? response * small - (1.0 - response) * large
                : response * large - (1.0 - response) * small;
        }
        @Override public double workingWeight(double response, double predictor, double mean, double weight) {
            return weight * meanDerivative(predictor);
        }
        @Override public double unitDevianceAtPredictor(double response, double predictor, double mean) {
            double logTail = -Math.log1p(Math.exp(-Math.abs(predictor)));
            double logP = predictor >= 0.0 ? logTail : predictor + logTail;
            double logQ = predictor >= 0.0 ? -predictor + logTail : logTail;
            double first = response == 0.0 ? 0.0 : response * (Math.log(response) - logP);
            double second = response == 1.0 ? 0.0 : (1.0-response) * (Math.log1p(-response) - logQ);
            return Math.max(0.0, 2.0 * (first + second));
        }
        @Override public double logLikelihoodAtPredictor(double response, double predictor, double mean,
                double weight, double dispersion) {
            // Dynamic dispatch retains quasi-binomial's undefined likelihood.
            double ordinary = logLikelihood(response, mean, weight, dispersion);
            if (Double.isNaN(ordinary) || predictor < 0.0) return ordinary;
            return Binomial.density(Math.rint(weight) - Math.rint(response * weight),
                Math.rint(weight), inverseLink(-predictor), true);
        }
        @Override public double variance(double mean) {
            return mean * (1.0 - mean);
        }
        @Override public double unitDeviance(double response, double mean) {
            if (mean == response) return 0.0;
            if (mean == 0.0 || mean == 1.0) return Double.POSITIVE_INFINITY;
            double bounded = clamp(mean, EPSILON, MAXIMUM);
            double first = response == 0.0 ? 0.0
                : response * Math.log(response / bounded);
            double second = response == 1.0 ? 0.0
                : (1.0 - response) * (Math.log1p(-response) - Math.log1p(-bounded));
            return 2.0 * (first + second);
        }
        @Override public double logLikelihood(
                double response, double mean, double priorWeight, double dispersion) {
            double successes = response * priorWeight;
            double roundedTrials = Math.rint(priorWeight);
            double roundedSuccesses = Math.rint(successes);
            if (Math.abs(priorWeight - roundedTrials) > 1e-9
                    || Math.abs(successes - roundedSuccesses) > 1e-9) {
                return Double.NaN;
            }
            return Binomial.density(
                roundedSuccesses, roundedTrials, mean, true);
        }
        @Override public boolean fixedDispersion() { return true; }
    }

    private static class PoissonLog implements GlmFamily {
        private static final double MINIMUM_MEAN = Double.MIN_NORMAL;
        private static final double MAXIMUM_PREDICTOR = 700.0;

        @Override public String name() { return "poisson(log)"; }
        @Override public void validateResponse(double response, double priorWeight) {
            if (response < 0.0 || response != Math.rint(response)) {
                throw new IllegalArgumentException(
                    "Poisson responses must be nonnegative integers");
            }
        }
        @Override public double initialMean(double response) {
            return Math.max(0.1, response + 0.1);
        }
        @Override public double link(double mean) {
            return Math.log(Math.max(MINIMUM_MEAN, mean));
        }
        @Override public double inverseLink(double predictor) {
            return Math.max(MINIMUM_MEAN,
                Math.exp(Math.min(MAXIMUM_PREDICTOR, predictor)));
        }
        @Override public double meanDerivative(double predictor) {
            return inverseLink(predictor);
        }
        @Override public double variance(double mean) {
            return Math.max(MINIMUM_MEAN, mean);
        }
        @Override public double unitDeviance(double response, double mean) {
            if (response == 0.0) {
                return 2.0 * mean;
            }
            return 2.0 * (response * Math.log(response / mean) - (response - mean));
        }
        @Override public double logLikelihood(
                double response, double mean, double priorWeight, double dispersion) {
            return priorWeight * Poisson.density(response, mean, true);
        }
        @Override public boolean fixedDispersion() { return true; }
    }

    private static final class QuasiBinomialLogit extends BinomialLogit {
        @Override public String name() { return "quasibinomial(logit)"; }
        @Override public double logLikelihood(
                double response, double mean, double priorWeight, double dispersion) {
            return Double.NaN;
        }
        @Override public boolean fixedDispersion() { return false; }
    }

    private static final class QuasiPoissonLog extends PoissonLog {
        @Override public String name() { return "quasipoisson(log)"; }
        @Override public void validateResponse(double response, double priorWeight) {
            if (response < 0.0) {
                throw new IllegalArgumentException(
                    "quasi-Poisson responses must be nonnegative");
            }
        }
        @Override public double logLikelihood(
                double response, double mean, double priorWeight, double dispersion) {
            return Double.NaN;
        }
        @Override public boolean fixedDispersion() { return false; }
    }

    private static final class GammaLog implements GlmFamily {
        private static final double MINIMUM = Double.MIN_NORMAL;
        @Override public String name() { return "Gamma(log)"; }
        @Override public void validateResponse(double response, double priorWeight) {
            if (!(response > 0.0)) {
                throw new IllegalArgumentException("Gamma responses must be positive");
            }
        }
        @Override public double initialMean(double response) {
            return Math.max(MINIMUM, response);
        }
        @Override public double link(double mean) {
            return Math.log(Math.max(MINIMUM, mean));
        }
        @Override public double inverseLink(double predictor) {
            return Math.max(MINIMUM, Math.exp(Math.min(700.0, predictor)));
        }
        @Override public double meanDerivative(double predictor) {
            return inverseLink(predictor);
        }
        @Override public double variance(double mean) {
            return Math.max(MINIMUM, mean * mean);
        }
        @Override public double unitDeviance(double response, double mean) {
            double ratio = response / mean;
            return 2.0 * (ratio - 1.0 - Math.log(ratio));
        }
        @Override public double logLikelihood(
                double response, double mean, double priorWeight, double dispersion) {
            double shape = priorWeight / dispersion;
            double scale = mean / shape;
            return Gamma.density(response, shape, scale, true);
        }
        @Override public boolean fixedDispersion() { return false; }
    }

    private static final class InverseGaussianLog implements GlmFamily {
        private static final double MINIMUM = Double.MIN_NORMAL;
        @Override public String name() { return "inverse.gaussian(log)"; }
        @Override public void validateResponse(double response, double priorWeight) {
            if (!(response > 0.0)) {
                throw new IllegalArgumentException(
                    "inverse-Gaussian responses must be positive");
            }
        }
        @Override public double initialMean(double response) {
            return Math.max(MINIMUM, response);
        }
        @Override public double link(double mean) {
            return Math.log(Math.max(MINIMUM, mean));
        }
        @Override public double inverseLink(double predictor) {
            return Math.max(MINIMUM, Math.exp(Math.min(700.0, predictor)));
        }
        @Override public double meanDerivative(double predictor) {
            return inverseLink(predictor);
        }
        @Override public double variance(double mean) {
            return Math.max(MINIMUM, mean * mean * mean);
        }
        @Override public double unitDeviance(double response, double mean) {
            double residual = response - mean;
            return residual * residual / (response * mean * mean);
        }
        @Override public double logLikelihood(
                double response, double mean, double priorWeight, double dispersion) {
            double residual = response - mean;
            return 0.5 * Math.log(priorWeight
                / (2.0 * Math.PI * dispersion * response * response * response))
                - priorWeight * residual * residual
                    / (2.0 * dispersion * response * mean * mean);
        }
        @Override public boolean fixedDispersion() { return false; }
    }

    private static final class NegativeBinomialLog implements GlmFamily {
        private static final double MINIMUM = Double.MIN_NORMAL;
        private final double size;

        NegativeBinomialLog(double size) {
            if (!(size > 0.0) || !Double.isFinite(size)) {
                throw new IllegalArgumentException(
                    "negative-binomial size must be finite and positive");
            }
            this.size = size;
        }

        @Override public String name() { return "negative.binomial(log; size=" + size + ")"; }
        @Override public void validateResponse(double response, double priorWeight) {
            if (response < 0.0 || response != Math.rint(response)) {
                throw new IllegalArgumentException(
                    "negative-binomial responses must be nonnegative integers");
            }
        }
        @Override public double initialMean(double response) {
            return Math.max(0.1, response + 0.1);
        }
        @Override public double link(double mean) {
            return Math.log(Math.max(MINIMUM, mean));
        }
        @Override public double inverseLink(double predictor) {
            return Math.max(MINIMUM, Math.exp(Math.min(700.0, predictor)));
        }
        @Override public double meanDerivative(double predictor) {
            return inverseLink(predictor);
        }
        @Override public double variance(double mean) {
            return Math.max(MINIMUM, mean + mean * mean / size);
        }
        @Override public double unitDeviance(double response, double mean) {
            double first = response == 0.0 ? 0.0 : response * Math.log(response / mean);
            double second = (response + size)
                * Math.log((response + size) / (mean + size));
            return 2.0 * (first - second);
        }
        @Override public double logLikelihood(
                double response, double mean, double priorWeight, double dispersion) {
            return priorWeight * NegBinomial.density_mu(response, size, mean, true);
        }
        @Override public boolean fixedDispersion() { return true; }
    }

    private static double clamp(double value, double minimum, double maximum) {
        return Math.max(minimum, Math.min(maximum, value));
    }

    // Pearson dispersion is for covariance, not a maximum-likelihood estimate.
    // Unknown/custom families retain their supplied-dispersion likelihood contract.
    static double likelihoodDispersion(GlmFamily family, double[] y, double[] mu,
            double[] weights, double pearson, double deviance) {
        if (family.fixedDispersion()) return 1.0;
        if (family instanceof GaussianIdentity || family instanceof InverseGaussianLog) {
            return Math.max(Double.MIN_NORMAL, deviance / y.length);
        }
        if (!(family instanceof GammaLog)) return pearson;
        double center = Math.log(pearson);
        double left = center - 4.0;
        double right = center + 4.0;
        for (int i = 0; i < 16; i++) {
            double middle = likelihood(family, y, mu, weights, center);
            if (likelihood(family, y, mu, weights, left) > middle) {
                right = center; center = left; left -= 4.0;
            } else if (likelihood(family, y, mu, weights, right) > middle) {
                left = center; center = right; right += 4.0;
            } else break;
        }
        double ratio = (Math.sqrt(5.0) - 1.0) / 2.0;
        double a = right - ratio * (right - left);
        double b = left + ratio * (right - left);
        double fa = likelihood(family, y, mu, weights, a);
        double fb = likelihood(family, y, mu, weights, b);
        for (int i = 0; i < 90; i++) {
            if (fa > fb) {
                right = b; b = a; fb = fa; a = right - ratio * (right - left);
                fa = likelihood(family, y, mu, weights, a);
            } else {
                left = a; a = b; fa = fb; b = left + ratio * (right - left);
                fb = likelihood(family, y, mu, weights, b);
            }
        }
        return Math.exp((left + right) / 2.0);
    }

    private static double likelihood(GlmFamily family, double[] y, double[] mu,
            double[] weights, double logDispersion) {
        double value = 0.0;
        double dispersion = Math.exp(logDispersion);
        for (int row = 0; row < y.length; row++) {
            value += family.logLikelihood(y[row], mu[row], weights[row], dispersion);
        }
        return Double.isFinite(value) ? value : Double.NEGATIVE_INFINITY;
    }
}
