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
        private static final double EPSILON = 1e-12;

        @Override public String name() { return "binomial(logit)"; }
        @Override public void validateResponse(double response, double priorWeight) {
            if (response < 0.0 || response > 1.0) {
                throw new IllegalArgumentException(
                    "binomial responses must lie between zero and one");
            }
        }
        @Override public double initialMean(double response) {
            return clamp((response + 0.5) / 2.0, EPSILON, 1.0 - EPSILON);
        }
        @Override public double link(double mean) {
            double bounded = clamp(mean, EPSILON, 1.0 - EPSILON);
            return Math.log(bounded / (1.0 - bounded));
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
            return clamp(result, EPSILON, 1.0 - EPSILON);
        }
        @Override public double meanDerivative(double predictor) {
            double mean = inverseLink(predictor);
            return mean * (1.0 - mean);
        }
        @Override public double variance(double mean) {
            return Math.max(EPSILON, mean * (1.0 - mean));
        }
        @Override public double unitDeviance(double response, double mean) {
            double bounded = clamp(mean, EPSILON, 1.0 - EPSILON);
            double first = response == 0.0 ? 0.0
                : response * Math.log(response / bounded);
            double second = response == 1.0 ? 0.0
                : (1.0 - response) * Math.log((1.0 - response) / (1.0 - bounded));
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
        private static final double MINIMUM_MEAN = 1e-12;
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
        private static final double MINIMUM = 1e-12;
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
        private static final double MINIMUM = 1e-12;
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
        private static final double MINIMUM = 1e-12;
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
}
