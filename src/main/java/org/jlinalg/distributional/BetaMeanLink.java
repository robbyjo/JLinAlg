/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import jdistlib.Normal;

/** Mean links supported by R's {@code betareg} package. */
public enum BetaMeanLink {
    LOGIT {
        @Override public double link(double mean) {
            requireMean(mean);
            return Math.log(mean) - Math.log1p(-mean);
        }
        @Override double inverseRaw(double predictor) {
            return predictor >= 0.0
                ? 1.0 / (1.0 + Math.exp(-predictor))
                : Math.exp(predictor) / (1.0 + Math.exp(predictor));
        }
        @Override double derivativeRaw(double predictor, double mean) {
            return mean * (1.0 - mean);
        }
    },
    PROBIT {
        @Override public double link(double mean) {
            requireMean(mean);
            return Normal.quantile(mean, 0.0, 1.0, true, false);
        }
        @Override double inverseRaw(double predictor) {
            return Normal.cumulative(predictor, 0.0, 1.0, true, false);
        }
        @Override double derivativeRaw(double predictor, double mean) {
            return Math.exp(-0.5 * predictor * predictor)
                / Math.sqrt(2.0 * Math.PI);
        }
    },
    CLOGLOG {
        @Override public double link(double mean) {
            requireMean(mean);
            return Math.log(-Math.log1p(-mean));
        }
        @Override double inverseRaw(double predictor) {
            return -Math.expm1(-Math.exp(Math.min(700.0, predictor)));
        }
        @Override double derivativeRaw(double predictor, double mean) {
            double exponential = Math.exp(Math.min(700.0, predictor));
            return exponential * Math.exp(-exponential);
        }
    },
    CAUCHIT {
        @Override public double link(double mean) {
            requireMean(mean);
            return Math.tan(Math.PI * (mean - 0.5));
        }
        @Override double inverseRaw(double predictor) {
            return 0.5 + Math.atan(predictor) / Math.PI;
        }
        @Override double derivativeRaw(double predictor, double mean) {
            return 1.0 / (Math.PI * (1.0 + predictor * predictor));
        }
    },
    LOG {
        @Override public double link(double mean) {
            requireMean(mean);
            return Math.log(mean);
        }
        @Override double inverseRaw(double predictor) {
            return Math.exp(Math.min(0.0, predictor));
        }
        @Override double derivativeRaw(double predictor, double mean) {
            return mean;
        }
    },
    LOGLOG {
        @Override public double link(double mean) {
            requireMean(mean);
            return -Math.log(-Math.log(mean));
        }
        @Override double inverseRaw(double predictor) {
            return Math.exp(-Math.exp(Math.min(700.0, -predictor)));
        }
        @Override double derivativeRaw(double predictor, double mean) {
            double exponential = Math.exp(Math.min(700.0, -predictor));
            return exponential * Math.exp(-exponential);
        }
    };

    private static final double EPSILON = 1e-14;

    /** Applies the inverse link while keeping the result numerically inside (0, 1). */
    public double inverse(double predictor) {
        return Math.max(EPSILON,
            Math.min(1.0 - EPSILON, inverseRaw(predictor)));
    }

    /** Returns d mu / d eta for the supplied predictor. */
    public double derivative(double predictor, double mean) {
        return Math.max(Double.MIN_NORMAL, derivativeRaw(predictor, mean));
    }

    /** Applies the link to a mean in the open unit interval. */
    public abstract double link(double mean);

    abstract double inverseRaw(double predictor);
    abstract double derivativeRaw(double predictor, double mean);

    private static void requireMean(double mean) {
        if (!(mean > 0.0 && mean < 1.0) || !Double.isFinite(mean)) {
            throw new IllegalArgumentException("mean must lie strictly between zero and one");
        }
    }
}
