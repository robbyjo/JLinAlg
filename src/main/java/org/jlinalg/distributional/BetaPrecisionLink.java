/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

/** Precision links supported by R's {@code betareg} package. */
public enum BetaPrecisionLink {
    IDENTITY {
        @Override public double link(double precision) {
            requirePrecision(precision);
            return precision;
        }
        @Override double inverseRaw(double predictor) { return predictor; }
        @Override double derivativeRaw(double predictor, double precision) {
            return 1.0;
        }
    },
    LOG {
        @Override public double link(double precision) {
            requirePrecision(precision);
            return Math.log(precision);
        }
        @Override double inverseRaw(double predictor) {
            return Math.exp(Math.min(700.0, predictor));
        }
        @Override double derivativeRaw(double predictor, double precision) {
            return precision;
        }
    },
    SQRT {
        @Override public double link(double precision) {
            requirePrecision(precision);
            return Math.sqrt(precision);
        }
        @Override double inverseRaw(double predictor) {
            return predictor * predictor;
        }
        @Override double derivativeRaw(double predictor, double precision) {
            return 2.0 * predictor;
        }
    };

    private static final double EPSILON = 1e-12;

    /** Applies the inverse link and returns a positive precision. */
    public double inverse(double predictor) {
        double precision = inverseRaw(predictor);
        return precision > EPSILON && Double.isFinite(precision)
            ? precision : Double.NaN;
    }

    /** Returns d phi / d eta for the supplied predictor. */
    public double derivative(double predictor, double precision) {
        double derivative = derivativeRaw(predictor, precision);
        if (Math.abs(derivative) < Double.MIN_NORMAL) {
            return Math.copySign(Double.MIN_NORMAL,
                derivative == 0.0 ? 1.0 : derivative);
        }
        return derivative;
    }

    /** Applies the link to a positive precision. */
    public abstract double link(double precision);

    abstract double inverseRaw(double predictor);
    abstract double derivativeRaw(double predictor, double precision);

    private static void requirePrecision(double precision) {
        if (!(precision > 0.0) || !Double.isFinite(precision)) {
            throw new IllegalArgumentException("precision must be finite and positive");
        }
    }
}
