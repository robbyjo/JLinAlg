/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

/** Link used by a GEE association or scale submodel. */
public enum GeeParameterLink {
    IDENTITY {
        @Override public double link(double value) { return value; }
        @Override public double inverse(double value) { return value; }
        @Override public double inverseDerivative(double value) { return 1.0; }
    },
    LOG {
        @Override public double link(double value) {
            return Math.log(Math.max(1e-12, value));
        }
        @Override public double inverse(double value) {
            return Math.exp(Math.max(-30.0, Math.min(30.0, value)));
        }
        @Override public double inverseDerivative(double value) {
            return inverse(value);
        }
    },
    LOGIT {
        @Override public double link(double value) {
            double bounded = Math.max(1e-10, Math.min(1.0 - 1e-10, value));
            return Math.log(bounded / (1.0 - bounded));
        }
        @Override public double inverse(double value) {
            if (value >= 0.0) {
                double exponential = Math.exp(-Math.min(30.0, value));
                return 1.0 / (1.0 + exponential);
            }
            double exponential = Math.exp(Math.max(-30.0, value));
            return exponential / (1.0 + exponential);
        }
        @Override public double inverseDerivative(double value) {
            double probability = inverse(value);
            return probability * (1.0 - probability);
        }
    },
    FISHER_Z {
        @Override public double link(double value) {
            double bounded = Math.max(-0.999999, Math.min(0.999999, value));
            return 0.5 * Math.log((1.0 + bounded) / (1.0 - bounded));
        }
        @Override public double inverse(double value) {
            return Math.tanh(value);
        }
        @Override public double inverseDerivative(double value) {
            double correlation = inverse(value);
            return 1.0 - correlation * correlation;
        }
    };

    public abstract double link(double value);
    public abstract double inverse(double value);
    public abstract double inverseDerivative(double value);
}
