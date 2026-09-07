/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.nonlinear;

/** User-supplied nonlinear mean and analytic fixed-parameter gradient. */
@FunctionalInterface
public interface NonlinearMeanFunction {
    Evaluation evaluate(double[] parameters, int row);

    /** One row of a nonlinear mean evaluation. */
    record Evaluation(double value, double[] gradient) {
        public Evaluation {
            if (!Double.isFinite(value) || gradient == null)
                throw new IllegalArgumentException("nonlinear evaluation is invalid");
            gradient = gradient.clone();
            for (double element : gradient) if (!Double.isFinite(element))
                throw new IllegalArgumentException("nonlinear gradient must be finite");
        }
        @Override public double[] gradient() { return gradient.clone(); }
    }
}
