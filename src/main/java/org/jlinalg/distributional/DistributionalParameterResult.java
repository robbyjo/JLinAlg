/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

/** Coefficient inference and fitted values for one distribution parameter. */
public final class DistributionalParameterResult {
    private final String name;
    private final double[] coefficients;
    private final double[] standardErrors;
    private final double[] statistics;
    private final double[] pValues;
    private final double[] fittedValues;
    private final double effectiveDegreesOfFreedom;

    DistributionalParameterResult(
            String name,
            double[] coefficients,
            double[] standardErrors,
            double[] statistics,
            double[] pValues,
            double[] fittedValues,
            double effectiveDegreesOfFreedom) {
        this.name = name;
        this.coefficients = coefficients.clone();
        this.standardErrors = standardErrors.clone();
        this.statistics = statistics.clone();
        this.pValues = pValues.clone();
        this.fittedValues = fittedValues.clone();
        this.effectiveDegreesOfFreedom = effectiveDegreesOfFreedom;
    }

    public String name() { return name; }
    public double[] coefficients() { return coefficients.clone(); }
    public double[] standardErrors() { return standardErrors.clone(); }
    public double[] statistics() { return statistics.clone(); }
    public double[] pValues() { return pValues.clone(); }
    public double[] fittedValues() { return fittedValues.clone(); }
    public double effectiveDegreesOfFreedom() {
        return effectiveDegreesOfFreedom;
    }
}
