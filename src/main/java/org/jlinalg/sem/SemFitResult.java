/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import java.util.List;
import org.jlinalg.compute.BackendProvenance;

/** Covariance-structure ML estimates and conventional SEM fit indices. */
public final class SemFitResult {
    private final List<SemParameterEstimate> parameters;
    private final double[] impliedCovariance;
    private final double logLikelihood;
    private final double chiSquare;
    private final int degreesOfFreedom;
    private final double pValue;
    private final double cfi;
    private final double tli;
    private final double rmsea;
    private final double srmr;
    private final double aic;
    private final double bic;
    private final int observations;
    private final int evaluations;
    private final boolean converged;
    private final BackendProvenance backend;

    SemFitResult(List<SemParameterEstimate> parameters, double[] impliedCovariance,
            double logLikelihood, double chiSquare, int degreesOfFreedom,
            double pValue, double cfi, double tli, double rmsea, double srmr,
            double aic, double bic, int observations, int evaluations,
            boolean converged, BackendProvenance backend) {
        this.parameters = List.copyOf(parameters);
        this.impliedCovariance = impliedCovariance.clone();
        this.logLikelihood = logLikelihood;
        this.chiSquare = chiSquare;
        this.degreesOfFreedom = degreesOfFreedom;
        this.pValue = pValue;
        this.cfi = cfi;
        this.tli = tli;
        this.rmsea = rmsea;
        this.srmr = srmr;
        this.aic = aic;
        this.bic = bic;
        this.observations = observations;
        this.evaluations = evaluations;
        this.converged = converged;
        this.backend = backend;
    }
    public List<SemParameterEstimate> parameters() { return parameters; }
    public SemParameterEstimate parameter(String label) {
        return parameters.stream().filter(value -> value.label().equals(label))
            .findFirst().orElseThrow(() -> new IllegalArgumentException("unknown SEM parameter: " + label));
    }
    public double[] impliedCovariance() { return impliedCovariance.clone(); }
    public double logLikelihood() { return logLikelihood; }
    public double chiSquare() { return chiSquare; }
    public int degreesOfFreedom() { return degreesOfFreedom; }
    public double pValue() { return pValue; }
    public double cfi() { return cfi; }
    public double tli() { return tli; }
    public double rmsea() { return rmsea; }
    public double srmr() { return srmr; }
    public double aic() { return aic; }
    public double bic() { return bic; }
    public int observations() { return observations; }
    public int functionEvaluations() { return evaluations; }
    public boolean converged() { return converged; }
    public BackendProvenance backend() { return backend; }
}
