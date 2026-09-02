/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Result of a multi-predictor penalized likelihood fit. */
public final class DistributionalResult {
    private final String family;
    private final List<DistributionalParameterResult> parameters;
    private final Map<String, DistributionalParameterResult> byName;
    private final double[] covariance;
    private final double logLikelihood;
    private final double penalizedLogLikelihood;
    private final int iterations;
    private final boolean converged;
    private final String convergenceMessage;

    DistributionalResult(
            String family,
            List<DistributionalParameterResult> parameters,
            double[] covariance,
            double logLikelihood,
            double penalizedLogLikelihood,
            int iterations,
            boolean converged,
            String convergenceMessage) {
        this.family = family;
        this.parameters = List.copyOf(parameters);
        Map<String, DistributionalParameterResult> map = new LinkedHashMap<>();
        for (DistributionalParameterResult parameter : parameters) {
            map.put(parameter.name(), parameter);
        }
        this.byName = Map.copyOf(map);
        this.covariance = covariance.clone();
        this.logLikelihood = logLikelihood;
        this.penalizedLogLikelihood = penalizedLogLikelihood;
        this.iterations = iterations;
        this.converged = converged;
        this.convergenceMessage = convergenceMessage;
    }

    public String family() { return family; }
    public List<DistributionalParameterResult> parameters() {
        return parameters;
    }
    public DistributionalParameterResult parameter(String name) {
        DistributionalParameterResult result = byName.get(name);
        if (result == null) {
            throw new IllegalArgumentException("unknown parameter: " + name);
        }
        return result;
    }
    public double[] covariance() { return covariance.clone(); }
    public double logLikelihood() { return logLikelihood; }
    public double penalizedLogLikelihood() { return penalizedLogLikelihood; }
    public int iterations() { return iterations; }
    public boolean converged() { return converged; }
    public String convergenceMessage() { return convergenceMessage; }
}
