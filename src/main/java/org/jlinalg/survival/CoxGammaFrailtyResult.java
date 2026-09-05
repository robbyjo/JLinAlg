/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;

/** Laplace-profiled multiplicative shared gamma-frailty Cox fit. */
public final class CoxGammaFrailtyResult {
    private final AssociationStatistics statistics;
    private final double[] covariance;
    private final double frailtyVariance;
    private final List<String> groupNames;
    private final double[] logFrailties;
    private final double[] frailties;
    private final Map<String, Double> frailtyByGroup;
    private final List<BaselineHazardPoint> baselineHazard;
    private final double partialLogLikelihood;
    private final double laplaceLogLikelihood;
    private final int varianceIterations;
    private final boolean converged;
    private final String convergenceMessage;
    private final BackendProvenance backend;

    CoxGammaFrailtyResult(double[] beta, double[] covariance,
            double frailtyVariance, List<String> groupNames,
            double[] logFrailties, List<BaselineHazardPoint> baselineHazard,
            double partialLogLikelihood, double laplaceLogLikelihood,
            int varianceIterations, boolean converged,
            String convergenceMessage, BackendProvenance backend) {
        this.covariance = covariance.clone();
        double[] standardErrors = new double[beta.length];
        for (int column = 0; column < beta.length; column++)
            standardErrors[column] = Math.sqrt(Math.max(0.0,
                covariance[column * beta.length + column]));
        statistics = AssociationStatistics.normal(beta, standardErrors);
        this.frailtyVariance = frailtyVariance;
        this.groupNames = List.copyOf(groupNames);
        this.logFrailties = logFrailties.clone();
        this.frailties = java.util.Arrays.stream(logFrailties)
            .map(Math::exp).toArray();
        Map<String, Double> indexed = new LinkedHashMap<>();
        for (int index = 0; index < groupNames.size(); index++)
            indexed.put(groupNames.get(index), frailties[index]);
        frailtyByGroup = Map.copyOf(indexed);
        this.baselineHazard = List.copyOf(baselineHazard);
        this.partialLogLikelihood = partialLogLikelihood;
        this.laplaceLogLikelihood = laplaceLogLikelihood;
        this.varianceIterations = varianceIterations;
        this.converged = converged;
        this.convergenceMessage = convergenceMessage;
        this.backend = backend;
    }

    public double[] beta() { return statistics.beta(); }
    public double[] standardErrors() { return statistics.standardErrors(); }
    public double[] zStatistics() { return statistics.statistics(); }
    public double[] pValues() { return statistics.pValues(); }
    public double[] covariance() { return covariance.clone(); }
    public double frailtyVariance() { return frailtyVariance; }
    public List<String> groupNames() { return groupNames; }
    public double[] logFrailties() { return logFrailties.clone(); }
    public double[] frailties() { return frailties.clone(); }
    public Map<String, Double> frailtyByGroup() { return frailtyByGroup; }
    public double frailty(String group) {
        Double result = frailtyByGroup.get(group);
        if (result == null)
            throw new IllegalArgumentException("unknown frailty group: " + group);
        return result;
    }
    public List<BaselineHazardPoint> baselineHazard() { return baselineHazard; }
    public double partialLogLikelihood() { return partialLogLikelihood; }
    public double laplaceLogLikelihood() { return laplaceLogLikelihood; }
    public int varianceIterations() { return varianceIterations; }
    public boolean converged() { return converged; }
    public String convergenceMessage() { return convergenceMessage; }
    public BackendProvenance backend() { return backend; }
}
