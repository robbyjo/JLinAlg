/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;

/** Pedigree-correlated Cox frailty estimates and fixed-effect hazard ratios. */
public final class CoxPedigreeResult {
    private final CoxMixedResult mixed;
    private final List<String> individualIds;
    private final Map<String, Double> frailties;

    CoxPedigreeResult(CoxMixedResult mixed, List<String> individualIds) {
        this.mixed = mixed;
        this.individualIds = List.copyOf(individualIds);
        double[] modes = mixed.randomEffects("pedigree").modes();
        Map<String, Double> values = new LinkedHashMap<>();
        for (int index = 0; index < modes.length; index++)
            values.put(this.individualIds.get(index), modes[index]);
        frailties = java.util.Collections.unmodifiableMap(values);
    }

    public CoxMixedResult mixedModel() { return mixed; }
    public double[] beta() { return mixed.beta(); }
    public double[] fixef() { return beta(); }
    public double[] hazardRatios() { return mixed.hazardRatios(); }
    public double[] standardErrors() { return mixed.standardErrors(); }
    public double[] zStatistics() { return mixed.zStatistics(); }
    public double[] pValues() { return mixed.pValues(); }
    public double[] negativeLog10PValues() {
        return mixed.negativeLog10PValues();
    }
    public AssociationStatistics associationStatistics() {
        return mixed.associationStatistics();
    }
    public double frailtyVariance() {
        return mixed.randomEffects("pedigree").variance();
    }
    public List<String> individualIds() { return individualIds; }
    public Map<String, Double> ranef() { return frailties; }
    public double frailty(String individualId) {
        Double value = frailties.get(individualId);
        if (value == null)
            throw new IllegalArgumentException(
                "individual is absent from Cox pedigree result: " + individualId);
        return value;
    }
    public List<BaselineHazardPoint> baselineHazard() {
        return mixed.baselineHazard();
    }
    public double laplaceLogLikelihood() {
        return mixed.laplaceLogLikelihood();
    }
    public boolean converged() { return mixed.converged(); }
    public String convergenceMessage() { return mixed.convergenceMessage(); }
    public BackendProvenance backend() { return mixed.backend(); }
}
