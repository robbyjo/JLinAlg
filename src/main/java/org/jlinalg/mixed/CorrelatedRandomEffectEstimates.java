/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

import java.util.List;

/** Covariance and conditional modes for one correlated random block. */
public final class CorrelatedRandomEffectEstimates {
    private final String name;
    private final List<String> groupNames;
    private final List<String> effectNames;
    private final double[] covariance;
    private final double[] modes;

    CorrelatedRandomEffectEstimates(
            String name, List<String> groupNames, List<String> effectNames,
            double[] covariance, double[] modes) {
        this.name = name;
        this.groupNames = List.copyOf(groupNames);
        this.effectNames = List.copyOf(effectNames);
        this.covariance = covariance.clone();
        this.modes = modes.clone();
    }

    public String name() { return name; }
    public List<String> groupNames() { return groupNames; }
    public List<String> effectNames() { return effectNames; }
    public double[] covariance() { return covariance.clone(); }
    /** Group-major modes: all effects for group 1, then group 2, and so on. */
    public double[] modes() { return modes.clone(); }
    public double correlation(int firstEffect, int secondEffect) {
        int effects = effectNames.size();
        double covarianceValue = covariance[firstEffect * effects + secondEffect];
        return covarianceValue / Math.sqrt(
            covariance[firstEffect * effects + firstEffect]
                * covariance[secondEffect * effects + secondEffect]);
    }
}
