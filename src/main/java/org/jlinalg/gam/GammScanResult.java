/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.gam;

import java.util.LinkedHashMap;
import java.util.Map;
import org.jlinalg.mixed.LinearMixedModelResult;
import org.jlinalg.reml.RemlResult;

/**
 * Sparse coefficient-space GAMM result for a changing-predictor scan.
 *
 * <p>The smooth EDF is reported exactly from its targeted prediction-error
 * variances. A total model EDF is intentionally not reported because doing so
 * for a large pedigree term would require selected inversion over every
 * pedigree coefficient.</p>
 */
public final class GammScanResult {
    private final LinearMixedModelResult mixedModel;
    private final SmoothTermEstimate smoothTerm;
    private final Map<String, double[]> randomContributions;

    GammScanResult(
            LinearMixedModelResult mixedModel,
            SmoothTermEstimate smoothTerm,
            Map<String, double[]> randomContributions) {
        this.mixedModel = mixedModel;
        this.smoothTerm = smoothTerm;
        Map<String, double[]> copied = new LinkedHashMap<>();
        randomContributions.forEach(
            (name, values) -> copied.put(name, values.clone()));
        this.randomContributions =
            java.util.Collections.unmodifiableMap(copied);
    }

    public LinearMixedModelResult mixedModel() { return mixedModel; }
    public RemlResult reml() { return mixedModel.reml(); }
    public SmoothTermEstimate smoothTerm() { return smoothTerm; }
    public Map<String, double[]> randomContributions() {
        Map<String, double[]> copied = new LinkedHashMap<>();
        randomContributions.forEach(
            (name, values) -> copied.put(name, values.clone()));
        return java.util.Collections.unmodifiableMap(copied);
    }
    public double[] randomContribution(String name) {
        double[] result = randomContributions.get(name);
        if (result == null) {
            throw new IllegalArgumentException(
                "unknown random-effect term: " + name);
        }
        return result.clone();
    }
    public double[] fittedValues() { return mixedModel.fittedValues(); }
    public double[] residuals() { return mixedModel.residuals(); }
}
