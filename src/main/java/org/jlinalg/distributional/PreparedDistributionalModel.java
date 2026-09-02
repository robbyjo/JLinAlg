/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.gam.PenalizedPredictor;

/** Reusable multi-predictor model with coefficient warm starts for repeated responses. */
public final class PreparedDistributionalModel {
    private final List<PenalizedPredictor> predictors;
    private final DistributionalFamily family;
    private final DistributionalOptions options;
    private final BackendPolicy backendPolicy;

    public PreparedDistributionalModel(
            List<PenalizedPredictor> predictors,
            DistributionalFamily family,
            DistributionalOptions options,
            BackendPolicy backendPolicy) {
        if (predictors == null || predictors.isEmpty()) {
            throw new IllegalArgumentException("predictors are required");
        }
        this.predictors = List.copyOf(predictors);
        this.family = Objects.requireNonNull(family, "family");
        this.options = Objects.requireNonNull(options, "options");
        this.backendPolicy = Objects.requireNonNull(backendPolicy, "backendPolicy");
        if (predictors.size() != family.parameterCount()) {
            throw new IllegalArgumentException("one predictor is required per parameter");
        }
    }

    public DistributionalResult fit(double[] response) {
        return DistributionalModel.fit(response, predictors, family,
            options, backendPolicy);
    }

    /** Reuses all bases and penalties and starts at a preceding fit's coefficients. */
    public DistributionalResult refit(
            DistributionalResult previous, double[] response) {
        Objects.requireNonNull(previous, "previous");
        List<double[]> starting = new ArrayList<>(previous.parameters().size());
        for (DistributionalParameterResult parameter : previous.parameters()) {
            starting.add(parameter.coefficients());
        }
        return DistributionalModel.fit(response, predictors, family,
            options, backendPolicy, starting);
    }

    public List<PenalizedPredictor> predictors() { return predictors; }
    public DistributionalFamily family() { return family; }
}
