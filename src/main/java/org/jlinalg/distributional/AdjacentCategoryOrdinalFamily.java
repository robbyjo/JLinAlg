/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** VGAM-style adjacent-category ordinal logits with one predictor per cut. */
public final class AdjacentCategoryOrdinalFamily implements DistributionalFamily {
    private final int categories;
    private final int predictors;
    private final List<String> names;

    public AdjacentCategoryOrdinalFamily(int categories) {
        if (categories < 2) {
            throw new IllegalArgumentException("ordinal family needs at least two categories");
        }
        this.categories = categories;
        this.predictors = categories - 1;
        List<String> values = new ArrayList<>(predictors);
        for (int cut = 0; cut < predictors; cut++) values.add("adjacentLogit" + cut);
        this.names = List.copyOf(values);
    }
    @Override public String name() { return "ordinal-acat(" + categories + ")"; }
    @Override public int parameterCount() { return predictors; }
    @Override public List<String> parameterNames() { return names; }
    @Override public void validateResponse(double response) {
        if (response != Math.rint(response) || response < 0.0
                || response >= categories) {
            throw new IllegalArgumentException(
                "ordinal responses must be integer category indices");
        }
    }
    @Override public double[] initialParameters(double[] response) {
        double[] counts = new double[categories];
        Arrays.fill(counts, 0.5);
        for (double value : response) counts[(int) value]++;
        double[] result = new double[predictors];
        for (int cut = 0; cut < predictors; cut++) {
            result[cut] = Math.log(counts[cut] / counts[cut + 1]);
        }
        return result;
    }
    @Override public double link(int parameter, double value) {
        check(parameter);
        return value;
    }
    @Override public double inverseLink(int parameter, double predictor) {
        check(parameter);
        return predictor;
    }
    @Override public double logLikelihood(double response, double[] parameters) {
        Probabilities probabilities = probabilities(parameters);
        return Math.log(probabilities.values()[(int) response]);
    }
    @Override public void derivatives(
            double response, double[] parameters,
            double[] score, double[] information) {
        double[] probabilities = probabilities(parameters).values();
        double[] cumulative = new double[predictors];
        double running = 0.0;
        for (int category = 0; category < predictors; category++) {
            running += probabilities[category];
            cumulative[category] = running;
        }
        int observed = (int) response;
        for (int first = 0; first < predictors; first++) {
            score[first] = (observed <= first ? 1.0 : 0.0) - cumulative[first];
            for (int second = 0; second < predictors; second++) {
                information[first * predictors + second] =
                    cumulative[Math.min(first, second)]
                        - cumulative[first] * cumulative[second];
            }
        }
    }
    /** Returns fitted category probabilities from adjacent logits. */
    @Override
    public double[] categoryProbabilities(double[] adjacentLogits) {
        return probabilities(adjacentLogits).values().clone();
    }
    private Probabilities probabilities(double[] parameters) {
        if (parameters.length != predictors) {
            throw new IllegalArgumentException("wrong number of ordinal predictors");
        }
        double[] logWeights = new double[categories];
        double running = 0.0;
        for (int category = predictors - 1; category >= 0; category--) {
            running += parameters[category];
            logWeights[category] = running;
        }
        double maximum = Arrays.stream(logWeights).max().orElse(0.0);
        double denominator = 0.0;
        double[] probabilities = new double[categories];
        for (int category = 0; category < categories; category++) {
            probabilities[category] = Math.exp(logWeights[category] - maximum);
            denominator += probabilities[category];
        }
        for (int category = 0; category < categories; category++) {
            probabilities[category] /= denominator;
        }
        return new Probabilities(probabilities);
    }
    private void check(int parameter) {
        if (parameter < 0 || parameter >= predictors) {
            throw new IllegalArgumentException("unknown ordinal predictor");
        }
    }
    private record Probabilities(double[] values) { }
}
