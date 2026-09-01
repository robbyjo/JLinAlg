/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pedigree;

import java.util.Objects;
import java.util.Random;
import jdistlib.accelerator.CholeskyFactor;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.mixed.MixedModelSimulationMode;

/** Gaussian animal-model response simulation with pedigree covariance. */
public final class PedigreeSimulation {
    private PedigreeSimulation() { }

    public static double[][] simulate(
            PreparedPedigreeReml model,
            PedigreeRemlResult fitted,
            int simulations,
            long randomSeed,
            MixedModelSimulationMode mode) {
        validate(model, fitted, simulations, mode);
        Random random = new Random(randomSeed);
        double[][] result = new double[simulations][];
        double[] relationshipLower = mode
            == MixedModelSimulationMode.MARGINAL
            ? relationshipLower(model) : null;
        for (int simulation = 0; simulation < simulations; simulation++)
            result[simulation] = simulateOne(model, fitted, random, mode,
                relationshipLower);
        return result;
    }

    static double[] simulateOne(
            PreparedPedigreeReml model,
            PedigreeRemlResult fitted,
            Random random,
            MixedModelSimulationMode mode) {
        return simulateOne(model, fitted, random, mode,
            mode == MixedModelSimulationMode.MARGINAL
                ? relationshipLower(model) : null);
    }

    static double[] simulateOne(
            PreparedPedigreeReml model,
            PedigreeRemlResult fitted,
            Random random,
            MixedModelSimulationMode mode,
            double[] relationshipLower) {
        validate(model, fitted, 1, mode);
        Objects.requireNonNull(random, "random");
        int observations = model.observations();
        int fixedColumns = model.fixedEffectColumns();
        double[] fixed = model.fixedEffectsView();
        double[] beta = fitted.beta();
        double[] response = new double[observations];
        for (int row = 0; row < observations; row++)
            for (int column = 0; column < fixedColumns; column++)
                response[row] += fixed[row * fixedColumns + column]
                    * beta[column];

        double[] breeding;
        if (mode == MixedModelSimulationMode.CONDITIONAL) {
            breeding = fitted.breedingValues();
        } else {
            int animals = model.pedigree().size();
            double[] independent = new double[animals];
            for (int animal = 0; animal < animals; animal++)
                independent[animal] = random.nextGaussian();
            breeding = lowerMultiply(relationshipLower, independent,
                animals, Math.sqrt(fitted.additiveGeneticVariance()));
        }
        int[] observedAnimals = model.observationIndicesView();
        double residualScale = Math.sqrt(fitted.residualVariance());
        for (int row = 0; row < observations; row++)
            response[row] += breeding[observedAnimals[row]]
                + residualScale * random.nextGaussian();
        return response;
    }

    static double[] relationshipLower(
            PreparedPedigreeReml model) {
        try (BackendContext context = BackendContext.select(
                model.backendPolicy())) {
            CholeskyFactor factor = context.backend().dpotrf(
                model.pedigree().relationshipMatrix(),
                model.pedigree().size());
            return factor.lower();
        }
    }

    private static double[] lowerMultiply(
            double[] lower, double[] values, int dimension, double scale) {
        double[] result = new double[dimension];
        for (int row = 0; row < dimension; row++)
            for (int column = 0; column <= row; column++)
                result[row] += scale * lower[row * dimension + column]
                    * values[column];
        return result;
    }

    private static void validate(
            PreparedPedigreeReml model,
            PedigreeRemlResult fitted,
            int simulations,
            MixedModelSimulationMode mode) {
        if (model == null || fitted == null || mode == null)
            throw new IllegalArgumentException(
                "prepared model, fitted result, and simulation mode are required");
        if (simulations < 1)
            throw new IllegalArgumentException("simulations must be positive");
        if (fitted.beta().length != model.fixedEffectColumns()
                || fitted.reml().observations() != model.observations()
                || !fitted.individualIds().equals(
                    model.pedigree().individualIds())
                || !fitted.reml().componentNames().equals(
                    java.util.List.of("additive genetic", "residual")))
            throw new IllegalArgumentException(
                "fitted pedigree result does not match prepared structure");
    }
}
