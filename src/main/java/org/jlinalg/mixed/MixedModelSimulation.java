/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mixed;

import java.util.List;
import java.util.Objects;
import java.util.Random;

/** Reproducible {@code simulate.merMod}-style Gaussian LMM simulation. */
public final class MixedModelSimulation {
    private MixedModelSimulation() { }

    public static double[][] simulate(
            PreparedLinearMixedModel model,
            LinearMixedModelResult fitted,
            int simulations,
            long randomSeed,
            MixedModelSimulationMode mode) {
        validate(model, fitted, simulations, mode);
        Random random = new Random(randomSeed);
        double[][] result = new double[simulations][];
        for (int simulation = 0; simulation < simulations; simulation++)
            result[simulation] = simulateOne(model, fitted, random, mode);
        return result;
    }

    static double[] simulateOne(
            PreparedLinearMixedModel model,
            LinearMixedModelResult fitted,
            Random random,
            MixedModelSimulationMode mode) {
        validate(model, fitted, 1, mode);
        Objects.requireNonNull(random, "random");
        int observations = model.observations();
        int fixedColumns = model.fixedEffectColumns();
        double[] fixed = model.fixedEffectsView();
        double[] beta = fitted.beta();
        double[] variances = fitted.reml().varianceComponents();
        double[] response = new double[observations];
        for (int row = 0; row < observations; row++)
            for (int column = 0; column < fixedColumns; column++)
                response[row] += fixed[row * fixedColumns + column]
                    * beta[column];

        List<RandomEffectTerm> terms = model.randomEffects();
        for (int termIndex = 0; termIndex < terms.size(); termIndex++) {
            RandomEffectTerm term = terms.get(termIndex);
            double[] coefficients;
            if (mode == MixedModelSimulationMode.CONDITIONAL) {
                coefficients = fitted.randomEffects(term.name()).estimates();
            } else {
                coefficients = new double[term.coefficients()];
                double scale = Math.sqrt(variances[termIndex]);
                for (int coefficient = 0;
                        coefficient < coefficients.length; coefficient++)
                    coefficients[coefficient] = scale * random.nextGaussian();
            }
            addTerm(response, term, coefficients);
        }
        double residualScale = Math.sqrt(variances[terms.size()]);
        for (int row = 0; row < observations; row++)
            response[row] += residualScale * random.nextGaussian();
        return response;
    }

    private static void validate(
            PreparedLinearMixedModel model,
            LinearMixedModelResult fitted,
            int simulations,
            MixedModelSimulationMode mode) {
        if (model == null || fitted == null || mode == null)
            throw new IllegalArgumentException(
                "prepared model, fitted result, and simulation mode are required");
        if (simulations < 1)
            throw new IllegalArgumentException("simulations must be positive");
        if (fitted.beta().length != model.fixedEffectColumns()
                || fitted.reml().observations() != model.observations()
                || fitted.reml().varianceComponents().length
                    != model.randomEffects().size() + 1)
            throw new IllegalArgumentException(
                "fitted result does not match the prepared model structure");
        List<String> names = fitted.reml().componentNames();
        for (int term = 0; term < model.randomEffects().size(); term++)
            if (!names.get(term).equals(model.randomEffects().get(term).name()))
                throw new IllegalArgumentException(
                    "fitted random-effect order does not match prepared model");
        if (!names.get(names.size() - 1).equals("residual"))
            throw new IllegalArgumentException(
                "fitted result has no matching residual component");
    }

    private static void addTerm(
            double[] response, RandomEffectTerm term, double[] coefficients) {
        if (coefficients.length != term.coefficients())
            throw new IllegalArgumentException(
                "random-effect coefficient count does not match term "
                    + term.name());
        if (term.sparse()) {
            int[] pointers = term.rowPointersView();
            int[] columns = term.columnIndicesView();
            double[] values = term.sparseValuesView();
            for (int row = 0; row < term.observations(); row++)
                for (int index = pointers[row]; index < pointers[row + 1]; index++)
                    response[row] += values[index] * coefficients[columns[index]];
        } else {
            double[] design = term.denseDesignView();
            for (int row = 0; row < term.observations(); row++)
                for (int column = 0; column < term.coefficients(); column++)
                    response[row] += design[row * term.coefficients() + column]
                        * coefficients[column];
        }
    }
}
