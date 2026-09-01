/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pedigree;

import java.util.List;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.Reml;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.RemlResult;
import org.jlinalg.reml.VarianceComponent;

/** Dense animal-model REML using a pedigree numerator relationship matrix. */
public final class PedigreeReml {
    private static final String GENETIC_COMPONENT = "additive genetic";
    private static final String RESIDUAL_COMPONENT = "residual";

    private PedigreeReml() {
    }

    /** Fits pedigree REML using default options and the preferred backend. */
    public static PedigreeRemlResult fit(
            double[] response,
            double[][] fixedEffects,
            List<String> observationIndividualIds,
            Pedigree pedigree) {
        return fit(response, fixedEffects, observationIndividualIds, pedigree,
            RemlOptions.defaults(), BackendPolicy.PREFERRED);
    }

    /** Fits pedigree REML from a conventional rectangular design matrix. */
    public static PedigreeRemlResult fit(
            double[] response,
            double[][] fixedEffects,
            List<String> observationIndividualIds,
            Pedigree pedigree,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        if (response == null) {
            throw new IllegalArgumentException("response is required");
        }
        double[] rowMajor = MatrixOps.rowMajor(fixedEffects, response.length);
        return fit(response, rowMajor, response.length, fixedEffects[0].length,
            observationIndividualIds, pedigree, options, backendPolicy);
    }

    /** Fits pedigree REML from a contiguous row-major design matrix. */
    public static PedigreeRemlResult fit(
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            List<String> observationIndividualIds,
            Pedigree pedigree,
            RemlOptions options,
            BackendPolicy backendPolicy) {
        MatrixOps.validateModelData(response, fixedEffects, rows, columns);
        if (pedigree == null || options == null || backendPolicy == null) {
            throw new IllegalArgumentException(
                "pedigree, options, and backendPolicy are required");
        }
        int[] observationIndices = observationIndices(
            observationIndividualIds, rows, pedigree);
        double[] observationRelationship = observationRelationship(
            observationIndices, pedigree);
        List<VarianceComponent> components = List.of(
            new VarianceComponent(GENETIC_COMPONENT, rows, observationRelationship),
            VarianceComponent.identity(RESIDUAL_COMPONENT, rows));

        try (BackendContext context = BackendContext.select(backendPolicy)) {
            RemlResult reml = Reml.fitWithKnownCovariance(
                response, fixedEffects, rows, columns, components,
                new double[rows * rows], options, context);
            return predictions(reml, response, fixedEffects, rows, columns,
                observationIndices, pedigree,
                context.backend());
        }
    }

    private static PedigreeRemlResult predictions(
            RemlResult reml,
            double[] response,
            double[] fixedEffects,
            int rows,
            int columns,
            int[] observationIndices,
            Pedigree pedigree,
            ComputeBackend backend) {
        double genetic = reml.varianceComponents()[0];
        double residual = reml.varianceComponents()[1];
        int animals = pedigree.size();
        int dimension = columns + animals;
        double inverseResidual = 1.0 / residual;
        double inverseGenetic = 1.0 / genetic;
        double[] equations = new double[dimension * dimension];
        double[] rightSide = new double[dimension];
        for (int observation = 0; observation < rows; observation++) {
            int animalColumn = columns + observationIndices[observation];
            double weightedResponse = inverseResidual * response[observation];
            rightSide[animalColumn] += weightedResponse;
            equations[animalColumn * dimension + animalColumn] += inverseResidual;
            for (int fixed = 0; fixed < columns; fixed++) {
                double value = fixedEffects[observation * columns + fixed];
                rightSide[fixed] += value * weightedResponse;
                double cross = inverseResidual * value;
                equations[fixed * dimension + animalColumn] += cross;
                equations[animalColumn * dimension + fixed] += cross;
                for (int second = 0; second < columns; second++) {
                    equations[fixed * dimension + second] += inverseResidual
                        * value * fixedEffects[observation * columns + second];
                }
            }
        }
        SparseSymmetricMatrix relationshipInverse =
            pedigree.sparseRelationshipMatrixInverse();
        int[] pointers = relationshipInverse.rowPointers();
        int[] indices = relationshipInverse.columnIndices();
        double[] values = relationshipInverse.values();
        for (int animal = 0; animal < animals; animal++) {
            for (int index = pointers[animal]; index < pointers[animal + 1]; index++) {
                equations[(columns + animal) * dimension
                    + columns + indices[index]] += inverseGenetic * values[index];
            }
        }
        CholeskyFactor factor = backend.dpotrf(equations, dimension);
        double[] solution = factor.solve(rightSide);
        double[] breedingValues = new double[animals];
        System.arraycopy(solution, columns, breedingValues, 0, animals);

        double[] predictionErrorVariances = new double[animals];
        double[] reliabilities = new double[animals];
        double[] coefficientCovariance = factor.solve(
            MatrixOps.identity(dimension), dimension);
        double[] inbreeding = pedigree.inbreedingCoefficients();
        for (int animal = 0; animal < animals; animal++) {
            int coefficient = columns + animal;
            double priorVariance = genetic * (1.0 + inbreeding[animal]);
            double pev = clamp(coefficientCovariance[
                coefficient * dimension + coefficient], 0.0, priorVariance);
            predictionErrorVariances[animal] = pev;
            reliabilities[animal] = priorVariance > 0.0
                ? clamp(1.0 - pev / priorVariance, 0.0, 1.0) : 0.0;
        }

        return new PedigreeRemlResult(reml, pedigree.individualIds(),
            breedingValues, predictionErrorVariances, reliabilities);
    }

    private static int[] observationIndices(
            List<String> individualIds, int rows, Pedigree pedigree) {
        if (individualIds == null || individualIds.size() != rows) {
            throw new IllegalArgumentException(
                "one pedigree individual identifier is required per observation");
        }
        int[] result = new int[rows];
        for (int index = 0; index < rows; index++) {
            String id = individualIds.get(index);
            if (id == null) {
                throw new IllegalArgumentException(
                    "observation individual identifiers must not be null");
            }
            result[index] = pedigree.indexOf(id);
        }
        return result;
    }

    private static double[] observationRelationship(
            int[] observationIndices, Pedigree pedigree) {
        int rows = observationIndices.length;
        double[] result = new double[rows * rows];
        for (int row = 0; row < rows; row++) {
            for (int column = 0; column <= row; column++) {
                double value = pedigree.relationship(
                    observationIndices[row], observationIndices[column]);
                result[row * rows + column] = value;
                result[column * rows + row] = value;
            }
        }
        return result;
    }

    private static double clamp(double value, double minimum, double maximum) {
        return Math.max(minimum, Math.min(maximum, value));
    }
}
