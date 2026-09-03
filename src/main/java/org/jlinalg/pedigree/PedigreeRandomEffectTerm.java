/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.pedigree;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparsePrecisionMatrix;

/** A pedigree incidence term paired with its sparse additive precision. */
public record PedigreeRandomEffectTerm(
        RandomEffectTerm randomEffect,
        SparsePrecisionMatrix precision) {

    public PedigreeRandomEffectTerm {
        Objects.requireNonNull(randomEffect, "randomEffect");
        Objects.requireNonNull(precision, "precision");
        if (randomEffect.coefficients() != precision.dimension())
            throw new IllegalArgumentException(
                "pedigree random-effect and precision dimensions must match");
    }

    /** Creates a pedigree term, retaining unobserved ancestors as coefficients. */
    public static PedigreeRandomEffectTerm of(
            String name, List<String> observationIndividualIds,
            Pedigree pedigree) {
        if (pedigree == null || observationIndividualIds == null
                || observationIndividualIds.isEmpty())
            throw new IllegalArgumentException(
                "pedigree and observation identifiers are required");
        int rows = observationIndividualIds.size();
        int[] rowStarts = new int[rows + 1];
        int[] columns = new int[rows];
        double[] values = new double[rows];
        for (int row = 0; row < rows; row++) {
            rowStarts[row] = row;
            columns[row] = pedigree.indexOf(observationIndividualIds.get(row));
            values[row] = 1.0;
        }
        rowStarts[rows] = rows;
        RandomEffectTerm term = RandomEffectTerm.ofSparseCsr(name, rows,
            pedigree.size(), rowStarts, columns, values, pedigree.individualIds());
        SparseSymmetricMatrix inverse = pedigree.sparseRelationshipMatrixInverse();
        return new PedigreeRandomEffectTerm(term, new SparsePrecisionMatrix(
            inverse.dimension(), inverse.rowPointers(),
            inverse.columnIndices(), inverse.values()));
    }

    /**
     * Creates a pedigree term directly from entries and known inbreeding
     * coefficients, without materializing the dense relationship matrix.
     */
    public static PedigreeRandomEffectTerm ofSparse(
            String name, List<String> observationIndividualIds,
            List<PedigreeIndividual> individuals,
            double[] inbreedingCoefficients) {
        if (individuals == null || individuals.isEmpty()
                || observationIndividualIds == null
                || observationIndividualIds.isEmpty()
                || inbreedingCoefficients == null
                || inbreedingCoefficients.length != individuals.size())
            throw new IllegalArgumentException(
                "pedigree entries, observations, and inbreeding are required");
        Map<String, Integer> indexById = new LinkedHashMap<>();
        List<String> coefficientNames = new ArrayList<>(individuals.size());
        for (int index = 0; index < individuals.size(); index++) {
            PedigreeIndividual value = individuals.get(index);
            if (value == null || indexById.putIfAbsent(value.id(), index) != null)
                throw new IllegalArgumentException(
                    "pedigree entries must be nonnull and uniquely named");
            coefficientNames.add(value.id());
            double inbreeding = inbreedingCoefficients[index];
            if (!Double.isFinite(inbreeding) || inbreeding < 0
                    || inbreeding >= 1)
                throw new IllegalArgumentException(
                    "inbreeding coefficients must be finite in [0, 1)");
        }
        int rows = observationIndividualIds.size();
        int[] incidenceStarts = new int[rows + 1];
        int[] incidenceColumns = new int[rows];
        double[] incidenceValues = new double[rows];
        for (int row = 0; row < rows; row++) {
            Integer coefficient = indexById.get(observationIndividualIds.get(row));
            if (coefficient == null)
                throw new IllegalArgumentException(
                    "observed individual is absent from pedigree: "
                        + observationIndividualIds.get(row));
            incidenceStarts[row] = row;
            incidenceColumns[row] = coefficient;
            incidenceValues[row] = 1.0;
        }
        incidenceStarts[rows] = rows;
        RandomEffectTerm term = RandomEffectTerm.ofSparseCsr(name, rows,
            individuals.size(), incidenceStarts, incidenceColumns,
            incidenceValues, coefficientNames);

        List<Map<Integer, Double>> precisionRows =
            new ArrayList<>(individuals.size());
        for (int row = 0; row < individuals.size(); row++)
            precisionRows.add(new TreeMap<>());
        for (int individual = 0; individual < individuals.size(); individual++) {
            PedigreeIndividual value = individuals.get(individual);
            int sire = parent(indexById, value.id(), value.sireId());
            int dam = parent(indexById, value.id(), value.damId());
            double mendelian;
            if (sire < 0 && dam < 0) mendelian = 1.0;
            else if (sire < 0 || dam < 0) {
                int known = sire >= 0 ? sire : dam;
                mendelian = 0.75 - 0.25 * inbreedingCoefficients[known];
            } else mendelian = 0.5 - 0.25
                * (inbreedingCoefficients[sire] + inbreedingCoefficients[dam]);
            double inverse = 1.0 / mendelian;
            int[] indices = {individual, sire, dam};
            double[] coefficients = {1.0, -0.5, -0.5};
            for (int left = 0; left < indices.length; left++) {
                if (indices[left] < 0) continue;
                for (int right = 0; right < indices.length; right++) {
                    if (indices[right] < 0) continue;
                    precisionRows.get(indices[left]).merge(indices[right],
                        coefficients[left] * coefficients[right] * inverse,
                        Double::sum);
                }
            }
        }
        int nonzeros = precisionRows.stream().mapToInt(Map::size).sum();
        int[] precisionStarts = new int[individuals.size() + 1];
        int[] precisionColumns = new int[nonzeros];
        double[] precisionValues = new double[nonzeros];
        int position = 0;
        for (int row = 0; row < individuals.size(); row++) {
            precisionStarts[row] = position;
            for (Map.Entry<Integer, Double> entry
                    : precisionRows.get(row).entrySet()) {
                precisionColumns[position] = entry.getKey();
                precisionValues[position++] = entry.getValue();
            }
        }
        precisionStarts[individuals.size()] = position;
        return new PedigreeRandomEffectTerm(term, new SparsePrecisionMatrix(
            individuals.size(), precisionStarts, precisionColumns,
            precisionValues));
    }

    /** Direct sparse construction for pedigrees known to be noninbred. */
    public static PedigreeRandomEffectTerm ofUninbred(
            String name, List<String> observationIndividualIds,
            List<PedigreeIndividual> individuals) {
        return ofSparse(name, observationIndividualIds, individuals,
            new double[individuals.size()]);
    }

    private static int parent(
            Map<String, Integer> indexById, String child, String parent) {
        if (parent == null) return -1;
        Integer index = indexById.get(parent);
        if (index == null)
            throw new IllegalArgumentException(
                "parent " + parent + " of " + child + " is absent");
        return index;
    }
}
