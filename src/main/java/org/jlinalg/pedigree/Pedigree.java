/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pedigree;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * An immutable pedigree and its numerator (additive) relationship matrix.
 * Unknown parents are represented by {@code null} and are assumed unrelated,
 * noninbred founders.
 */
public final class Pedigree {
    private final List<PedigreeIndividual> individuals;
    private final List<String> individualIds;
    private final Map<String, Integer> indexById;
    private final double[] relationshipMatrix;
    private final double[] inbreedingCoefficients;

    private Pedigree(List<PedigreeIndividual> entries) {
        if (entries == null || entries.isEmpty()) {
            throw new IllegalArgumentException("at least one pedigree individual is required");
        }
        this.individuals = List.copyOf(entries);

        LinkedHashMap<String, Integer> inputIndex = new LinkedHashMap<>();
        for (int index = 0; index < individuals.size(); index++) {
            PedigreeIndividual individual = individuals.get(index);
            if (individual == null) {
                throw new IllegalArgumentException("pedigree individuals must not be null");
            }
            if (inputIndex.putIfAbsent(individual.id(), index) != null) {
                throw new IllegalArgumentException(
                    "duplicate pedigree individual: " + individual.id());
            }
        }
        for (PedigreeIndividual individual : individuals) {
            requireKnownParent(inputIndex, individual.id(), individual.sireId());
            requireKnownParent(inputIndex, individual.id(), individual.damId());
        }

        List<Integer> topological = topologicalOrder(individuals, inputIndex);
        double[] topologicalMatrix = tabularRelationship(
            individuals, inputIndex, topological);
        int size = individuals.size();
        int[] topologicalPosition = new int[size];
        for (int position = 0; position < size; position++) {
            topologicalPosition[topological.get(position)] = position;
        }

        this.relationshipMatrix = new double[size * size];
        this.inbreedingCoefficients = new double[size];
        for (int row = 0; row < size; row++) {
            int sourceRow = topologicalPosition[row];
            for (int column = 0; column < size; column++) {
                int sourceColumn = topologicalPosition[column];
                relationshipMatrix[row * size + column] =
                    topologicalMatrix[sourceRow * size + sourceColumn];
            }
            inbreedingCoefficients[row] = relationshipMatrix[row * size + row] - 1.0;
        }
        this.indexById = Map.copyOf(inputIndex);
        this.individualIds = List.copyOf(inputIndex.keySet());
    }

    /** Constructs and validates a pedigree in any input order. */
    public static Pedigree of(List<PedigreeIndividual> individuals) {
        return new Pedigree(individuals);
    }

    /** Returns the entries in caller-supplied order. */
    public List<PedigreeIndividual> individuals() {
        return individuals;
    }

    /** Returns individual identifiers in matrix row/column order. */
    public List<String> individualIds() {
        return individualIds;
    }

    /** Returns the number of individuals. */
    public int size() {
        return individuals.size();
    }

    /** Returns an individual's zero-based row/column index. */
    public int indexOf(String individualId) {
        Integer index = indexById.get(individualId);
        if (index == null) {
            throw new IllegalArgumentException(
                "individual is absent from the pedigree: " + individualId);
        }
        return index;
    }

    /** Returns the additive relationship between two individuals. */
    public double relationship(String firstId, String secondId) {
        return relationshipMatrix[indexOf(firstId) * size() + indexOf(secondId)];
    }

    /** Returns a defensive row-major copy of the numerator relationship matrix. */
    public double[] relationshipMatrix() {
        return relationshipMatrix.clone();
    }

    /** Returns inbreeding coefficients in {@link #individualIds()} order. */
    public double[] inbreedingCoefficients() {
        return inbreedingCoefficients.clone();
    }

    /**
     * Constructs the inverse numerator relationship matrix directly in sparse
     * CSR form using the pedigree {@code T D T'} decomposition.
     */
    public SparseSymmetricMatrix sparseRelationshipMatrixInverse() {
        int size = size();
        List<Map<Integer, Double>> rows = new ArrayList<>(size);
        for (int row = 0; row < size; row++) {
            rows.add(new TreeMap<>());
        }
        for (int individual = 0; individual < size; individual++) {
            PedigreeIndividual entry = individuals.get(individual);
            int sire = entry.sireId() == null ? -1 : indexOf(entry.sireId());
            int dam = entry.damId() == null ? -1 : indexOf(entry.damId());
            double diagonal = mendelianSamplingVariance(sire, dam);
            double inverseDiagonal = 1.0 / diagonal;
            int[] indices = {individual, sire, dam};
            double[] coefficients = {1.0, -0.5, -0.5};
            for (int left = 0; left < indices.length; left++) {
                if (indices[left] < 0) {
                    continue;
                }
                for (int right = 0; right < indices.length; right++) {
                    if (indices[right] < 0) {
                        continue;
                    }
                    double contribution = coefficients[left]
                        * coefficients[right] * inverseDiagonal;
                    rows.get(indices[left]).merge(
                        indices[right], contribution, Double::sum);
                }
            }
        }

        int nonzeros = 0;
        for (Map<Integer, Double> row : rows) {
            nonzeros += row.size();
        }
        int[] rowPointers = new int[size + 1];
        int[] columnIndices = new int[nonzeros];
        double[] values = new double[nonzeros];
        int position = 0;
        for (int row = 0; row < size; row++) {
            rowPointers[row] = position;
            for (Map.Entry<Integer, Double> entry : rows.get(row).entrySet()) {
                columnIndices[position] = entry.getKey();
                values[position] = entry.getValue();
                position++;
            }
        }
        rowPointers[size] = position;
        return new SparseSymmetricMatrix(
            size, rowPointers, columnIndices, values);
    }

    double relationship(int firstIndex, int secondIndex) {
        return relationshipMatrix[firstIndex * size() + secondIndex];
    }

    private double mendelianSamplingVariance(int sire, int dam) {
        if (sire < 0 && dam < 0) {
            return 1.0;
        }
        if (sire < 0 || dam < 0) {
            int parent = sire >= 0 ? sire : dam;
            return 0.75 - 0.25 * inbreedingCoefficients[parent];
        }
        return 0.5 - 0.25 * (
            inbreedingCoefficients[sire] + inbreedingCoefficients[dam]);
    }

    private static void requireKnownParent(
            Map<String, Integer> inputIndex, String child, String parent) {
        if (parent != null && !inputIndex.containsKey(parent)) {
            throw new IllegalArgumentException(
                "parent " + parent + " of " + child + " is absent from the pedigree");
        }
    }

    private static List<Integer> topologicalOrder(
            List<PedigreeIndividual> individuals,
            Map<String, Integer> inputIndex) {
        int size = individuals.size();
        int[] indegree = new int[size];
        List<List<Integer>> children = new ArrayList<>(size);
        for (int index = 0; index < size; index++) {
            children.add(new ArrayList<>());
        }
        for (int child = 0; child < size; child++) {
            PedigreeIndividual individual = individuals.get(child);
            if (individual.sireId() != null) {
                int parent = inputIndex.get(individual.sireId());
                indegree[child]++;
                children.get(parent).add(child);
            }
            if (individual.damId() != null) {
                int parent = inputIndex.get(individual.damId());
                indegree[child]++;
                children.get(parent).add(child);
            }
        }

        ArrayDeque<Integer> ready = new ArrayDeque<>();
        for (int index = 0; index < size; index++) {
            if (indegree[index] == 0) {
                ready.add(index);
            }
        }
        List<Integer> order = new ArrayList<>(size);
        while (!ready.isEmpty()) {
            int parent = ready.remove();
            order.add(parent);
            for (int child : children.get(parent)) {
                indegree[child]--;
                if (indegree[child] == 0) {
                    ready.add(child);
                }
            }
        }
        if (order.size() != size) {
            throw new IllegalArgumentException("pedigree contains an ancestry cycle");
        }
        return order;
    }

    private static double[] tabularRelationship(
            List<PedigreeIndividual> individuals,
            Map<String, Integer> inputIndex,
            List<Integer> topological) {
        int size = individuals.size();
        Map<Integer, Integer> positionByInputIndex = new HashMap<>();
        double[] matrix = new double[size * size];
        for (int position = 0; position < size; position++) {
            int input = topological.get(position);
            PedigreeIndividual individual = individuals.get(input);
            int sire = parentPosition(
                individual.sireId(), inputIndex, positionByInputIndex);
            int dam = parentPosition(
                individual.damId(), inputIndex, positionByInputIndex);

            for (int previous = 0; previous < position; previous++) {
                double value = 0.5 * (
                    relationshipOrZero(matrix, size, sire, previous)
                    + relationshipOrZero(matrix, size, dam, previous));
                matrix[position * size + previous] = value;
                matrix[previous * size + position] = value;
            }
            matrix[position * size + position] = 1.0
                + 0.5 * relationshipOrZero(matrix, size, sire, dam);
            positionByInputIndex.put(input, position);
        }
        return matrix;
    }

    private static int parentPosition(
            String parentId,
            Map<String, Integer> inputIndex,
            Map<Integer, Integer> positionByInputIndex) {
        if (parentId == null) {
            return -1;
        }
        Integer position = positionByInputIndex.get(inputIndex.get(parentId));
        if (position == null) {
            throw new IllegalStateException("parent was not topologically ordered");
        }
        return position;
    }

    private static double relationshipOrZero(
            double[] matrix, int size, int row, int column) {
        if (row < 0 || column < 0) {
            return 0.0;
        }
        return matrix[row * size + column];
    }
}
