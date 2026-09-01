/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.internal;

import java.util.Arrays;
import org.jlinalg.model.MissingDataPolicy;

/** Complete-case compaction shared by model front ends. */
public final class CompleteCases {
    private CompleteCases() { }

    public static Data prepare(
            double[] response,
            double[] design,
            int rows,
            int columns,
            double[] weights,
            double[] offset,
            MissingDataPolicy policy) {
        if (response == null || design == null || rows < 1 || columns < 1
                || response.length != rows || design.length != rows * columns) {
            throw new IllegalArgumentException("model data dimensions are invalid");
        }
        if (weights != null && weights.length != rows) {
            throw new IllegalArgumentException("weight length must equal rows");
        }
        if (offset != null && offset.length != rows) {
            throw new IllegalArgumentException("offset length must equal rows");
        }
        if (policy == null) {
            throw new IllegalArgumentException("missing-data policy is required");
        }

        int[] retained = new int[rows];
        int count = 0;
        for (int row = 0; row < rows; row++) {
            boolean finite = Double.isFinite(response[row])
                && (weights == null || Double.isFinite(weights[row]))
                && (offset == null || Double.isFinite(offset[row]));
            for (int column = 0; finite && column < columns; column++) {
                finite = Double.isFinite(design[row * columns + column]);
            }
            if (finite) {
                retained[count++] = row;
            } else if (policy == MissingDataPolicy.ERROR) {
                throw new IllegalArgumentException(
                    "model data contain a non-finite value in row " + row);
            }
        }
        if (count == 0) {
            throw new IllegalArgumentException("no complete observations remain");
        }
        if (count == rows) {
            return new Data(response, design, weights, offset,
                Arrays.copyOf(retained, count), rows);
        }

        double[] compactResponse = new double[count];
        double[] compactDesign = new double[count * columns];
        double[] compactWeights = weights == null ? null : new double[count];
        double[] compactOffset = offset == null ? null : new double[count];
        for (int target = 0; target < count; target++) {
            int source = retained[target];
            compactResponse[target] = response[source];
            System.arraycopy(design, source * columns,
                compactDesign, target * columns, columns);
            if (weights != null) compactWeights[target] = weights[source];
            if (offset != null) compactOffset[target] = offset[source];
        }
        return new Data(compactResponse, compactDesign,
            compactWeights, compactOffset,
            Arrays.copyOf(retained, count), rows);
    }

    public record Data(
            double[] response,
            double[] design,
            double[] weights,
            double[] offset,
            int[] retainedRows,
            int originalRows) {
    }
}
