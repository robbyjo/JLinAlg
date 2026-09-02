/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.Arrays;
import java.util.Objects;
import org.jlinalg.internal.MatrixOps;

/** Unique-row P-spline storage for BAM-like repeated-covariate execution. */
public final class DiscretePSplineBasis {
    private final PSplineTerm term;
    private final double[] uniqueValues;
    private final int[] rowMap;
    private final double[] uniqueBasis;

    private DiscretePSplineBasis(
            PSplineTerm term,
            double[] uniqueValues,
            int[] rowMap,
            double[] uniqueBasis) {
        this.term = term;
        this.uniqueValues = uniqueValues;
        this.rowMap = rowMap;
        this.uniqueBasis = uniqueBasis;
    }

    /** Discretizes an already compiled P-spline without changing its knots. */
    public static DiscretePSplineBasis compile(PSplineTerm term) {
        Objects.requireNonNull(term, "term");
        double[] values = term.covariate();
        double[] unique = Arrays.stream(values).distinct().sorted().toArray();
        int[] rowMap = new int[values.length];
        for (int row = 0; row < values.length; row++) {
            rowMap[row] = Arrays.binarySearch(unique, values[row]);
        }
        return new DiscretePSplineBasis(term, unique, rowMap,
            term.basis(unique));
    }

    public int observations() { return rowMap.length; }
    public int uniqueRows() { return uniqueValues.length; }
    public int columns() { return term.basisDimension(); }
    public double compressionRatio() {
        return uniqueValues.length / (double) rowMap.length;
    }
    public double[] uniqueValues() { return uniqueValues.clone(); }
    public int[] rowMap() { return rowMap.clone(); }
    public double[] uniqueBasis() { return uniqueBasis.clone(); }

    /** Expands to the ordinary observation-by-basis row-major matrix. */
    public double[] expand() {
        int columns = columns();
        double[] result = new double[observations() * columns];
        for (int row = 0; row < observations(); row++) {
            System.arraycopy(uniqueBasis, rowMap[row] * columns,
                result, row * columns, columns);
        }
        return result;
    }

    /** Forms B'WB by aggregating weights at unique covariate values. */
    public double[] crossProduct(double[] weights) {
        double[] checked = weights(weights);
        int columns = columns();
        double[] aggregated = new double[uniqueRows()];
        for (int row = 0; row < observations(); row++) {
            aggregated[rowMap[row]] += checked[row];
        }
        double[] result = new double[columns * columns];
        for (int unique = 0; unique < uniqueRows(); unique++) {
            int offset = unique * columns;
            for (int first = 0; first < columns; first++) {
                double left = uniqueBasis[offset + first] * aggregated[unique];
                for (int second = 0; second <= first; second++) {
                    result[first * columns + second] +=
                        left * uniqueBasis[offset + second];
                }
            }
        }
        for (int first = 0; first < columns; first++) {
            for (int second = 0; second < first; second++) {
                result[second * columns + first] =
                    result[first * columns + second];
            }
        }
        return result;
    }

    /** Forms B'Wy without expanding repeated basis rows. */
    public double[] transposeMultiply(double[] response, double[] weights) {
        if (response == null || response.length != observations()) {
            throw new IllegalArgumentException("response length must match observations");
        }
        MatrixOps.requireFinite(response, "response");
        double[] checked = weights(weights);
        double[] aggregated = new double[uniqueRows()];
        for (int row = 0; row < observations(); row++) {
            aggregated[rowMap[row]] += checked[row] * response[row];
        }
        int columns = columns();
        double[] result = new double[columns];
        for (int unique = 0; unique < uniqueRows(); unique++) {
            for (int column = 0; column < columns; column++) {
                result[column] += uniqueBasis[unique * columns + column]
                    * aggregated[unique];
            }
        }
        return result;
    }

    private double[] weights(double[] supplied) {
        if (supplied == null) {
            double[] result = new double[observations()];
            Arrays.fill(result, 1.0);
            return result;
        }
        if (supplied.length != observations()) {
            throw new IllegalArgumentException("weight length must match observations");
        }
        double[] result = MatrixOps.finiteCopy(supplied, "weights");
        for (double value : result) {
            if (value < 0.0) {
                throw new IllegalArgumentException("weights must be nonnegative");
            }
        }
        return result;
    }
}
