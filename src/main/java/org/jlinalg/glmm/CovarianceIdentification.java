/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.function.UnaryOperator;
import org.jlinalg.mixed.RandomEffectTerm;

/** Conservative matrix-free rank certificate for induced covariance operators. */
final class CovarianceIdentification {
    private final int rows;
    private final int terms;
    private final double[][] probes;
    private final List<double[]> orthogonal = new ArrayList<>();

    CovarianceIdentification(int rows, int terms) {
        this.rows = rows;
        this.terms = terms;
        // Small problems use the complete coordinate basis. Larger problems
        // use reproducible continuous probes, never an observation covariance.
        int count = rows <= 8 ? rows : 4;
        probes = terms == 1 ? new double[0][] : new double[count][rows];
        Random random = new Random(0x434f56415249414eL);
        for (int p = 0; p < probes.length; p++) {
            if (rows <= 8) probes[p][p] = 1;
            else for (int row = 0; row < rows; row++) probes[p][row] = random.nextGaussian();
        }
    }

    void add(RandomEffectTerm term, UnaryOperator<double[]> precisionSolve) {
        boolean sparse = term.sparse();
        double[] values = sparse ? term.sparseValues() : term.design();
        int[] starts = sparse ? term.rowPointers() : null;
        int[] columns = sparse ? term.columnIndices() : null;
        double scale = 0;
        for (double value : values) scale = Math.max(scale, Math.abs(value));
        if (!(scale > 0)) throw unidentified(term);
        // SPD Q makes a single nonzero Z covariance identifiable structurally.
        if (terms == 1) return;
        for (int i = 0; i < values.length; i++) values[i] /= scale;
        double[] image = new double[Math.multiplyExact(rows, probes.length)];
        int dimension = term.coefficients();
        for (int p = 0; p < probes.length; p++) {
            double[] right = new double[dimension];
            for (int row = 0; row < rows; row++) {
                int begin = sparse ? starts[row] : row * dimension;
                int end = sparse ? starts[row + 1] : begin + dimension;
                for (int i = begin; i < end; i++)
                    right[sparse ? columns[i] : i - begin] += values[i] * probes[p][row];
            }
            double[] solution = precisionSolve.apply(right);
            for (int row = 0; row < rows; row++) {
                int begin = sparse ? starts[row] : row * dimension;
                int end = sparse ? starts[row + 1] : begin + dimension;
                for (int i = begin; i < end; i++)
                    image[p * rows + row] += values[i] * solution[sparse ? columns[i] : i - begin];
            }
        }
        double norm = norm(image);
        if (!(norm > 0) || !Double.isFinite(norm)) throw unidentified(term);
        for (int i = 0; i < image.length; i++) image[i] /= norm;
        // Two-pass modified Gram-Schmidt avoids Gram-matrix cancellation on
        // exactly proportional bases and also detects multi-term dependence.
        for (int pass = 0; pass < 2; pass++) for (double[] previous : orthogonal) {
            double projection = 0;
            for (int i = 0; i < image.length; i++) projection += image[i] * previous[i];
            for (int i = 0; i < image.length; i++) image[i] -= projection * previous[i];
        }
        double remainder = norm(image);
        if (!(remainder > 1e-8)) throw unidentified(term);
        for (int i = 0; i < image.length; i++) image[i] /= remainder;
        orthogonal.add(image);
    }

    private static double norm(double[] values) {
        double result = 0;
        for (double value : values) result = Math.hypot(result, value);
        return result;
    }

    private static IllegalArgumentException unidentified(RandomEffectTerm term) {
        return new IllegalArgumentException("random-effect covariance bases are not independently identifiable "
            + "at numerical precision (term " + term.name() + "); remove redundant terms or reparameterize");
    }
}
