/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.Arrays;
import java.util.List;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.internal.MatrixOps;

/** Compiled additive predictor with a diagonalized quadratic penalty. */
public final class PenalizedPredictor {
    private final double[] design;
    private final double[] penaltyDiagonal;
    private final int observations;
    private final int columns;
    private final int parametricColumns;

    PenalizedPredictor(
            double[] design,
            double[] penaltyDiagonal,
            int observations,
            int columns,
            int parametricColumns) {
        this.design = design;
        this.penaltyDiagonal = penaltyDiagonal;
        this.observations = observations;
        this.columns = columns;
        this.parametricColumns = parametricColumns;
    }

    /** Creates an unpenalized linear predictor. */
    public static PenalizedPredictor linear(double[][] design) {
        if (design == null || design.length == 0 || design[0] == null) {
            throw new IllegalArgumentException("predictor design is required");
        }
        int rows = design.length;
        int columns = design[0].length;
        return new PenalizedPredictor(MatrixOps.rowMajor(design, rows),
            new double[columns], rows, columns, columns);
    }

    /**
     * Compiles P-spline terms using fixed positive smoothing parameters.
     * Null-space columns are unpenalized and penalized columns are whitened.
     */
    public static PenalizedPredictor additive(
            double[][] parametricDesign,
            List<PSplineTerm> smoothTerms,
            double[] smoothingParameters,
            BackendPolicy backendPolicy) {
        if (parametricDesign == null || parametricDesign.length == 0
                || parametricDesign[0] == null || smoothTerms == null
                || smoothingParameters == null
                || smoothingParameters.length != smoothTerms.size()
                || backendPolicy == null) {
            throw new IllegalArgumentException(
                "parametric design, smooths, smoothing parameters, and backend are required");
        }
        int rows = parametricDesign.length;
        int parametricColumns = parametricDesign[0].length;
        double[] parametric = MatrixOps.rowMajor(parametricDesign, rows);
        for (double value : smoothingParameters) {
            if (!(value > 0.0) || !Double.isFinite(value)) {
                throw new IllegalArgumentException(
                    "smoothing parameters must be finite and positive");
            }
        }
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            PSplineMixedModelCompiler.Compiled compiled =
                PSplineMixedModelCompiler.compile(parametric, rows,
                    parametricColumns, smoothTerms, context.backend());
            int columns = compiled.fixedColumns();
            for (PSplineMixedModelCompiler.Term term : compiled.terms()) {
                columns += term.randomColumns();
            }
            double[] design = new double[rows * columns];
            for (int row = 0; row < rows; row++) {
                System.arraycopy(compiled.fixedDesign(),
                    row * compiled.fixedColumns(), design, row * columns,
                    compiled.fixedColumns());
            }
            double[] penalty = new double[columns];
            int destination = compiled.fixedColumns();
            for (int index = 0; index < compiled.terms().size(); index++) {
                PSplineMixedModelCompiler.Term term = compiled.terms().get(index);
                for (int row = 0; row < rows; row++) {
                    System.arraycopy(term.randomDesign(),
                        row * term.randomColumns(),
                        design, row * columns + destination,
                        term.randomColumns());
                }
                Arrays.fill(penalty, destination,
                    destination + term.randomColumns(),
                    smoothingParameters[index]);
                destination += term.randomColumns();
            }
            return new PenalizedPredictor(design, penalty,
                rows, columns, parametricColumns);
        }
    }

    public int observations() { return observations; }
    public int columns() { return columns; }
    public int parametricColumns() { return parametricColumns; }
    public double[] design() { return design.clone(); }
    public double[] penaltyDiagonal() { return penaltyDiagonal.clone(); }

    /** Returns the sum of effective degrees of freedom for this block. */
    double effectiveDegreesOfFreedom(
            double[] covariance, int totalColumns, int start) {
        double penaltyTrace = 0.0;
        for (int column = 0; column < columns; column++) {
            penaltyTrace += covariance[(start + column) * totalColumns
                + start + column] * penaltyDiagonal[column];
        }
        return Math.max(0.0, Math.min(columns, columns - penaltyTrace));
    }

    double[] designView() { return design; }
    double[] penaltyDiagonalView() { return penaltyDiagonal; }
}
