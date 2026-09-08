/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

import java.util.ArrayList;
import java.util.List;
import jdistlib.accelerator.CholeskyFactor;
import jdistlib.accelerator.ComputeBackend;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;

/** GLS meta-regression for a full study-by-study covariance matrix. */
public final class MetaMultilevelRegression {
    private MetaMultilevelRegression() { }

    public static MetaMultilevelRegressionResult fit(List<MetaStudy> studies, double[][] moderators,
                                                     List<String> moderatorNames, double[][] covariance,
                                                     boolean includeIntercept, BackendPolicy backendPolicy) {
        if (studies == null || studies.isEmpty() || moderators == null || moderatorNames == null || covariance == null || backendPolicy == null || moderators.length != studies.size()) throw new IllegalArgumentException("multilevel meta-regression inputs are invalid");
        int rows = studies.size(), moderatorCount = moderators[0] == null ? 0 : moderators[0].length, columns = moderatorCount + (includeIntercept ? 1 : 0);
        if (moderatorCount < 1 || moderatorNames.size() != moderatorCount || rows <= columns) throw new IllegalArgumentException("meta-regression design is invalid");
        double[] y = new double[rows], v = new double[rows * rows], x = new double[rows * columns];
        for (int row = 0; row < rows; row++) { y[row] = studies.get(row).effectSize(); if (!Double.isFinite(y[row])) throw new IllegalArgumentException("effects must be finite"); if (covariance[row] == null || covariance[row].length != rows) throw new IllegalArgumentException("covariance must be square"); for (int column = 0; column < rows; column++) { v[row * rows + column] = covariance[row][column]; if (!Double.isFinite(v[row * rows + column]) || Math.abs(v[row * rows + column] - covariance[column][row]) > 1e-10) throw new IllegalArgumentException("covariance must be finite and symmetric"); } int offset = 0; if (includeIntercept) x[row * columns + offset++] = 1.0; if (moderators[row] == null || moderators[row].length != moderatorCount) throw new IllegalArgumentException("moderator rows must have equal widths"); for (int column = 0; column < moderatorCount; column++) x[row * columns + offset + column] = moderators[row][column]; }
        try (BackendContext context = BackendContext.select(backendPolicy)) {
            ComputeBackend backend = context.backend(); CholeskyFactor factor = backend.dpotrf(v, rows); double[] inverseX = factor.solve(x, columns), inverseY = factor.solve(y); double[] information = new double[columns * columns], right = new double[columns];
            for (int left = 0; left < columns; left++) { for (int rightColumn = 0; rightColumn < columns; rightColumn++) for (int row = 0; row < rows; row++) information[left * columns + rightColumn] += x[row * columns + left] * inverseX[row * columns + rightColumn]; for (int row = 0; row < rows; row++) right[left] += x[row * columns + left] * inverseY[row]; }
            CholeskyFactor informationFactor = backend.dpotrf(information, columns); double[] beta = informationFactor.solve(right), covarianceBeta = informationFactor.solve(identity(columns), columns), fitted = new double[rows], inverseResidual = new double[rows];
            for (int row = 0; row < rows; row++) { for (int column = 0; column < columns; column++) fitted[row] += x[row * columns + column] * beta[column]; }
            double q = 0.0; for (int row = 0; row < rows; row++) inverseResidual[row] = y[row] - fitted[row]; double[] solvedResidual = factor.solve(inverseResidual); for (int row = 0; row < rows; row++) q += inverseResidual[row] * solvedResidual[row];
            List<String> names = new ArrayList<>(columns); if (includeIntercept) names.add("(Intercept)"); names.addAll(moderatorNames);
            return new MetaMultilevelRegressionResult(names, beta, covarianceBeta, q, rows - columns, -0.5 * (factor.logDeterminant() + q), true);
        }
    }
    private static double[] identity(int dimension) { double[] result = new double[dimension * dimension]; for (int index = 0; index < dimension; index++) result[index * dimension + index] = 1.0; return result; }
}
