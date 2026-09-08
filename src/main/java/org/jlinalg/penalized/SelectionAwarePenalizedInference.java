/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.penalized;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.ols.Ols;
import org.jlinalg.ols.OlsOptions;

/** Honest selection-aware inference through a deterministic held-out split. */
public final class SelectionAwarePenalizedInference {
    private SelectionAwarePenalizedInference() { }

    public static SelectionAwareInferenceResult fit(double[] response, double[][] predictors,
                                                    double lambda, ElasticNetOptions options,
                                                    double selectionFraction, double activeTolerance,
                                                    BackendPolicy backendPolicy) {
        if (response == null || predictors == null || options == null || backendPolicy == null
                || predictors.length != response.length || !(selectionFraction > 0.0) || !(selectionFraction < 1.0)
                || activeTolerance < 0.0 || !Double.isFinite(activeTolerance)) throw new IllegalArgumentException("selection-aware inputs are invalid");
        int rows = response.length, columns = predictors[0].length, selectionRows = Math.max(columns + 2, (int) Math.floor(rows * selectionFraction)); if (selectionRows >= rows) throw new IllegalArgumentException("selection split leaves no inference rows");
        double[] selectionResponse = new double[selectionRows]; double[][] selectionPredictors = new double[selectionRows][columns]; for (int row = 0; row < selectionRows; row++) { selectionResponse[row] = response[row]; selectionPredictors[row] = predictors[row].clone(); }
        PenalizedRegressionResult selected = PenalizedRegression.fit(selectionResponse, selectionPredictors, lambda, options);
        List<Integer> active = new ArrayList<>(); for (int column = 0; column < columns; column++) if (Math.abs(selected.coefficients()[column]) > activeTolerance) active.add(column);
        int inferenceRows = rows - selectionRows, inferenceColumns = active.size() + (options.fitIntercept() ? 1 : 0); if (inferenceRows <= inferenceColumns) throw new IllegalArgumentException("inference split has too few observations"); double[][] inferenceDesign = new double[inferenceRows][inferenceColumns]; double[] inferenceResponse = new double[inferenceRows];
        for (int row = 0; row < inferenceRows; row++) { inferenceResponse[row] = response[selectionRows + row]; int destination = 0; if (options.fitIntercept()) inferenceDesign[row][destination++] = 1.0; for (int column : active) inferenceDesign[row][destination++] = predictors[selectionRows + row][column]; }
        int[] activeIndices = new int[active.size()]; for (int index = 0; index < activeIndices.length; index++) activeIndices[index] = active.get(index);
        return new SelectionAwareInferenceResult(selected, activeIndices, selectionRows, inferenceRows, Ols.fit(inferenceResponse, inferenceDesign, OlsOptions.defaults(), backendPolicy));
    }
}
