/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.penalized;

import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;
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
                || predictors.length != response.length || response.length < 4 || predictors[0] == null
                || !(selectionFraction > 0.0) || !(selectionFraction < 1.0)
                || activeTolerance < 0.0 || !Double.isFinite(activeTolerance)) throw new IllegalArgumentException("selection-aware inputs are invalid");
        int rows = response.length, columns = predictors[0].length, selectionRows = (int) Math.floor(rows * selectionFraction);
        if(selectionRows<2||selectionRows>=rows-1)throw new IllegalArgumentException("each split needs at least two observations");
        for(int row=0;row<rows;row++) {
            if(predictors[row]==null||predictors[row].length!=columns||!Double.isFinite(response[row]))throw new IllegalArgumentException("finite rectangular data required");
            for(double value:predictors[row])if(!Double.isFinite(value))throw new IllegalArgumentException("finite predictors required");
        }
        double[] weights=options.observationWeights();
        if(weights!=null&&weights.length!=rows)throw new IllegalArgumentException("weights must match the full sample");
        var selectionOptions=ElasticNetOptions.builder().alpha(options.alpha()).fitIntercept(options.fitIntercept())
            .standardize(options.standardize()).maximumIterations(options.maximumIterations())
            .parallelism(options.parallelism()).relativeTolerance(options.relativeTolerance());
        if(weights!=null)selectionOptions.observationWeights(Arrays.copyOf(weights,selectionRows));
        if(options.penaltyFactors()!=null)selectionOptions.penaltyFactors(options.penaltyFactors());
        double[] selectionResponse = new double[selectionRows]; double[][] selectionPredictors = new double[selectionRows][columns]; for (int row = 0; row < selectionRows; row++) { selectionResponse[row] = response[row]; selectionPredictors[row] = predictors[row].clone(); }
        PenalizedRegressionResult selected = PenalizedRegression.fit(selectionResponse, selectionPredictors, lambda, selectionOptions.build());
        if(!selected.converged())throw new IllegalStateException("selection fit did not converge");
        List<Integer> active = new ArrayList<>(); for (int column = 0; column < columns; column++) if (Math.abs(selected.coefficients()[column]) > activeTolerance) active.add(column);
        int inferenceRows = rows - selectionRows, inferenceColumns = active.size() + (options.fitIntercept() ? 1 : 0); if (inferenceRows <= inferenceColumns) throw new IllegalArgumentException("inference split has too few observations"); double[][] inferenceDesign = new double[inferenceRows][inferenceColumns]; double[] inferenceResponse = new double[inferenceRows];
        for (int row = 0; row < inferenceRows; row++) { inferenceResponse[row] = response[selectionRows + row]; int destination = 0; if (options.fitIntercept()) inferenceDesign[row][destination++] = 1.0; for (int column : active) inferenceDesign[row][destination++] = predictors[selectionRows + row][column]; }
        int[] activeIndices = new int[active.size()]; for (int index = 0; index < activeIndices.length; index++) activeIndices[index] = active.get(index);
        return new SelectionAwareInferenceResult(selected, activeIndices, selectionRows, inferenceRows,
            Ols.fit(inferenceResponse, inferenceDesign, weights==null?null:Arrays.copyOfRange(weights,selectionRows,rows),
                null, OlsOptions.defaults(), backendPolicy));
    }
}
