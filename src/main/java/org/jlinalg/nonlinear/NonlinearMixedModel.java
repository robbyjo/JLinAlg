/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.nonlinear;

import java.util.ArrayList;
import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparseLinearMixedModel;
import org.jlinalg.mixed.SparseLinearMixedModelResult;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.jlinalg.pedigree.SparsePedigreeMixedModel;
import org.jlinalg.reml.RemlOptions;

/** Gauss-Newton nonlinear mixed-effects fitting on prepared sparse equations. */
public final class NonlinearMixedModel {
    private NonlinearMixedModel() { }

    public static NonlinearMixedFitResult fit(double[] response,
            double[] initialParameters, NonlinearMeanFunction mean,
            List<RandomEffectTerm> randomEffects, RemlOptions remlOptions,
            NonlinearModelOptions options, BackendPolicy backendPolicy) {
        if (randomEffects == null || randomEffects.isEmpty())
            throw new IllegalArgumentException("at least one ordinary random effect is required");
        try (SparseLinearMixedModel.Prepared prepared = SparseLinearMixedModel.prepare(
                response.length, randomEffects, remlOptions, backendPolicy)) {
            return iterate(response, initialParameters, mean, options, backendPolicy,
                prepared, null, randomEffects, List.of());
        }
    }

    public static NonlinearMixedFitResult fitPedigree(double[] response,
            double[] initialParameters, NonlinearMeanFunction mean,
            List<PedigreeRandomEffectTerm> pedigreeEffects,
            List<RandomEffectTerm> ordinaryEffects, RemlOptions remlOptions,
            NonlinearModelOptions options, BackendPolicy backendPolicy) {
        if (pedigreeEffects == null || pedigreeEffects.isEmpty())
            throw new IllegalArgumentException("at least one pedigree random effect is required");
        try (SparseLinearMixedModel.Prepared prepared = SparsePedigreeMixedModel.prepare(
                response.length, pedigreeEffects, ordinaryEffects == null
                    ? List.of() : ordinaryEffects, remlOptions, backendPolicy)) {
            List<RandomEffectTerm> terms = new ArrayList<>();
            for (PedigreeRandomEffectTerm term : pedigreeEffects) terms.add(term.randomEffect());
            if (ordinaryEffects != null) terms.addAll(ordinaryEffects);
            return iterate(response, initialParameters, mean, options, backendPolicy,
                prepared, null, terms, pedigreeEffects);
        }
    }

    private static NonlinearMixedFitResult iterate(double[] response,
            double[] initial, NonlinearMeanFunction mean,
            NonlinearModelOptions options, BackendPolicy backendPolicy,
            SparseLinearMixedModel.Prepared prepared,
            Object ignored, List<RandomEffectTerm> terms,
            List<PedigreeRandomEffectTerm> pedigreeTerms) {
        if (response == null || initial == null || mean == null || options == null
                || response.length < 2 || initial.length == 0)
            throw new IllegalArgumentException("response, parameters, mean, and options are required");
        double[] parameters = initial.clone();
        SparseLinearMixedModelResult fit = null;
        double objective = Double.POSITIVE_INFINITY; boolean converged = false;
        int iteration = 0;
        for (; iteration < options.maximumIterations(); iteration++) {
            NonlinearFixedModel.Evaluation evaluation = NonlinearFixedModel.evaluate(
                response, parameters, mean);
            double[] pseudo = response.clone();
            for (int row = 0; row < response.length; row++) {
                double linear = 0.0;
                for (int column = 0; column < parameters.length; column++)
                    linear += evaluation.gradient()[row * parameters.length + column] * parameters[column];
                pseudo[row] += linear - evaluation.fitted()[row];
            }
            fit = prepared.fit(pseudo, evaluation.gradient(), parameters.length);
            double[] candidate = fit.beta();
            double[] actualFitted = actualFitted(candidate, mean, fit, terms);
            double candidateObjective = squaredResidual(response, actualFitted);
            double previousObjective = objective;
            if (candidateObjective > objective && objective < Double.POSITIVE_INFINITY) {
                double[] damped = parameters.clone();
                for (int column = 0; column < damped.length; column++)
                    damped[column] = 0.5 * (damped[column] + candidate[column]);
                candidate = damped;
                actualFitted = actualFitted(candidate, mean, fit, terms);
                candidateObjective = squaredResidual(response, actualFitted);
            }
            double change = distance(parameters, candidate);
            if (candidateObjective <= objective) { parameters = candidate; objective = candidateObjective; }
            double improvement = previousObjective - candidateObjective;
            if (change <= options.parameterTolerance() * (1.0 + distance(new double[parameters.length], parameters))
                    || (Double.isFinite(previousObjective)
                        && improvement <= options.objectiveTolerance() * (1.0 + objective))) {
                converged = true; break;
            }
        }
        if (fit == null) throw new IllegalStateException("nonlinear mixed model did not evaluate");
        double[] fitted = actualFitted(parameters, mean, fit, terms);
        return new NonlinearMixedFitResult(parameters, fitted,
            residuals(response, fitted), squaredResidual(response, fitted),
            iteration + 1, converged, fit, fit.backend());
    }

    private static double[] actualFitted(double[] parameters,
            NonlinearMeanFunction mean, SparseLinearMixedModelResult fit,
            List<RandomEffectTerm> terms) {
        int rows = fit.conditionalFittedValues().length;
        double[] result = new double[rows];
        for (int row = 0; row < rows; row++) result[row] = mean.evaluate(parameters, row).value();
        for (RandomEffectTerm term : terms) {
            double[] coefficients = fit.randomEffects(term.name()).estimates();
            double[] design = term.design();
            for (int row = 0; row < rows; row++) for (int column = 0; column < term.coefficients(); column++)
                result[row] += design[row * term.coefficients() + column] * coefficients[column];
        }
        return result;
    }

    private static double squaredResidual(double[] y, double[] fitted) {
        double result = 0.0; for (int i = 0; i < y.length; i++) { double value = y[i] - fitted[i]; result += value * value; } return result;
    }
    private static double[] residuals(double[] y, double[] fitted) { double[] result = y.clone(); for (int i = 0; i < result.length; i++) result[i] -= fitted[i]; return result; }
    private static double distance(double[] left, double[] right) { double result = 0.0; for (int i = 0; i < left.length; i++) { double value = left[i] - right[i]; result += value * value; } return Math.sqrt(result); }
}
