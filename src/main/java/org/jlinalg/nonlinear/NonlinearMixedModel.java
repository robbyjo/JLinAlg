/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.nonlinear;

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
        NonlinearFixedModel.validate(response, initialParameters, mean, options, backendPolicy);
        try (SparseLinearMixedModel.Prepared prepared = SparseLinearMixedModel.prepare(
                response.length, randomEffects, remlOptions, backendPolicy)) {
            return iterate(response, initialParameters, mean, options, prepared);
        }
    }

    public static NonlinearMixedFitResult fitPedigree(double[] response,
            double[] initialParameters, NonlinearMeanFunction mean,
            List<PedigreeRandomEffectTerm> pedigreeEffects,
            List<RandomEffectTerm> ordinaryEffects, RemlOptions remlOptions,
            NonlinearModelOptions options, BackendPolicy backendPolicy) {
        if (pedigreeEffects == null || pedigreeEffects.isEmpty())
            throw new IllegalArgumentException("at least one pedigree random effect is required");
        NonlinearFixedModel.validate(response, initialParameters, mean, options, backendPolicy);
        try (SparseLinearMixedModel.Prepared prepared = SparsePedigreeMixedModel.prepare(
                response.length, pedigreeEffects, ordinaryEffects == null
                    ? List.of() : ordinaryEffects, remlOptions, backendPolicy)) {
            return iterate(response, initialParameters, mean, options, prepared);
        }
    }

    private static NonlinearMixedFitResult iterate(double[] response,
            double[] initial, NonlinearMeanFunction mean,
            NonlinearModelOptions options, SparseLinearMixedModel.Prepared prepared) {
        double[] parameters = initial.clone();
        SparseLinearMixedModelResult fit = null;
        boolean converged = false, needsFinalFit = false;
        int iteration = 0;
        double[] additive = null;
        for (; iteration < options.maximumIterations(); iteration++) {
            NonlinearFixedModel.Evaluation evaluation = NonlinearFixedModel.evaluate(
                response, parameters, mean);
            fit = linearize(evaluation, parameters, prepared);
            needsFinalFit = false;
            additive = additive(fit, evaluation.gradient(), parameters.length);
            double[] currentFitted = actualFitted(parameters, mean, additive);
            double[] step = fit.beta();
            for(int j=0;j<step.length;j++)step[j]-=parameters[j];
            var conditional = new NonlinearFixedModel.Evaluation(currentFitted,
                residuals(response,currentFitted),evaluation.gradient());
            if (NonlinearFixedModel.stationary(conditional,step,response,options)) {
                converged = fit.converged(); iteration++; break;
            }
            // Compare nonlinear means holding the current BLUP fixed. Comparing
            // raw conditional SSE across different variance/BLUP fits is not a
            // valid descent test. The mixed normal equations make this a local
            // Gauss-Newton descent direction for the conditional mean.
            double objective = squaredResidual(response,currentFitted);
            double scale = options.initialStep();
            double[] candidate = parameters.clone();
            boolean accepted = false;
            for(int halving=0;halving<=options.maximumStepHalvings();halving++) {
                for(int j=0;j<step.length;j++)candidate[j]=parameters[j]+scale*step[j];
                try {
                    double trial=squaredResidual(response,actualFitted(candidate,mean,additive));
                    if(Double.isFinite(trial)&&trial<=objective) {accepted=true;break;}
                } catch(IllegalArgumentException invalidTrial) {
                    // A trial may leave the supplied nonlinear mean's domain.
                }
                scale*=.5;
            }
            if(!accepted) {iteration++;break;}
            parameters=candidate;
            needsFinalFit=true;
            if(fit.converged())prepared.warmStart(fit.varianceComponents());
        }
        if (fit == null) throw new IllegalStateException("nonlinear mixed model did not evaluate");
        if(needsFinalFit) {
            var evaluation=NonlinearFixedModel.evaluate(response,parameters,mean);
            fit=linearize(evaluation,parameters,prepared);
            additive=additive(fit,evaluation.gradient(),parameters.length);
            double[] step=fit.beta();for(int j=0;j<step.length;j++)step[j]-=parameters[j];
            double[] currentFitted=actualFitted(parameters,mean,additive);
            converged=fit.converged()&&NonlinearFixedModel.stationary(
                new NonlinearFixedModel.Evaluation(currentFitted,residuals(response,currentFitted),evaluation.gradient()),
                step,response,options);
        }
        double[] fitted = actualFitted(parameters, mean, additive);
        return new NonlinearMixedFitResult(parameters, fitted,
            residuals(response, fitted), squaredResidual(response, fitted),
            iteration, converged && fit.converged(), fit, fit.backend());
    }

    private static SparseLinearMixedModelResult linearize(NonlinearFixedModel.Evaluation evaluation,
            double[] parameters,SparseLinearMixedModel.Prepared prepared) {
        double[] pseudo=evaluation.residual().clone();
        for(int i=0;i<pseudo.length;i++)for(int j=0;j<parameters.length;j++)
            pseudo[i]+=evaluation.gradient()[i*parameters.length+j]*parameters[j];
        return prepared.fit(pseudo,evaluation.gradient(),parameters.length);
    }

    /** Recover Zb from the fitted linearization without densifying any sparse term. */
    private static double[] additive(SparseLinearMixedModelResult fit,double[] gradient,int columns) {
        double[] result=fit.conditionalFittedValues(), beta=fit.beta();
        for(int i=0;i<result.length;i++)for(int j=0;j<columns;j++)result[i]-=gradient[i*columns+j]*beta[j];
        return result;
    }

    private static double[] actualFitted(double[] parameters,
            NonlinearMeanFunction mean, double[] additive) {
        int rows = additive.length;
        double[] result = new double[rows];
        for (int row = 0; row < rows; row++) result[row] = mean.evaluate(parameters, row).value()+additive[row];
        return result;
    }

    private static double squaredResidual(double[] y, double[] fitted) {
        double result = 0.0; for (int i = 0; i < y.length; i++) { double value = y[i] - fitted[i]; result += value * value; } return result;
    }
    private static double[] residuals(double[] y, double[] fitted) { double[] result = y.clone(); for (int i = 0; i < result.length; i++) result[i] -= fitted[i]; return result; }
}
