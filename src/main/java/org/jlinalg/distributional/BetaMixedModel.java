/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.distributional;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.LaplaceTunableFamily;
import org.jlinalg.glmm.GlmmLaplaceResult;
import org.jlinalg.glmm.LaplaceFamilyDerivatives;
import org.jlinalg.glmm.SparseGlmmLaplace;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparsePrecisionMatrix;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;

/**
 * Sparse first-order Laplace beta mixed models with logit-linked conditional
 * means and a profiled constant precision.
 *
 * <p>Random effects enter the mean predictor. Arbitrary sparse
 * coefficient-space precisions make grouped, spatial, and pedigree effects the
 * same computational problem.</p>
 */
public final class BetaMixedModel {
    private BetaMixedModel() { }

    /** Fits grouped random effects with default controls and backend selection. */
    public static BetaMixedModelResult fit(
            double[] response, double[][] fixedEffects,
            List<RandomEffectTerm> randomEffects) {
        if (response == null || fixedEffects == null || fixedEffects.length == 0
                || fixedEffects[0] == null)
            throw new IllegalArgumentException(
                "response and fixed effects are required");
        return fit(response, MatrixOps.rowMajor(fixedEffects, response.length),
            response.length, fixedEffects[0].length, randomEffects, null,
            null, null, BetaMixedModelOptions.defaults(),
            BackendPolicy.PREFERRED);
    }

    /** Fits random effects whose coefficients have independent unit precision. */
    public static BetaMixedModelResult fit(
            double[] response, double[] fixedEffects, int rows, int columns,
            List<RandomEffectTerm> randomEffects,
            double[] priorWeights, double[] offset,
            BetaMixedModelOptions options, BackendPolicy backendPolicy) {
        return fit(response, fixedEffects, rows, columns, randomEffects, null,
            priorWeights, offset, options, backendPolicy);
    }

    /** Fits random effects with caller-supplied sparse precision bases. */
    public static BetaMixedModelResult fit(
            double[] response, double[] fixedEffects, int rows, int columns,
            List<RandomEffectTerm> randomEffects,
            List<SparsePrecisionMatrix> precisionBases,
            double[] priorWeights, double[] offset,
            BetaMixedModelOptions options, BackendPolicy backendPolicy) {
        if (options == null) options = BetaMixedModelOptions.defaults();
        LaplaceBetaFamily family = new LaplaceBetaFamily(
            options.initialPrecision(), options.minimumPrecision(),
            options.maximumPrecision());
        GlmmLaplaceResult result = SparseGlmmLaplace.fitWithPrecision(response,
            fixedEffects, rows, columns, family, randomEffects, precisionBases,
            priorWeights, offset, options.laplace(), backendPolicy);
        return new BetaMixedModelResult(result, family.precision());
    }

    /** Convenience overload for one additive genetic pedigree effect. */
    public static BetaMixedModelResult fitPedigree(
            double[] response, double[] fixedEffects, int rows, int columns,
            PedigreeRandomEffectTerm pedigree,
            double[] priorWeights, double[] offset,
            BetaMixedModelOptions options, BackendPolicy backendPolicy) {
        if (pedigree == null)
            throw new IllegalArgumentException("pedigree term is required");
        return fit(response, fixedEffects, rows, columns,
            List.of(pedigree.randomEffect()), List.of(pedigree.precision()),
            priorWeights, offset, options, backendPolicy);
    }

    private static final class LaplaceBetaFamily
            implements LaplaceTunableFamily, LaplaceFamilyDerivatives {
        private static final double MINIMUM_MEAN = 1e-14;
        private final double minimumLogPrecision;
        private final double maximumLogPrecision;
        private double logPrecision;

        LaplaceBetaFamily(double initialPrecision, double minimumPrecision,
                double maximumPrecision) {
            this.logPrecision = Math.log(initialPrecision);
            this.minimumLogPrecision = Math.log(minimumPrecision);
            this.maximumLogPrecision = Math.log(maximumPrecision);
        }

        double precision() { return Math.exp(logPrecision); }
        @Override public String name() { return "beta(logit)"; }
        @Override public void validateResponse(double response, double weight) {
            if (!(response > 0.0 && response < 1.0)
                    || !Double.isFinite(response))
                throw new IllegalArgumentException(
                    "beta responses must lie strictly between zero and one");
        }
        @Override public double initialMean(double response) { return response; }
        @Override public double link(double mean) {
            return Math.log(mean) - Math.log1p(-mean);
        }
        @Override public double inverseLink(double predictor) {
            double value = predictor >= 0.0
                ? 1.0 / (1.0 + Math.exp(-predictor))
                : Math.exp(predictor) / (1.0 + Math.exp(predictor));
            return Math.max(MINIMUM_MEAN,
                Math.min(1.0 - MINIMUM_MEAN, value));
        }
        @Override public double meanDerivative(double predictor) {
            double mean = inverseLink(predictor);
            return mean * (1.0 - mean);
        }
        @Override public double variance(double mean) {
            return mean * (1.0 - mean) / (1.0 + precision());
        }
        @Override public double unitDeviance(double response, double mean) {
            return -2.0 * logLikelihood(response, mean, 1.0, 1.0);
        }
        @Override public double logLikelihood(double response, double mean,
                double priorWeight, double dispersion) {
            double phi = precision();
            return priorWeight * betaLogLikelihood(response, mean, phi);
        }
        @Override public boolean fixedDispersion() { return true; }

        @Override public double[] laplaceParameters() {
            return new double[] {logPrecision};
        }
        @Override public void setLaplaceParameters(double[] parameters) {
            if (parameters == null || parameters.length != 1)
                throw new IllegalArgumentException(
                    "beta family requires one log-precision parameter");
            logPrecision = parameters[0];
        }
        @Override public double minimumLaplaceParameter(int index) {
            if (index != 0) throw new IndexOutOfBoundsException(index);
            return minimumLogPrecision;
        }
        @Override public double maximumLaplaceParameter(int index) {
            if (index != 0) throw new IndexOutOfBoundsException(index);
            return maximumLogPrecision;
        }

        @Override public double workingWeight(double response,
                double predictor, double mean, double priorWeight) {
            double phi = precision();
            double derivative = mean * (1.0 - mean);
            double information = phi * phi
                * (SpecialFunctions.trigamma(mean * phi)
                    + SpecialFunctions.trigamma((1.0 - mean) * phi))
                * derivative * derivative;
            return priorWeight * information;
        }

        @Override public double linearPredictorScore(double response, double predictor, double priorWeight) {
            double mean = inverseLink(predictor), phi = precision();
            double transformedResidual = Math.log(response) - Math.log1p(-response)
                - SpecialFunctions.digamma(mean * phi) + SpecialFunctions.digamma((1 - mean) * phi);
            return priorWeight * phi * mean * (1 - mean) * transformedResidual;
        }

        @Override public double linearPredictorInformation(double response, double predictor, double priorWeight) {
            double mean = inverseLink(predictor);
            // The second derivative of the logit mean contributes a score term.
            // Unlike Fisher information this observed curvature can be negative.
            return workingWeight(response, predictor, mean, priorWeight)
                - (1 - 2 * mean) * linearPredictorScore(response, predictor, priorWeight);
        }

        @Override public double workingResponse(double response,
                double predictor, double mean, double priorWeight,
                double offset) {
            double phi = precision();
            double derivative = mean * (1.0 - mean);
            double transformedResidual = Math.log(response)
                - Math.log1p(-response)
                - SpecialFunctions.digamma(mean * phi)
                + SpecialFunctions.digamma((1.0 - mean) * phi);
            double score = phi * transformedResidual * derivative;
            double information = phi * phi
                * (SpecialFunctions.trigamma(mean * phi)
                    + SpecialFunctions.trigamma((1.0 - mean) * phi))
                * derivative * derivative;
            return predictor - offset + score / information;
        }

        private static double betaLogLikelihood(
                double response, double mean, double phi) {
            double alpha = mean * phi;
            double beta = (1.0 - mean) * phi;
            return jdistlib.Beta.density(response, alpha, beta, true);
        }
    }
}
