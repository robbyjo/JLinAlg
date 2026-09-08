/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.association;

import java.util.List;
import org.jlinalg.glm.Glm;
import org.jlinalg.glm.GlmFamily;
import org.jlinalg.glm.GlmOptions;
import org.jlinalg.glmm.GlmmPql;
import org.jlinalg.glmm.GlmmPqlOptions;
import org.jlinalg.mixed.LinearMixedModel;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparseLinearMixedModel;
import org.jlinalg.mixed.CorrelatedLinearMixedModel;
import org.jlinalg.mixed.CorrelatedRandomEffectBlock;
import org.jlinalg.ols.Ols;
import org.jlinalg.ols.OlsOptions;
import org.jlinalg.penalized.ElasticNetOptions;
import org.jlinalg.penalized.PenalizedRegressionInference;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeGlmmPql;
import org.jlinalg.pedigree.PedigreeReml;
import org.jlinalg.pedigree.SparsePedigreeReml;
import org.jlinalg.reml.Reml;
import org.jlinalg.reml.RemlOptions;
import org.jlinalg.reml.VarianceComponent;
import org.jlinalg.timeseries.ArimaErrorLinearMixedModel;
import org.jlinalg.timeseries.ArimaErrorLmmOptions;
import org.jlinalg.timeseries.ArimaOrder;

/** Thread-safe built-in adapters for the parallel association engine. */
public final class AssociationModels {
    private AssociationModels() { }

    public static AssociationFitter ols(OlsOptions options) {
        if (options == null) throw new IllegalArgumentException("options are required");
        return (response, design, rows, columns, backend) ->
            Ols.fit(response, design, rows, columns, options, backend)
                .associationStatistics();
    }

    public static AssociationFitter weightedOls(
            double[] weights, double[] offset, OlsOptions options) {
        double[] retainedWeights = weights == null ? null : weights.clone();
        double[] retainedOffset = offset == null ? null : offset.clone();
        return (response, design, rows, columns, backend) ->
            Ols.fit(response, design, rows, columns,
                retainedWeights, retainedOffset, options, backend)
                .associationStatistics();
    }

    public static AssociationFitter glm(
            GlmFamily family, GlmOptions options) {
        return glm(family, null, null, options);
    }

    public static AssociationFitter glm(
            GlmFamily family, double[] weights, double[] offset,
            GlmOptions options) {
        if (family == null || options == null)
            throw new IllegalArgumentException("family and options are required");
        double[] retainedWeights = weights == null ? null : weights.clone();
        double[] retainedOffset = offset == null ? null : offset.clone();
        return (response, design, rows, columns, backend) -> {
            var fit = Glm.fit(response, design, rows, columns, family,
                retainedWeights, retainedOffset, options, backend);
            requireConverged(fit.converged());
            return fit.associationStatistics();
        };
    }

    public static AssociationFitter reml(
            List<VarianceComponent> components, RemlOptions options) {
        List<VarianceComponent> retained = List.copyOf(components);
        return (response, design, rows, columns, backend) -> {
            var fit = Reml.fit(response, design, rows, columns, retained, options, backend);
            requireConverged(fit.converged());
            return fit.associationStatistics();
        };
    }

    public static AssociationFitter linearMixedModel(
            List<RandomEffectTerm> randomEffects, RemlOptions options) {
        List<RandomEffectTerm> retained = List.copyOf(randomEffects);
        return (response, design, rows, columns, backend) -> {
            var fit = LinearMixedModel.fit(response, design, rows, columns, retained, options, backend);
            requireConverged(fit.reml().converged());
            return fit.associationStatistics();
        };
    }

    /** Sparse-equation exact-refit LMM association adapter. */
    public static AssociationFitter sparseLinearMixedModel(
            List<RandomEffectTerm> randomEffects, RemlOptions options) {
        List<RandomEffectTerm> retained = List.copyOf(randomEffects);
        return (response, design, rows, columns, backend) -> {
            var fit = SparseLinearMixedModel.fit(response, design, rows, columns, retained, options, backend);
            requireConverged(fit.converged());
            return fit.associationStatistics();
        };
    }

    public static AssociationFitter correlatedLinearMixedModel(
            List<CorrelatedRandomEffectBlock> randomEffects,
            RemlOptions options) {
        List<CorrelatedRandomEffectBlock> retained = List.copyOf(randomEffects);
        return (response, design, rows, columns, backend) -> {
            var fit = CorrelatedLinearMixedModel.fit(response, design, rows, columns, retained, options, backend);
            requireConverged(fit.converged());
            return fit.associationStatistics();
        };
    }

    public static AssociationFitter linearMixedModel(
            List<RandomEffectTerm> randomEffects,
            double[] residualCorrelation,
            RemlOptions options) {
        List<RandomEffectTerm> retained = List.copyOf(randomEffects);
        if (residualCorrelation == null || options == null)
            throw new IllegalArgumentException(
                "residual correlation and options are required");
        double[] correlation = residualCorrelation.clone();
        return (response, design, rows, columns, backend) -> {
            var fit = LinearMixedModel.fitWithResidualCorrelation(
                response, design, rows, columns, retained, correlation,
                options, backend);
            requireConverged(fit.reml().converged());
            return fit.associationStatistics();
        };
    }

    public static AssociationFitter pedigreeReml(
            List<String> observationIndividualIds,
            Pedigree pedigree,
            RemlOptions options) {
        List<String> ids = List.copyOf(observationIndividualIds);
        return (response, design, rows, columns, backend) -> {
            var fit = PedigreeReml.fit(response, design, rows, columns, ids, pedigree, options, backend);
            requireConverged(fit.reml().converged());
            return fit.associationStatistics();
        };
    }

    /** Sparse A-inverse pedigree REML association adapter. */
    public static AssociationFitter sparsePedigreeReml(
            List<String> observationIndividualIds,
            Pedigree pedigree,
            RemlOptions options) {
        List<String> ids = List.copyOf(observationIndividualIds);
        return (response, design, rows, columns, backend) -> {
            var fit = SparsePedigreeReml.fit(response, design, rows, columns, ids, pedigree, options, backend);
            requireConverged(fit.mixedModel().converged());
            return fit.associationStatistics();
        };
    }

    public static AssociationFitter glmmPql(
            GlmFamily family,
            List<VarianceComponent> randomComponents,
            double[] weights,
            double[] offset,
            GlmmPqlOptions options) {
        List<VarianceComponent> retained = List.copyOf(randomComponents);
        double[] retainedWeights = weights == null ? null : weights.clone();
        double[] retainedOffset = offset == null ? null : offset.clone();
        return (response, design, rows, columns, backend) -> {
            var fit = GlmmPql.fit(response, design, rows, columns,
                family, retained, retainedWeights, retainedOffset,
                options, backend);
            requireConverged(fit.converged());
            return fit.associationStatistics();
        };
    }

    public static AssociationFitter pedigreeGlmmPql(
            GlmFamily family,
            List<String> observationIndividualIds,
            Pedigree pedigree,
            double[] weights,
            double[] offset,
            GlmmPqlOptions options) {
        List<String> ids = List.copyOf(observationIndividualIds);
        double[] retainedWeights = weights == null ? null : weights.clone();
        double[] retainedOffset = offset == null ? null : offset.clone();
        return (response, design, rows, columns, backend) -> {
            var fit = PedigreeGlmmPql.fit(response, design, rows, columns,
                family, ids, pedigree, retainedWeights, retainedOffset,
                options, backend);
            requireConverged(fit.glmm().converged());
            return fit.associationStatistics();
        };
    }

    /** Ridge association with effective-DF model-based inference. */
    public static AssociationFitter ridge(
            double lambda, ElasticNetOptions options) {
        if (!Double.isFinite(lambda) || lambda < 0.0 || options == null)
            throw new IllegalArgumentException(
                "lambda must be finite and nonnegative and options are required");
        return (response, design, rows, columns, backend) ->
            PenalizedRegressionInference.ridge(response,
                matrix(design, rows, columns), lambda, options, backend)
                .associationStatistics();
    }

    /** Profile-REML LMM with ARIMA-correlated errors. */
    public static AssociationFitter arimaErrorLinearMixedModel(
            List<RandomEffectTerm> randomEffects,
            ArimaOrder errorOrder,
            ArimaErrorLmmOptions options) {
        List<RandomEffectTerm> retained = List.copyOf(randomEffects);
        if (errorOrder == null || options == null)
            throw new IllegalArgumentException("order and options are required");
        return (response, design, rows, columns, backend) -> {
            var fit = ArimaErrorLinearMixedModel.fit(response,
                matrix(design, rows, columns), retained, errorOrder,
                options, backend);
            requireConverged(fit.converged());
            return fit.associationStatistics();
        };
    }

    private static void requireConverged(boolean converged) {
        if (!converged) throw new IllegalArgumentException("association model did not converge");
    }

    private static double[][] matrix(
            double[] values, int rows, int columns) {
        double[][] result = new double[rows][columns];
        for (int row = 0; row < rows; row++)
            System.arraycopy(values, row * columns,
                result[row], 0, columns);
        return result;
    }
}
