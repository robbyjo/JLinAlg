/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mixed.LinearMixedModel;
import org.jlinalg.mixed.LinearMixedModelResult;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparseLinearMixedModel;
import org.jlinalg.mixed.SparseLinearMixedModelResult;
import org.jlinalg.mixed.CorrelatedLinearMixedModel;
import org.jlinalg.mixed.CorrelatedLinearMixedModelResult;
import org.jlinalg.mixed.CorrelatedRandomEffectBlock;
import org.jlinalg.mixed.SparseUnstructuredCorrelatedModel;
import org.jlinalg.reml.RemlOptions;

/** A mixed formula compiled once into dense fixed and sparse random designs. */
public final class CompiledMixedFormula {
    private final CompiledFormula fixed;
    private final List<RandomEffectTerm> randomEffects;
    private final List<CorrelatedRandomEffectBlock> correlatedRandomEffects;

    CompiledMixedFormula(
            CompiledFormula fixed, List<RandomEffectTerm> randomEffects,
            List<CorrelatedRandomEffectBlock> correlatedRandomEffects) {
        this.fixed = fixed;
        this.randomEffects = List.copyOf(randomEffects);
        this.correlatedRandomEffects = List.copyOf(correlatedRandomEffects);
    }

    public CompiledFormula fixed() { return fixed; }
    public List<RandomEffectTerm> randomEffects() { return randomEffects; }
    public List<CorrelatedRandomEffectBlock> correlatedRandomEffects() {
        return correlatedRandomEffects;
    }

    /** Fits without reparsing or rebuilding fixed/random design matrices. */
    public LinearMixedModelResult fit(
            RemlOptions options, BackendPolicy backendPolicy) {
        if (!correlatedRandomEffects.isEmpty())
            throw new IllegalArgumentException(
                "formula contains correlated random blocks; use fitCorrelated");
        if (fixed.weightsView() != null || fixed.offsetView() != null)
            return LinearMixedModelResult.fromSparse(fitSparse(options, backendPolicy),
                fixed.responseView(), fixed.designView(), fixed.rows(), fixed.columns(), fixed.offsetView());
        return LinearMixedModel.fit(fixed.responseView(), fixed.designView(),
            fixed.rows(), fixed.columns(), randomEffects, options, backendPolicy);
    }

    /**
     * Fits sparse equations. For correlated terms this exposes latent-coordinate
     * diagnostics; use fitSparseUnstructured().correlatedFit() for original
     * random coefficients and complete block covariance estimates.
     */
    public SparseLinearMixedModelResult fitSparse(
            RemlOptions options, BackendPolicy backendPolicy) {
        if (!correlatedRandomEffects.isEmpty())
            return fitSparseUnstructured(options, backendPolicy).fit();
        List<RandomEffectTerm> terms = new java.util.ArrayList<>();
        double[] roots = weightRoots();
        for (RandomEffectTerm term : randomEffects) {
            double[] values = term.sparseValues();
            int[] starts = term.rowPointers();
            for (int r = 0; r < fixed.rows(); r++)
                for (int k = starts[r]; k < starts[r + 1]; k++) values[k] *= roots[r];
            terms.add(RandomEffectTerm.ofSparseCsr(term.name(), term.observations(), term.coefficients(),
                starts, term.columnIndices(), values, term.coefficientNames()));
        }
        return SparseLinearMixedModel.fit(whitenedResponse(), whitenedFixed(), fixed.rows(),
            fixed.columns(), terms, options, backendPolicy)
            .withObservationScale(fixed.weightsView(), fixed.offsetView());
    }

    /** Fits single-bar correlated blocks with Cholesky covariance parameters. */
    public CorrelatedLinearMixedModelResult fitCorrelated(
            RemlOptions options, BackendPolicy backendPolicy) {
        return fitSparseUnstructured(options, backendPolicy).correlatedFit();
    }

    /** Fits all parsed blocks, estimating one covariance per block, with sparse equations. */
    public SparseUnstructuredCorrelatedModel.Result fitSparseUnstructured(
            RemlOptions options, BackendPolicy backendPolicy) {
        return SparseUnstructuredCorrelatedModel.fit(whitenedResponse(), whitenedFixed(),
            fixed.rows(), fixed.columns(), whitenedBlocks(), options, backendPolicy)
            .withObservationScale(fixed.weightsView(), fixed.offsetView());
    }

    /** ML profile of a fixed coefficient, refitting all covariance and nuisance fixed parameters. */
    public org.jlinalg.mixed.ProfileLikelihoodInterval profileFixedEffect(int coefficient,
            double confidence, double lowerBound, double upperBound,
            RemlOptions options, BackendPolicy backendPolicy) {
        return SparseUnstructuredCorrelatedModel.profileFixedEffect(whitenedResponse(), whitenedFixed(),
            fixed.rows(), fixed.columns(), whitenedBlocks(), coefficient, confidence,
            lowerBound, upperBound, options, backendPolicy);
    }

    private List<CorrelatedRandomEffectBlock> whitenedBlocks() {
        List<CorrelatedRandomEffectBlock> blocks = new java.util.ArrayList<>(correlatedRandomEffects);
        for (RandomEffectTerm term : randomEffects) blocks.add(asScalarBlock(term));
        List<CorrelatedRandomEffectBlock> result = new java.util.ArrayList<>();
        double[] roots = weightRoots();
        for (CorrelatedRandomEffectBlock block : blocks) {
            List<String> groups = new java.util.ArrayList<>();
            for (int index : block.groupIndices()) groups.add(block.groupNames().get(index));
            double[][] design = rows(block.effectDesign(), fixed.rows(), block.effectCount());
            for (int r = 0; r < design.length; r++) for (int c = 0; c < design[r].length; c++)
                design[r][c] *= roots[r];
            result.add(CorrelatedRandomEffectBlock.of(block.name(), groups, block.effectNames(), design));
        }
        return result;
    }

    public org.jlinalg.mixed.ProfileLikelihoodInterval profileResidualSd(double confidence,
            double lowerBound, double upperBound, RemlOptions options, BackendPolicy policy) {
        return SparseUnstructuredCorrelatedModel.profileResidualSd(whitenedResponse(), whitenedFixed(),
            fixed.rows(), fixed.columns(), whitenedBlocks(), confidence, lowerBound, upperBound, options, policy);
    }

    public org.jlinalg.mixed.ProfileLikelihoodInterval profileCorrelation(int block, double confidence,
            RemlOptions options, BackendPolicy policy) {
        return SparseUnstructuredCorrelatedModel.profileCorrelation(whitenedResponse(), whitenedFixed(),
            fixed.rows(), fixed.columns(), whitenedBlocks(), block, confidence, options, policy);
    }

    /** Block order is correlatedRandomEffects followed by scalar randomEffects. */
    public org.jlinalg.mixed.ProfileLikelihoodInterval profileRandomSd(int block, int effect, double confidence,
            double lowerBound, double upperBound, RemlOptions options, BackendPolicy policy) {
        return SparseUnstructuredCorrelatedModel.profileRandomSd(whitenedResponse(), whitenedFixed(),
            fixed.rows(), fixed.columns(), whitenedBlocks(), block, effect, confidence,
            lowerBound, upperBound, options, policy);
    }

    private double[] weightRoots() {
        double[] roots = new double[fixed.rows()], weights = fixed.weightsView();
        for (int r = 0; r < roots.length; r++) {
            double w = weights == null ? 1 : weights[r];
            if (!(w > 0) || !Double.isFinite(w))
                throw new IllegalArgumentException("mixed formula weights must be finite and positive");
            roots[r] = Math.sqrt(w);
        }
        return roots;
    }
    private double[] adjustedResponse() {
        double[] response = fixed.responseView().clone(), offset = fixed.offsetView();
        if (offset != null) for (int r = 0; r < response.length; r++) response[r] -= offset[r];
        return response;
    }
    private double[] whitenedResponse() {
        double[] response = adjustedResponse(), roots = weightRoots();
        for (int r = 0; r < response.length; r++) response[r] *= roots[r];
        return response;
    }
    private double[] whitenedFixed() {
        double[] design = fixed.designView().clone(), roots = weightRoots();
        for (int r = 0; r < fixed.rows(); r++) for (int c = 0; c < fixed.columns(); c++)
            design[r * fixed.columns() + c] *= roots[r];
        return design;
    }

    private static double[][] rows(double[] values, int observations, int columns) {
        double[][] result = new double[observations][columns];
        for (int row = 0; row < observations; row++)
            System.arraycopy(values, row * columns, result[row], 0, columns);
        return result;
    }

    private static CorrelatedRandomEffectBlock asScalarBlock(
            RandomEffectTerm term) {
        if (!term.sparse())
            throw new IllegalArgumentException(
                "dense independent terms cannot be converted to grouped blocks");
        int[] starts = term.rowPointers();
        int[] columns = term.columnIndices();
        double[] values = term.sparseValues();
        java.util.ArrayList<String> groups =
            new java.util.ArrayList<>(term.observations());
        double[][] design = new double[term.observations()][1];
        for (int row = 0; row < term.observations(); row++) {
            if (starts[row + 1] - starts[row] != 1)
                throw new IllegalArgumentException(
                    "formula scalar term must have one group entry per row");
            int entry = starts[row];
            groups.add(term.coefficientNames().get(columns[entry]));
            design[row][0] = values[entry];
        }
        return CorrelatedRandomEffectBlock.of(term.name(), groups,
            List.of(term.name().startsWith("1|")
                ? "(Intercept)" : term.name()), design);
    }
}
