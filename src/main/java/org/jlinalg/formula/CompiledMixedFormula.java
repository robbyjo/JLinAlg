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
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;
import org.jlinalg.pedigree.SparsePedigreeMixedModel;
import org.jlinalg.reml.RemlOptions;

/** A mixed formula compiled once into dense fixed and sparse random designs. */
public final class CompiledMixedFormula {
    private final CompiledFormula fixed;
    private final List<RandomEffectTerm> randomEffects;
    private final List<CorrelatedRandomEffectBlock> correlatedRandomEffects;
    private final List<PedigreeRandomEffectTerm> pedigreeRandomEffects;
    private final int[] retainedRows;
    private final int originalRows;

    CompiledMixedFormula(
            CompiledFormula fixed, List<RandomEffectTerm> randomEffects,
            List<CorrelatedRandomEffectBlock> correlatedRandomEffects,
            List<PedigreeRandomEffectTerm> pedigreeRandomEffects,
            int[] retainedRows, int originalRows) {
        this.fixed = fixed;
        this.randomEffects = List.copyOf(randomEffects);
        this.correlatedRandomEffects = List.copyOf(correlatedRandomEffects);
        this.pedigreeRandomEffects = List.copyOf(pedigreeRandomEffects);
        this.retainedRows = retainedRows.clone();
        this.originalRows = originalRows;
    }

    public CompiledFormula fixed() { return fixed; }
    public List<RandomEffectTerm> randomEffects() { return randomEffects; }
    public List<CorrelatedRandomEffectBlock> correlatedRandomEffects() {
        return correlatedRandomEffects;
    }
    public List<PedigreeRandomEffectTerm> pedigreeRandomEffects() {
        return pedigreeRandomEffects;
    }
    /** Original table rows retained by the joint complete-case selection. */
    public int[] retainedRows() { return retainedRows.clone(); }
    public int originalRows() { return originalRows; }
    public int omittedRows() { return originalRows - retainedRows.length; }

    /** Fits without reparsing or rebuilding fixed/random design matrices. */
    public LinearMixedModelResult fit(
            RemlOptions options, BackendPolicy backendPolicy) {
        if (!correlatedRandomEffects.isEmpty())
            throw new IllegalArgumentException(
                "formula contains correlated random blocks; use fitCorrelated");
        if (fixed.weightsView() != null || fixed.offsetView() != null
                || !pedigreeRandomEffects.isEmpty())
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
        List<RandomEffectTerm> terms = whitenedRandomEffects(randomEffects);
        SparseLinearMixedModelResult result = pedigreeRandomEffects.isEmpty()
            ? SparseLinearMixedModel.fit(whitenedResponse(), whitenedFixed(),
                fixed.rows(), fixed.columns(), terms, options, backendPolicy)
            : SparsePedigreeMixedModel.fit(whitenedResponse(), whitenedFixed(),
                fixed.rows(), fixed.columns(), whitenedPedigreeEffects(),
                terms, options, backendPolicy);
        return result.withObservationScale(
            fixed.weightsView(), fixed.offsetView());
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
            fixed.rows(), fixed.columns(), whitenedBlocks(),
            whitenedPedigreeEffects(), options, backendPolicy)
            .withObservationScale(fixed.weightsView(), fixed.offsetView());
    }

    /** ML profile of a fixed coefficient, refitting all covariance and nuisance fixed parameters. */
    public org.jlinalg.mixed.ProfileLikelihoodInterval profileFixedEffect(int coefficient,
            double confidence, double lowerBound, double upperBound,
            RemlOptions options, BackendPolicy backendPolicy) {
        requireNoPedigreeCovarianceOptimization();
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
        requireNoPedigreeCovarianceOptimization();
        return SparseUnstructuredCorrelatedModel.profileResidualSd(whitenedResponse(), whitenedFixed(),
            fixed.rows(), fixed.columns(), whitenedBlocks(), confidence, lowerBound, upperBound, options, policy);
    }

    public org.jlinalg.mixed.ProfileLikelihoodInterval profileCorrelation(int block, double confidence,
            RemlOptions options, BackendPolicy policy) {
        requireNoPedigreeCovarianceOptimization();
        return SparseUnstructuredCorrelatedModel.profileCorrelation(whitenedResponse(), whitenedFixed(),
            fixed.rows(), fixed.columns(), whitenedBlocks(), block, confidence, options, policy);
    }

    /** Profiles any pairwise correlation in a correlated random-effect block. */
    public org.jlinalg.mixed.ProfileLikelihoodInterval profileCorrelation(
            int block, int firstEffect, int secondEffect, double confidence,
            RemlOptions options, BackendPolicy policy) {
        requireNoPedigreeCovarianceOptimization();
        return SparseUnstructuredCorrelatedModel.profileCorrelation(
            whitenedResponse(), whitenedFixed(), fixed.rows(), fixed.columns(),
            whitenedBlocks(), block, firstEffect, secondEffect,
            confidence, options, policy);
    }

    /** Block order is correlatedRandomEffects followed by scalar randomEffects. */
    public org.jlinalg.mixed.ProfileLikelihoodInterval profileRandomSd(int block, int effect, double confidence,
            double lowerBound, double upperBound, RemlOptions options, BackendPolicy policy) {
        requireNoPedigreeCovarianceOptimization();
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

    private List<RandomEffectTerm> whitenedRandomEffects(
            List<RandomEffectTerm> source) {
        List<RandomEffectTerm> result = new java.util.ArrayList<>();
        double[] roots = weightRoots();
        for (RandomEffectTerm term : source) {
            if (!term.sparse()) {
                throw new IllegalArgumentException(
                    "formula random effects must use sparse designs");
            }
            double[] values = term.sparseValues();
            int[] starts = term.rowPointers();
            for (int row = 0; row < fixed.rows(); row++) {
                for (int entry = starts[row]; entry < starts[row + 1]; entry++) {
                    values[entry] *= roots[row];
                }
            }
            result.add(RandomEffectTerm.ofSparseCsr(
                term.name(), term.observations(), term.coefficients(),
                starts, term.columnIndices(), values, term.coefficientNames()));
        }
        return result;
    }

    private List<PedigreeRandomEffectTerm> whitenedPedigreeEffects() {
        List<PedigreeRandomEffectTerm> result = new java.util.ArrayList<>();
        for (PedigreeRandomEffectTerm term : pedigreeRandomEffects) {
            RandomEffectTerm whitened = whitenedRandomEffects(
                List.of(term.randomEffect())).get(0);
            result.add(new PedigreeRandomEffectTerm(
                whitened, term.precision()));
        }
        return result;
    }

    private void requireNoPedigreeCovarianceOptimization() {
        if (!pedigreeRandomEffects.isEmpty()) {
            throw new UnsupportedOperationException(
                "joint unstructured-plus-pedigree covariance optimization is not implemented");
        }
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
