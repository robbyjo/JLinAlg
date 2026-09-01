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
        if (fixed.weightsView() != null || fixed.offsetView() != null) {
            throw new IllegalArgumentException(
                "mixed formula weights and offsets require a weighted LMM fit");
        }
        return LinearMixedModel.fit(fixed.responseView(), fixed.designView(),
            fixed.rows(), fixed.columns(), randomEffects, options, backendPolicy);
    }

    /** Fits through sparse mixed-model equations without dense n-by-n covariance. */
    public SparseLinearMixedModelResult fitSparse(
            RemlOptions options, BackendPolicy backendPolicy) {
        if (!correlatedRandomEffects.isEmpty())
            throw new IllegalArgumentException(
                "correlated random blocks currently use fitCorrelated");
        if (fixed.weightsView() != null || fixed.offsetView() != null)
            throw new IllegalArgumentException(
                "sparse mixed formulas do not yet accept weights or offsets");
        return SparseLinearMixedModel.fit(fixed.responseView(),
            fixed.designView(), fixed.rows(), fixed.columns(), randomEffects,
            options, backendPolicy);
    }

    /** Fits single-bar correlated blocks with Cholesky covariance parameters. */
    public CorrelatedLinearMixedModelResult fitCorrelated(
            RemlOptions options, BackendPolicy backendPolicy) {
        if (fixed.weightsView() != null || fixed.offsetView() != null)
            throw new IllegalArgumentException(
                "correlated mixed formulas do not yet accept weights or offsets");
        List<CorrelatedRandomEffectBlock> blocks =
            new java.util.ArrayList<>(correlatedRandomEffects);
        for (RandomEffectTerm term : randomEffects)
            blocks.add(asScalarBlock(term));
        return CorrelatedLinearMixedModel.fit(fixed.responseView(),
            fixed.designView(), fixed.rows(), fixed.columns(), blocks,
            options, backendPolicy);
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
