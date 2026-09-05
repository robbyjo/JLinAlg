/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.formula;

import java.util.List;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.distributional.SparseZeroInflatedMixedModel;
import org.jlinalg.distributional.ZeroInflatedMixedOptions;
import org.jlinalg.distributional.ZeroInflatedMixedResult;
import org.jlinalg.mixed.RandomEffectTerm;

/** Reusable numerical arrays compiled from zero-inflated mixed formulas. */
public final class CompiledZeroInflatedMixedFormula {
    private final CompiledFormula count;
    private final List<RandomEffectTerm> countRandom;
    private final CompiledFormula zero;
    private final List<RandomEffectTerm> zeroRandom;
    private final CompiledFormula size;

    CompiledZeroInflatedMixedFormula(
            CompiledFormula count, List<RandomEffectTerm> countRandom,
            CompiledFormula zero, List<RandomEffectTerm> zeroRandom,
            CompiledFormula size) {
        this.count = count;
        this.countRandom = List.copyOf(countRandom);
        this.zero = zero;
        this.zeroRandom = List.copyOf(zeroRandom);
        this.size = size;
    }

    public CompiledFormula count() { return count; }
    public CompiledFormula zero() { return zero; }
    public CompiledFormula size() { return size; }
    public List<RandomEffectTerm> countRandomEffects() { return countRandom; }
    public List<RandomEffectTerm> zeroRandomEffects() { return zeroRandom; }

    /** Fits the compiled ZIP model without reparsing or rebuilding designs. */
    public ZeroInflatedMixedResult fitPoisson(
            ZeroInflatedMixedOptions options, BackendPolicy backendPolicy) {
        if (size != null)
            throw new IllegalStateException("compiled formula includes NB2 size");
        return SparseZeroInflatedMixedModel.fitPoisson(count.responseView(),
            count.designView(), count.columns(), zero.designView(),
            zero.columns(), countRandom, null, zeroRandom, null, List.of(),
            count.offsetView(), options, backendPolicy);
    }

    /** Fits the compiled NB2-ZINB model without reparsing or rebuilding designs. */
    public ZeroInflatedMixedResult fitNegativeBinomial(
            ZeroInflatedMixedOptions options, BackendPolicy backendPolicy) {
        if (size == null)
            throw new IllegalStateException("compiled formula has no NB2 size");
        return SparseZeroInflatedMixedModel.fitNegativeBinomial(
            count.responseView(), count.designView(), count.columns(),
            zero.designView(), zero.columns(), size.designView(),
            size.columns(), countRandom, null, zeroRandom, null, List.of(),
            count.offsetView(), options, backendPolicy);
    }
}
