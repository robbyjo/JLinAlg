/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.formula;

import java.util.List;
import org.jlinalg.mixed.RandomEffectTerm;

/** Compiler for separate count, structural-zero, and NB2 size formulas. */
public final class ZeroInflatedMixedFormula {
    private ZeroInflatedMixedFormula() { }

    /** Compiles a ZIP pair; independent random terms may occur in either formula. */
    public static CompiledZeroInflatedMixedFormula compilePoisson(
            String countFormula, String zeroFormula, ModelTable table) {
        return compile(countFormula, zeroFormula, null, table);
    }

    /** Compiles an NB2-ZINB triple; size is currently fixed effects only. */
    public static CompiledZeroInflatedMixedFormula compileNegativeBinomial(
            String countFormula, String zeroFormula,
            String sizeFormula, ModelTable table) {
        if (sizeFormula == null)
            throw new IllegalArgumentException("NB2 size formula is required");
        return compile(countFormula, zeroFormula, sizeFormula, table);
    }

    private static CompiledZeroInflatedMixedFormula compile(
            String countFormula, String zeroFormula,
            String sizeFormula, ModelTable table) {
        Predictor count = predictor(countFormula, table);
        Predictor zero = predictor(zeroFormula, table);
        CompiledFormula size = sizeFormula == null
            ? null : Formula.compile(sizeFormula, table);
        if (!java.util.Arrays.equals(count.fixed().responseView(),
                zero.fixed().responseView())
                || size != null && !java.util.Arrays.equals(
                    count.fixed().responseView(), size.responseView()))
            throw new IllegalArgumentException(
                "all zero-inflated formulas must use the same response");
        if (count.random().isEmpty() && zero.random().isEmpty())
            throw new IllegalArgumentException(
                "at least one count or zero random-effect term is required");
        if (count.fixed().weightsView() != null
                || zero.fixed().weightsView() != null
                || size != null && size.weightsView() != null
                || zero.fixed().offsetView() != null
                || size != null && size.offsetView() != null)
            throw new IllegalArgumentException(
                "zero-inflated formulas accept an offset only in the count predictor and do not accept weights");
        return new CompiledZeroInflatedMixedFormula(count.fixed(),
            count.random(), zero.fixed(), zero.random(), size);
    }

    private static Predictor predictor(String formula, ModelTable table) {
        if (formula == null)
            throw new IllegalArgumentException("predictor formula is required");
        if (!formula.contains("|"))
            return new Predictor(Formula.compile(formula, table), List.of());
        CompiledMixedFormula mixed = MixedFormula.compile(formula, table);
        if (!mixed.correlatedRandomEffects().isEmpty())
            throw new IllegalArgumentException(
                "zero-inflated formulas currently require || for multi-coefficient random terms");
        return new Predictor(mixed.fixed(), mixed.randomEffects());
    }

    private record Predictor(
            CompiledFormula fixed, List<RandomEffectTerm> random) { }
}
