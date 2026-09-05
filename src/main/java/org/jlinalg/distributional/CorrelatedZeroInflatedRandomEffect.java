/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.distributional;

import java.util.Objects;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.mixed.SparsePrecisionMatrix;
import org.jlinalg.pedigree.PedigreeRandomEffectTerm;

/**
 * Matched count and structural-zero random coefficients sharing a two-process
 * covariance matrix and one sparse coefficient precision.
 *
 * <p>If the coefficient precision is {@code A^-1}, the joint precision is
 * {@code G^-1 (x) A^-1}. The fitted {@code G} contains a count variance, a
 * structural-zero variance, and a guarded correlation.</p>
 */
public record CorrelatedZeroInflatedRandomEffect(
        String name,
        RandomEffectTerm countTerm,
        RandomEffectTerm zeroTerm,
        SparsePrecisionMatrix coefficientPrecision) {

    public CorrelatedZeroInflatedRandomEffect {
        if (name == null || name.isBlank())
            throw new IllegalArgumentException("correlated effect name is required");
        Objects.requireNonNull(countTerm, "countTerm");
        Objects.requireNonNull(zeroTerm, "zeroTerm");
        Objects.requireNonNull(coefficientPrecision, "coefficientPrecision");
        if (countTerm.observations() != zeroTerm.observations()
                || countTerm.coefficients() != zeroTerm.coefficients()
                || countTerm.coefficients() != coefficientPrecision.dimension()
                || !countTerm.coefficientNames().equals(
                    zeroTerm.coefficientNames()))
            throw new IllegalArgumentException(
                "correlated count/zero terms and precision must have matching coefficients");
    }

    /** Creates matched count/zero effects from one pedigree incidence term. */
    public static CorrelatedZeroInflatedRandomEffect pedigree(
            String name, PedigreeRandomEffectTerm pedigree) {
        Objects.requireNonNull(pedigree, "pedigree");
        RandomEffectTerm source = pedigree.randomEffect();
        RandomEffectTerm count = copy(name + ":count", source);
        RandomEffectTerm zero = copy(name + ":zero", source);
        return new CorrelatedZeroInflatedRandomEffect(
            name, count, zero, pedigree.precision());
    }

    private static RandomEffectTerm copy(String name, RandomEffectTerm source) {
        if (source.sparse())
            return RandomEffectTerm.ofSparseCsr(name, source.observations(),
                source.coefficients(), source.rowPointers(),
                source.columnIndices(), source.sparseValues(),
                source.coefficientNames());
        return RandomEffectTerm.of(name, source.design(), source.observations(),
            source.coefficients(), source.coefficientNames());
    }
}
