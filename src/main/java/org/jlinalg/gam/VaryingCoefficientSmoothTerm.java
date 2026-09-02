/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import org.jlinalg.internal.MatrixOps;

/** Factor-by or numeric varying-coefficient P-spline construction. */
public final class VaryingCoefficientSmoothTerm {
    private VaryingCoefficientSmoothTerm() { }

    /** Multiplies every marginal basis column by the supplied by-variable. */
    public static QuadraticSmoothTerm of(
            String name, PSplineTerm marginal, double[] by) {
        double[] multiplier = MatrixOps.finiteCopy(by, "by variable");
        if (multiplier.length != marginal.observations()) {
            throw new IllegalArgumentException("by variable must match smooth observations");
        }
        QuadraticSmoothTerm base = QuadraticSmoothTerm.from(marginal);
        double[] design = base.design();
        for (int row = 0; row < base.observations(); row++) {
            for (int column = 0; column < base.columns(); column++) {
                design[row * base.columns() + column] *= multiplier[row];
            }
        }
        return new QuadraticSmoothTerm(name, base.observations(), base.columns(),
            design, base.penalties());
    }
}
