/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.formula.CompiledFormula;
import org.jlinalg.glm.GlmFamily;

/** Formula adapter for fitting GEE from a reusable compiled model matrix. */
public final class GeeFormula {
    private GeeFormula() { }

    /** Fits a GEE without reparsing the supplied compiled formula. */
    public static GeeResult fit(
            CompiledFormula formula,
            int[] cluster,
            int[] repeated,
            GlmFamily family,
            GeeOptions options,
            BackendPolicy backendPolicy) {
        if (formula == null) {
            throw new IllegalArgumentException("formula is required");
        }
        return Gee.fit(formula.response(), formula.design(),
            formula.rows(), formula.columns(), cluster, repeated,
            family, formula.weights(), formula.offset(), options, backendPolicy);
    }
}
