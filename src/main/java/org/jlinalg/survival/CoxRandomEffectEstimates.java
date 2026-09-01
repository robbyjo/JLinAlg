/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.survival;

import java.util.List;

/** Conditional Gaussian frailty modes for one Cox random-effect term. */
public record CoxRandomEffectEstimates(
        String termName,
        List<String> coefficientNames,
        double variance,
        double[] modes) {
    public CoxRandomEffectEstimates {
        coefficientNames = List.copyOf(coefficientNames);
        modes = modes.clone();
        if (termName == null || coefficientNames.size() != modes.length
                || !(variance >= 0))
            throw new IllegalArgumentException("invalid Cox frailty estimates");
    }
    @Override public double[] modes() { return modes.clone(); }
}
