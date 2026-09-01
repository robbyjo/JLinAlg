/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.susie;

import java.util.List;

/** One single-effect credible set and its minimum absolute LD. */
public record CredibleSet(
        int effectIndex,
        List<String> variables,
        double posteriorCoverage,
        double purity) {
    public CredibleSet { variables = List.copyOf(variables); }
}
