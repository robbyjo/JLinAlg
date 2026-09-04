/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.coloc;

/** Posterior probabilities for one pair of SuSiE credible signals. */
public record ColocSignalPair(
        int trait1EffectIndex,
        int trait2EffectIndex,
        int variants,
        String trait1LeadVariant,
        String trait2LeadVariant,
        double posteriorH0,
        double posteriorH1,
        double posteriorH2,
        double posteriorH3,
        double posteriorH4) {

    /** Returns H0 through H4 in hypothesis order. */
    public double[] hypothesisPosteriors() {
        return new double[] {
            posteriorH0, posteriorH1, posteriorH2, posteriorH3, posteriorH4
        };
    }
}
