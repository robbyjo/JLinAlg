/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

import java.util.ArrayList;
import java.util.List;

/** Robust sensitivity analyses performed from one retained xWAS hit. */
public record XwasMrFollowUp(
        XwasMrHit hit,
        MrRapsResult raps,
        ContaminationMixtureResult contaminationMixture,
        MrPressoResult presso,
        List<String> warnings) {
    public XwasMrFollowUp { warnings = List.copyOf(warnings); }

    static XwasMrFollowUp analyze(XwasMrHit hit, int gridPoints,
            double pressoAlpha) {
        List<String> warnings = new ArrayList<>();
        MrRapsResult raps = null;
        ContaminationMixtureResult contamination = null;
        MrPressoResult presso = null;
        try {
            raps = RobustMendelianRandomization.raps(
                hit.harmonizedInstruments());
            if (!raps.converged()) warnings.add("RAPS did not converge");
        } catch (RuntimeException exception) {
            warnings.add("RAPS: " + message(exception));
        }
        try {
            contamination = ContaminationMixture.fit(
                hit.harmonizedInstruments(), gridPoints);
        } catch (RuntimeException exception) {
            warnings.add("contamination mixture: " + message(exception));
        }
        try {
            presso = MrPresso.analyze(hit.harmonizedInstruments(), pressoAlpha);
        } catch (RuntimeException exception) {
            warnings.add("PRESSO: " + message(exception));
        }
        return new XwasMrFollowUp(hit, raps, contamination, presso, warnings);
    }

    private static String message(RuntimeException exception) {
        return exception.getMessage() == null
            ? exception.getClass().getSimpleName() : exception.getMessage();
    }
}
