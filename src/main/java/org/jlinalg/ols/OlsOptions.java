/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.ols;

import java.util.Objects;
import org.jlinalg.model.MissingDataPolicy;

/** Immutable OLS fitting options. */
public record OlsOptions(
        RankDeficiencyStrategy rankDeficiencyStrategy,
        double confidenceLevel,
        MissingDataPolicy missingDataPolicy) {

    /** Source-compatible constructor using the default strict missing policy. */
    public OlsOptions(
            RankDeficiencyStrategy rankDeficiencyStrategy,
            double confidenceLevel) {
        this(rankDeficiencyStrategy, confidenceLevel, MissingDataPolicy.ERROR);
    }

    /** Creates validated options. */
    public OlsOptions {
        Objects.requireNonNull(rankDeficiencyStrategy, "rankDeficiencyStrategy");
        Objects.requireNonNull(missingDataPolicy, "missingDataPolicy");
        if (!(confidenceLevel > 0.0 && confidenceLevel < 1.0)) {
            throw new IllegalArgumentException(
                "confidenceLevel must be strictly between zero and one");
        }
    }

    /** Returns conservative default options. */
    public static OlsOptions defaults() {
        return new OlsOptions(RankDeficiencyStrategy.ERROR, 0.95,
            MissingDataPolicy.ERROR);
    }
}
