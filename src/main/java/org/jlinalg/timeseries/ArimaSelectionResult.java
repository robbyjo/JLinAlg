/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

import java.util.List;

/** Best conditional ARIMA model and all evaluated candidates. */
public record ArimaSelectionResult(
        ArimaResult bestModel, List<ArimaCandidate> candidates) {
    public ArimaSelectionResult { candidates = List.copyOf(candidates); }
}
