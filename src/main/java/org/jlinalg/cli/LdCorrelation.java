/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;

/** Canonicalizes roundoff at the mathematical correlation bounds, not LD thresholds. */
final class LdCorrelation {
    private LdCorrelation() { }
    static double bounded(double value) throws IOException {
        if (!Double.isFinite(value) || Math.abs(value) > 1.0 + 1e-12)
            throw new IOException("LD correlation must be finite and bounded by one");
        return Math.max(-1.0, Math.min(1.0, value));
    }
}
