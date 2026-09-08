/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

import org.jlinalg.compute.BackendPolicy;

/** The portable RAM engine can honor CPU and policies permitting CPU fallback. */
final class SemBackendPolicy {
    private SemBackendPolicy() { }

    static void requireSupported(BackendPolicy policy) {
        switch (policy) {
            case CPU, PREFERRED, AUTO -> { }
            default -> throw new UnsupportedOperationException(
                "SEM RAM engine does not support explicit backend policy " + policy
                    + "; use CPU, PREFERRED, or AUTO for portable Java CPU execution");
        }
    }
}
