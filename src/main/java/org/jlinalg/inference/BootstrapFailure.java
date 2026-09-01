/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.inference;

/** Structured failed or nonconverged bootstrap replicate. */
public record BootstrapFailure(
        int simulation,
        String exceptionType,
        String message) {
    public BootstrapFailure {
        if (simulation < 0 || exceptionType == null || message == null)
            throw new IllegalArgumentException("invalid bootstrap failure");
    }
}
