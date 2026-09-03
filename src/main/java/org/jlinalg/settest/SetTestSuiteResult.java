/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.settest;

/** Burden, SKAT, and SKAT-O results computed from one shared score projection. */
public record SetTestSuiteResult(
        SetTestResult burden, SetTestResult skat, SkatOResult skatO) {
    public SetTestSuiteResult {
        if (burden == null || skat == null || skatO == null)
            throw new IllegalArgumentException("all set-test results are required");
    }
}
