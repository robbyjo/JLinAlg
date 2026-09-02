/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

/** Final convergence metrics for the mean, association, scale, and score updates. */
public record GeeConvergenceDiagnostics(
        int iterations,
        boolean converged,
        String message,
        double coefficientChange,
        double associationChange,
        double scaleChange,
        double estimatingScoreNorm,
        double acceptedStepMultiplier) { }
