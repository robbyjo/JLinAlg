/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

/**
 * Final convergence metrics for mean, association, scale, and score updates.
 * The score norm is max(abs(score[j]) / sqrt(sensitivity[j,j])), so changing
 * coefficient or response units does not change the score convergence test.
 */
public record GeeConvergenceDiagnostics(
        int iterations,
        boolean converged,
        String message,
        double coefficientChange,
        double associationChange,
        double scaleChange,
        double estimatingScoreNorm,
        double acceptedStepMultiplier) { }
