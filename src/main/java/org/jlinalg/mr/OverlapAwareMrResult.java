/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

/** Iterative IVW estimating-equation result using exposure/outcome sampling covariance. */
public record OverlapAwareMrResult(MrEstimate estimate, int iterations,
                                   boolean converged) { }
