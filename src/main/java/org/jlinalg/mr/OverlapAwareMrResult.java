/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

/** Errors-in-variables IVW estimate using exposure/outcome sampling covariance. */
public record OverlapAwareMrResult(MrEstimate estimate, int iterations,
                                   boolean converged) { }
