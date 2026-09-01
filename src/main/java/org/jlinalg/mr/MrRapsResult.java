/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

/** Robust adjusted-profile-score MR estimate and overdispersion. */
public record MrRapsResult(MrEstimate estimate, double overdispersion,
                           int iterations, boolean converged) { }
