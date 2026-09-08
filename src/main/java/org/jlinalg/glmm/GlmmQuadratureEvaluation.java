/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.glmm;

/** Likelihood and diagnostics at fixed parameters. Error is a refinement
 * estimate, not a certified bound. Nodes is the largest rule used by any group.
 * A single rule cannot establish convergence. */
public record GlmmQuadratureEvaluation(double logLikelihood, int nodes,
        double estimatedError, boolean converged) { }
