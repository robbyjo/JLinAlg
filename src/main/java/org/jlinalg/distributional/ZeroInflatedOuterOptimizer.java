/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.distributional;

/** Outer optimizers for sparse zero-inflated Laplace mixed models. */
public enum ZeroInflatedOuterOptimizer {
    /** Uses bounded BFGS for small modes and BOBYQA for large sparse modes. */
    AUTO,
    /** Uses bounded BFGS with a numerical gradient of the Laplace objective. */
    BOUNDED_BFGS,
    /** Uses the derivative-free bound-constrained BOBYQA optimizer. */
    BOBYQA
}
