/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.loess;

/** Residual-weighting rule for LOESS. */
public enum LoessFamily {
    /** Ordinary locally weighted least squares. */
    GAUSSIAN,
    /** Tukey-bisquare residual reweighting, matching R's symmetric family. */
    SYMMETRIC
}
