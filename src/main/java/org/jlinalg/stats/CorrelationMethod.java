/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.stats;

/** Correlation coefficient and test family. */
public enum CorrelationMethod {
    /** Pearson product-moment correlation. */
    PEARSON,
    /** Kendall rank correlation tau-b. */
    KENDALL,
    /** Spearman rank correlation rho. */
    SPEARMAN
}
