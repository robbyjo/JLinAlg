/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.glm;

/** Controls whether GLM dispersion follows the family or Pearson estimate. */
public enum DispersionMode {
    /** Unit dispersion for binomial/Poisson; Pearson estimate for Gaussian. */
    FAMILY_DEFAULT,
    /** Estimate dispersion as Pearson chi-square divided by residual degrees of freedom. */
    PEARSON
}
