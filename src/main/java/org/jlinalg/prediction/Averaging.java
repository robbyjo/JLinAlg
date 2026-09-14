/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.prediction;

/** How a nonlinear response-scale estimand is formed from a design matrix. */
public enum Averaging {
    /** Transform each row to the response scale and then average the results. */
    POPULATION_AVERAGE,
    /** Average each design column first and transform that one synthetic row. */
    AT_AVERAGE_COVARIATES
}
