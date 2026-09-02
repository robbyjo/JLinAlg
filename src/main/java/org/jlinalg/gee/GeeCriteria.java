/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

/** Quasi-likelihood criteria for mean and working-correlation selection. */
public record GeeCriteria(
        double quasiLikelihood,
        double qic,
        double qicu,
        double cic,
        double qicc,
        int parameters) {
}
