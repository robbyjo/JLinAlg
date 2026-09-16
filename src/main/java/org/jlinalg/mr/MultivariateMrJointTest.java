/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mr;

/** Wald test that a named set of multivariate MR coefficients is jointly zero. */
public record MultivariateMrJointTest(
        String name, double chiSquare, int degreesOfFreedom, double pValue) { }
