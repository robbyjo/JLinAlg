/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

/** Reference distribution for coefficient and joint Wald inference. */
public enum GeeInference {
    /** Standard-normal coefficient tests and chi-square joint tests. */
    ASYMPTOTIC,
    /** Student-t and F tests using {@code clusters - parameters} denominator DF. */
    CLUSTER_T
}
