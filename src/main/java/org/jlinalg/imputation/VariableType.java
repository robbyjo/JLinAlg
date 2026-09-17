/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.imputation;

/** Chained-equations imputer for one numeric or integer-coded variable. */
public enum VariableType {
    /** Predictive mean matching with observed donors. */
    CONTINUOUS,
    /** Logistic regression followed by a Bernoulli draw. */
    BINARY,
    /** Integer-coded categorical probability draws from one-vs-rest logits. */
    CATEGORICAL
}
