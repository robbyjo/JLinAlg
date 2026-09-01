/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.meta;

/** Coefficient inference policy for meta-analysis and meta-regression. */
public enum MetaInferenceMethod {
    /** Wald z inference. */
    NORMAL,
    /** Unadjusted Wald inference using residual Student t degrees of freedom. */
    STUDENT_T,
    /** Knapp-Hartung covariance scaling with residual Student t inference. */
    HARTUNG_KNAPP,
    /** Knapp-Hartung scaling restricted to never reduce standard errors. */
    MODIFIED_HARTUNG_KNAPP
}
