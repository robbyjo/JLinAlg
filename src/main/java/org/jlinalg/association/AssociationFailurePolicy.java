/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.association;

/** Behavior when one repeated association fit fails. */
public enum AssociationFailurePolicy {
    /** Cancel the scan and propagate the first failure. */
    FAIL_FAST,
    /** Preserve ordering, record the failure, and return NaN statistics. */
    RECORD_NAN
}
