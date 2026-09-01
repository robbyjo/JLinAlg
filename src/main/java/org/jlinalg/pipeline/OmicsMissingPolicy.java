/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

/** Fast-path handling for missing values in an omics feature row. */
public enum OmicsMissingPolicy {
    ERROR,
    MEAN_IMPUTE
}
