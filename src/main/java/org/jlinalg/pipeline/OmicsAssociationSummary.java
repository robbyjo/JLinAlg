/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

/** Counts from a bounded-memory generic omics scan. */
public record OmicsAssociationSummary(
        long sourceFeatures, long testedFeatures, long failedFeatures) { }
