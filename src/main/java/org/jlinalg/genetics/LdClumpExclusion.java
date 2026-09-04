/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.genetics;

/** Auditable exclusion from an LD-clumping result. */
public record LdClumpExclusion(
        LdClumpCandidate candidate,
        LdClumpExclusionReason reason,
        String indexVariantId,
        double rSquared) { }
