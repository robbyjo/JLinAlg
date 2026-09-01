/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

/** Structured per-variant fitting failure from a long-running scan. */
public record AssociationPipelineFailure(
        String variantId, String exceptionType, String message) { }
