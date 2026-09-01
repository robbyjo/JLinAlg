/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mr;

/** Fixed-effect IVW estimate after omitting one named instrument. */
public record LeaveOneOutEstimate(String omittedVariantId, MrEstimate estimate) {
}
