/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

/** One free SEM parameter with Wald inference. */
public record SemParameterEstimate(String label, double estimate,
        double standardError, double zStatistic, double pValue) { }
