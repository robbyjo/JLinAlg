/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.ewas;

/** Coordinate, signed effect, and two-sided p-value for one EWAS probe. */
public record EwasProbe(String id, String chromosome, long position,
        double effect, double pValue) { }
