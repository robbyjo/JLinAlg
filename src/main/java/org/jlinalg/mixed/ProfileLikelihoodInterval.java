/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.mixed;

/** Profile-likelihood confidence interval for one scalar mixed-model parameter. */
public record ProfileLikelihoodInterval(double estimate, double lower,
                                       double upper, double cutoff,
                                       boolean lowerFound, boolean upperFound) { }
