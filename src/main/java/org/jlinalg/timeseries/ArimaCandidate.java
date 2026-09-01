/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.timeseries;

/** One automatically evaluated nonseasonal ARIMA candidate. */
public record ArimaCandidate(ArimaOrder order, double aicc, double bic,
                             boolean converged) { }
