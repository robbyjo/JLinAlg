/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import java.util.Arrays;
import java.util.Objects;

/** Cumulative-link ordinal GEE thresholds, slopes, and underlying fit. */
public final class OrdinalGeeResult {
    private final int categories;
    private final GeeResult fit;

    OrdinalGeeResult(int categories, GeeResult fit) {
        this.categories = categories;
        this.fit = Objects.requireNonNull(fit, "fit");
    }

    public int categories() { return categories; }
    public double[] thresholds() {
        return Arrays.copyOf(fit.coefficients(), categories - 1);
    }
    public double[] coefficients() {
        double[] all = fit.coefficients();
        return Arrays.copyOfRange(all, categories - 1, all.length);
    }
    public GeeResult fit() { return fit; }
}
