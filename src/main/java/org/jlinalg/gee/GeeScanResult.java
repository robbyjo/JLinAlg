/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gee;

import java.util.List;
import org.jlinalg.compute.BackendProvenance;

/** Compact coefficient-level output from a prepared marginal GEE scan. */
public final class GeeScanResult {
    private final List<String> names;
    private final double[] coefficients;
    private final double[] standardErrors;
    private final double[] statistics;
    private final double[] pValues;
    private final double[] associationParameters;
    private final int[] iterations;
    private final boolean[] converged;
    private final int parallelism;
    private final long elapsedNanoseconds;
    private final BackendProvenance backend;

    GeeScanResult(
            List<String> names,
            double[] coefficients,
            double[] standardErrors,
            double[] statistics,
            double[] pValues,
            double[] associationParameters,
            int[] iterations,
            boolean[] converged,
            int parallelism,
            long elapsedNanoseconds,
            BackendProvenance backend) {
        this.names = List.copyOf(names);
        this.coefficients = coefficients.clone();
        this.standardErrors = standardErrors.clone();
        this.statistics = statistics.clone();
        this.pValues = pValues.clone();
        this.associationParameters = associationParameters.clone();
        this.iterations = iterations.clone();
        this.converged = converged.clone();
        this.parallelism = parallelism;
        this.elapsedNanoseconds = elapsedNanoseconds;
        this.backend = backend;
    }

    public List<String> names() { return names; }
    public double[] coefficients() { return coefficients.clone(); }
    public double[] standardErrors() { return standardErrors.clone(); }
    public double[] statistics() { return statistics.clone(); }
    public double[] pValues() { return pValues.clone(); }
    /** First working-association estimate per scan, or NaN for independence. */
    public double[] associationParameters() {
        return associationParameters.clone();
    }
    public int[] iterations() { return iterations.clone(); }
    public boolean[] converged() { return converged.clone(); }
    public int parallelism() { return parallelism; }
    public long elapsedNanoseconds() { return elapsedNanoseconds; }
    public BackendProvenance backend() { return backend; }
    public int size() { return coefficients.length; }
}
