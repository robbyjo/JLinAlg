/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

/** Intercept/mean estimates for a SEM data matrix. */
public final class SemMeanStructure {
    private SemMeanStructure() { }

    public static Result fit(double[][] data) {
        if (data == null || data.length < 2 || data[0] == null) throw new IllegalArgumentException("mean-structure data are invalid");
        int observations = data.length, variables = data[0].length; double[] means = new double[variables];
        for (double[] row : data) { if (row == null || row.length != variables) throw new IllegalArgumentException("mean-structure data must be rectangular"); for (int variable = 0; variable < variables; variable++) { if (!Double.isFinite(row[variable])) throw new IllegalArgumentException("mean-structure data require complete rows"); means[variable] += row[variable]; } }
        for (int variable = 0; variable < variables; variable++) means[variable] /= observations;
        double[] covariance = new double[variables * variables]; for (double[] row : data) for (int left = 0; left < variables; left++) for (int right = 0; right < variables; right++) covariance[left * variables + right] += (row[left] - means[left]) * (row[right] - means[right]) / Math.max(1, observations - 1);
        double[] standardErrors = new double[variables]; for (int variable = 0; variable < variables; variable++) standardErrors[variable] = Math.sqrt(Math.max(0.0, covariance[variable * variables + variable] / observations));
        return new Result(means, standardErrors, covariance, observations);
    }

    public record Result(double[] intercepts, double[] standardErrors, double[] covariance, int observations) {
        public Result { intercepts = intercepts.clone(); standardErrors = standardErrors.clone(); covariance = covariance.clone(); }
        public double[] intercepts() { return intercepts.clone(); }
        public double[] standardErrors() { return standardErrors.clone(); }
        public double[] covariance() { return covariance.clone(); }
    }
}
