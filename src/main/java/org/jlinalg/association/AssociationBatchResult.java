/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.association;

import java.util.List;
import java.util.Objects;
import org.jlinalg.inference.PValueScale;

/** Compact, input-ordered association results for a repeated model. */
public final class AssociationBatchResult {
    private final List<String> names;
    private final double[] beta;
    private final double[] standardErrors;
    private final double[] statistics;
    private final double[] pValues;
    private final double[] negativeLog10PValues;
    private final double[] degreesOfFreedom;
    private final List<AssociationFailure> failures;
    private final int coefficientIndex;
    private final int parallelism;
    private final long elapsedNanoseconds;

    AssociationBatchResult(List<String> names, double[] beta, double[] standardErrors,
            double[] statistics, double[] pValues, double[] degreesOfFreedom,
            double[] negativeLog10PValues,
            List<AssociationFailure> failures, int coefficientIndex,
            int parallelism, long elapsedNanoseconds) {
        this.names = List.copyOf(names);
        this.beta = beta.clone();
        this.standardErrors = standardErrors.clone();
        this.statistics = statistics.clone();
        this.pValues = pValues.clone();
        this.negativeLog10PValues = negativeLog10PValues.clone();
        this.degreesOfFreedom = degreesOfFreedom.clone();
        this.failures = List.copyOf(failures);
        this.coefficientIndex = coefficientIndex;
        this.parallelism = parallelism;
        this.elapsedNanoseconds = elapsedNanoseconds;
    }
    public List<String> names() { return names; }
    public double[] beta() { return beta.clone(); }
    public double[] effectSizes() { return beta.clone(); }
    public double[] standardErrors() { return standardErrors.clone(); }
    public double[] statistics() { return statistics.clone(); }
    public double[] tOrZStatistics() { return statistics.clone(); }
    public double[] pValues() { return pValues.clone(); }
    public double[] pValues(PValueScale scale) {
        Objects.requireNonNull(scale, "scale");
        return switch (scale) {
            case REGULAR -> pValues();
            case LOG10 -> log10PValues();
            case NEGATIVE_LOG10 -> negativeLog10PValues();
        };
    }
    public double[] log10PValues() {
        double[] result = negativeLog10PValues.clone();
        for (int index = 0; index < result.length; index++)
            result[index] = -result[index];
        return result;
    }
    public double[] negativeLog10PValues() {
        return negativeLog10PValues.clone();
    }
    public double[] degreesOfFreedom() { return degreesOfFreedom.clone(); }
    public List<AssociationFailure> failures() { return failures; }
    public int coefficientIndex() { return coefficientIndex; }
    public int parallelism() { return parallelism; }
    public long elapsedNanoseconds() { return elapsedNanoseconds; }
    public int size() { return names.size(); }
    public boolean successful() { return failures.isEmpty(); }

    public AssociationEstimate estimate(int index) {
        return new AssociationEstimate(names.get(index), beta[index],
            standardErrors[index], statistics[index], pValues[index],
            -negativeLog10PValues[index], negativeLog10PValues[index],
            degreesOfFreedom[index]);
    }
}
