/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.distributional;

import java.util.List;

/** Deterministic parametric-bootstrap outer estimates and fit accounting. */
public record ZeroInflatedBootstrapResult(
        List<String> parameterNames,
        double[] estimates,
        int replicates,
        int successfulReplicates) {
    public ZeroInflatedBootstrapResult {
        parameterNames = List.copyOf(parameterNames);
        if (estimates == null
                || estimates.length != replicates * parameterNames.size()
                || successfulReplicates < 0
                || successfulReplicates > replicates)
            throw new IllegalArgumentException("bootstrap dimensions are invalid");
        estimates = estimates.clone();
    }

    @Override public double[] estimates() { return estimates.clone(); }

    /** One optimizer-scale parameter across replicates, with NaN on failure. */
    public double[] estimates(String parameter) {
        int column = parameterNames.indexOf(parameter);
        if (column < 0)
            throw new IllegalArgumentException("unknown parameter: " + parameter);
        double[] result = new double[replicates];
        for (int row = 0; row < replicates; row++)
            result[row] = estimates[row * parameterNames.size() + column];
        return result;
    }
}
