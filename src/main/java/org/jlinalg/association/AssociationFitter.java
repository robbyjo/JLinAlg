/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.association;

import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.AssociationStatistics;

/**
 * Thread-safe adapter from row-major model data to coefficient association
 * statistics. Implementations must not retain or mutate engine-owned arrays.
 */
@FunctionalInterface
public interface AssociationFitter {
    AssociationStatistics fit(
        double[] response, double[] design, int rows, int columns,
        BackendPolicy backendPolicy);
}
