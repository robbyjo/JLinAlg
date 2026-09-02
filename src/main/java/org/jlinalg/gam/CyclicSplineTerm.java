/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.gam;

import java.util.List;
import org.jlinalg.internal.MatrixOps;

/** Periodic Fourier representation of a cyclic cubic smooth. */
public final class CyclicSplineTerm {
    private CyclicSplineTerm() { }

    /** Creates a periodic basis whose penalty is proportional to integrated squared curvature. */
    public static QuadraticSmoothTerm of(
            String name, double[] covariate, int harmonics, double period) {
        double[] values = MatrixOps.finiteCopy(covariate, "covariate");
        if (harmonics < 1 || !(period > 0.0) || !Double.isFinite(period)) {
            throw new IllegalArgumentException("positive harmonics and period are required");
        }
        int columns = 1 + 2 * harmonics;
        double[] design = new double[values.length * columns];
        double[] penalty = new double[columns * columns];
        for (int row = 0; row < values.length; row++) {
            design[row * columns] = 1.0;
            for (int harmonic = 1; harmonic <= harmonics; harmonic++) {
                double angle = 2.0 * Math.PI * harmonic * values[row] / period;
                design[row * columns + 2 * harmonic - 1] = Math.sin(angle);
                design[row * columns + 2 * harmonic] = Math.cos(angle);
            }
        }
        for (int harmonic = 1; harmonic <= harmonics; harmonic++) {
            double frequency = 2.0 * Math.PI * harmonic / period;
            double curvature = frequency * frequency * frequency * frequency;
            penalty[(2 * harmonic - 1) * columns + 2 * harmonic - 1] = curvature;
            penalty[(2 * harmonic) * columns + 2 * harmonic] = curvature;
        }
        return new QuadraticSmoothTerm(name, values.length, columns,
            design, List.of(penalty));
    }
}
