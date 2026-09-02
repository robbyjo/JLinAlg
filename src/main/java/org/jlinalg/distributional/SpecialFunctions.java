/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.distributional;

/** Stable scalar special functions needed by distributional likelihoods. */
final class SpecialFunctions {
    private static final double[] LANCZOS = {
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    };

    private SpecialFunctions() { }

    static double logGamma(double value) {
        if (!(value > 0.0)) return Double.NaN;
        if (value < 0.5) {
            return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * value))
                - logGamma(1.0 - value);
        }
        double shifted = value - 1.0;
        double sum = 0.99999999999980993;
        for (int index = 0; index < LANCZOS.length; index++) {
            sum += LANCZOS[index] / (shifted + index + 1.0);
        }
        double temporary = shifted + LANCZOS.length - 0.5;
        return 0.5 * Math.log(2.0 * Math.PI)
            + (shifted + 0.5) * Math.log(temporary)
            - temporary + Math.log(sum);
    }

    static double digamma(double value) {
        if (!(value > 0.0)) return Double.NaN;
        double result = 0.0;
        double x = value;
        while (x < 8.0) {
            result -= 1.0 / x;
            x += 1.0;
        }
        double inverse = 1.0 / x;
        double squared = inverse * inverse;
        return result + Math.log(x) - 0.5 * inverse
            - squared * (1.0 / 12.0
                - squared * (1.0 / 120.0
                    - squared * (1.0 / 252.0
                        - squared * (1.0 / 240.0))));
    }

    static double trigamma(double value) {
        if (!(value > 0.0)) return Double.NaN;
        double result = 0.0;
        double x = value;
        while (x < 8.0) {
            result += 1.0 / (x * x);
            x += 1.0;
        }
        double inverse = 1.0 / x;
        double squared = inverse * inverse;
        return result + inverse + 0.5 * squared
            + inverse * squared / 6.0
            - inverse * squared * squared / 30.0
            + inverse * squared * squared * squared / 42.0
            - inverse * squared * squared * squared * squared / 30.0;
    }
}
