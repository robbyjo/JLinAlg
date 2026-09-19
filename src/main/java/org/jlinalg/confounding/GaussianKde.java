/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.confounding;

import java.util.Arrays;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

/** Gaussian KDE using linear binning, zero-padded FFT convolution, and linear
 * interpolation. Grid spacing is at most bandwidth/32. */
final class GaussianKde {
    private GaussianKde() { }
    static double[] atObservations(double[] x, double bandwidth) {
        double min = Arrays.stream(x).min().orElseThrow();
        double max = Arrays.stream(x).max().orElseThrow();
        if (x.length <= 256 || min == max) return direct(x, bandwidth);
        double requested = Math.ceil((max - min) / bandwidth * 32);
        if (!(requested <= 131072))
            throw new IllegalArgumentException("SVA KDE bandwidth is too narrow for the bounded grid; rescale or inspect degenerate feature probabilities");
        int intervals = Math.max(1024, (int) requested);
        double step = (max - min) / intervals;
        int size = 1;
        while (size < 2 * (intervals + 1)) size *= 2;
        double[] mass = new double[size], kernel = new double[size];
        for (double value : x) {
            double position = Math.min(intervals, (value - min) / step);
            int left = Math.min(intervals - 1, (int) position);
            double fraction = position - left;
            mass[left] += (1 - fraction) / x.length;
            mass[left + 1] += fraction / x.length;
        }
        double normalizer = bandwidth * Math.sqrt(2 * Math.PI);
        kernel[0] = 1 / normalizer;
        for (int i = 1; i <= intervals; i++) {
            double z = i * step / bandwidth;
            kernel[i] = kernel[size - i] = Math.exp(-0.5 * z * z) / normalizer;
        }
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] a = fft.transform(mass, TransformType.FORWARD);
        Complex[] b = fft.transform(kernel, TransformType.FORWARD);
        for (int i = 0; i < size; i++) a[i] = a[i].multiply(b[i]);
        Complex[] density = fft.transform(a, TransformType.INVERSE);
        double[] result = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            double position = Math.min(intervals, (x[i] - min) / step);
            int left = Math.min(intervals - 1, (int) position);
            double fraction = position - left;
            result[i] = (1 - fraction) * density[left].getReal()
                + fraction * density[left + 1].getReal();
            if (!(result[i] > 0) || !Double.isFinite(result[i]))
                throw new IllegalArgumentException("SVA KDE lost numerical precision");
        }
        return result;
    }
    private static double[] direct(double[] x, double bandwidth) {
        double[] result = new double[x.length];
        if (Arrays.stream(x).allMatch(v -> v == x[0])) {
            Arrays.fill(result, 1 / (bandwidth * Math.sqrt(2 * Math.PI))); return result;
        }
        for (int i = 0; i < x.length; i++) {
            for (double value : x) {
                double z = (x[i] - value) / bandwidth;
                result[i] += Math.exp(-0.5 * z * z);
            }
            result[i] /= x.length * bandwidth * Math.sqrt(2 * Math.PI);
        }
        return result;
    }
}
