/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.xwas;

import jdistlib.Normal;
import org.jlinalg.genetics.GeneticCovarianceValidation;

/** Small dense matrix operations for summary models; no implicit ridge repair. */
public final class SummaryMath {
    private SummaryMath() { }
    public static void finite(double[] x) {
        if (x == null) throw new IllegalArgumentException("array is required");
        for (double v : x) if (!Double.isFinite(v))
            throw new IllegalArgumentException("all numeric inputs must be finite");
    }
    public static double dot(double[] a, double[] b) {
        if (a.length != b.length) throw new IllegalArgumentException("vector dimensions differ");
        double s = 0;
        for (int i = 0; i < a.length; i++) s = Math.fma(a[i], b[i], s);
        return s;
    }
    public static double[] multiply(double[] a, double[] x) {
        int n = x.length;
        if (a.length != n * n) throw new IllegalArgumentException("matrix dimensions differ");
        double[] y = new double[n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
            y[i] = Math.fma(a[i*n+j], x[j], y[i]);
        return y;
    }
    /** Diagonally scaled Cholesky inverse, rejecting numerically singular inputs. */
    public static double[] inverse(double[] a, int n) {
        GeneticCovarianceValidation.requirePositiveSemidefinite(a, n);
        double[] l = new double[n*n], s = new double[n];
        for (int i = 0; i < n; i++) {
            s[i] = Math.sqrt(a[i*n+i]);
            if (!(s[i] > 0)) throw new IllegalArgumentException("matrix must be positive definite");
        }
        for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) {
            double v = a[i*n+j] / s[i] / s[j];
            for (int k = 0; k < j; k++) v -= l[i*n+k] * l[j*n+k];
            if (i == j) {
                if (!(v > 1e-10)) throw new IllegalArgumentException("matrix is singular or ill-conditioned");
                l[i*n+j] = Math.sqrt(v);
            } else l[i*n+j] = v / l[j*n+j];
        }
        double[] result = new double[n*n];
        for (int c = 0; c < n; c++) {
            double[] x = new double[n];
            for (int i = 0; i < n; i++) {
                double v = i == c ? 1 : 0;
                for (int k = 0; k < i; k++) v -= l[i*n+k]*x[k];
                x[i] = v/l[i*n+i];
            }
            for (int i = n-1; i >= 0; i--) {
                double v = x[i];
                for (int k = i+1; k < n; k++) v -= l[k*n+i]*x[k];
                x[i] = v/l[i*n+i];
            }
            for (int i = 0; i < n; i++) result[i*n+c] = x[i]/s[i]/s[c];
        }
        return result;
    }
    public static double p(double z) { return 2*Normal.cumulative(-Math.abs(z), 0, 1, true, false); }
    public static double logP(double z) { return Math.log(2)+Normal.cumulative(-Math.abs(z), 0, 1, true, true); }
    /** Covariance of equal-sized delete-block estimates (blocks may differ by one row). */
    public static double[] jackknifeCovariance(double[][] deletes) {
        int b = deletes.length, p = deletes[0].length;
        double[] mean = new double[p], v = new double[p*p];
        for (double[] row : deletes) for (int j = 0; j < p; j++) mean[j] += row[j]/b;
        for (double[] row : deletes) for (int i = 0; i < p; i++) for (int j = 0; j < p; j++)
            v[i*p+j] += (b-1.0)/b*(row[i]-mean[i])*(row[j]-mean[j]);
        finite(v);
        return v;
    }
}
