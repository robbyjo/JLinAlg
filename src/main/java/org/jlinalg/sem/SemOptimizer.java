/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.sem;

/** BFGS on a per-observation objective; convergence requires a small score. */
final class SemOptimizer {
    interface Objective { Value at(double[] point); }
    record Value(double value, double[] gradient) { }
    record Optimum(double[] point, int evaluations, boolean converged, double scoreNorm) { }
    static Optimum minimize(Objective f, double[] start, int maximum, double tolerance) {
        double[] x = start.clone(); int k = x.length, evaluations = 1;
        Value v = f.at(x);
        if (!Double.isFinite(v.value()) || v.gradient() == null)
            throw new IllegalArgumentException("SEM starting values do not define a valid distribution");
        double[] h = RamFit.identity(k);
        while (evaluations < maximum && norm(v.gradient()) > tolerance) {
            double[] d = RamFit.mv(h, v.gradient());
            for (int j = 0; j < k; j++) d[j] = -d[j];
            double slope = dot(d, v.gradient());
            if (!(slope < 0)) {
                h = RamFit.identity(k);
                for (int j = 0; j < k; j++) d[j] = -v.gradient()[j];
                slope = -dot(d, d);
            }
            Value next = null; double[] y = null;
            for (double alpha = 1; alpha >= 1e-14 && evaluations < maximum; alpha *= .5) {
                y = x.clone(); for (int j = 0; j < k; j++) y[j] += alpha * d[j];
                Value trial = f.at(y); evaluations++;
                if (Double.isFinite(trial.value()) && trial.gradient() != null
                        && trial.value() <= v.value() + 1e-4 * alpha * slope) {
                    next = trial; break;
                }
            }
            if (next == null) break;
            double[] s = new double[k], z = new double[k];
            for (int j = 0; j < k; j++) { s[j] = y[j] - x[j]; z[j] = next.gradient()[j] - v.gradient()[j]; }
            double sz = dot(s, z);
            if (sz > 1e-12 * Math.sqrt(dot(s,s) * dot(z,z))) {
                double[] hz = RamFit.mv(h, z); double c = (sz + dot(z,hz)) / (sz*sz);
                for (int i=0;i<k;i++) for(int j=0;j<k;j++)
                    h[i*k+j] += c*s[i]*s[j] - (hz[i]*s[j]+s[i]*hz[j])/sz;
            } else h = RamFit.identity(k);
            x=y; v=next;
        }
        return new Optimum(x, evaluations, norm(v.gradient()) <= tolerance, norm(v.gradient()));
    }
    static double dot(double[] a,double[] b) { double s=0; for(int i=0;i<a.length;i++) s+=a[i]*b[i]; return s; }
    static double norm(double[] a) { double s=0; for(double v:a) s=Math.max(s,Math.abs(v)); return s; }
}
