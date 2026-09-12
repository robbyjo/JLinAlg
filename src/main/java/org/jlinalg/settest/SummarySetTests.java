/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.settest;

import java.util.List;
import jdistlib.Normal;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;

/** Tests of unweighted, nuisance-adjusted score summaries on the information scale.
 * Weights are linear burden coefficients; SKAT uses their squares as kernel weights.
 * Burden estimates use null information (one-step inference), with normal calibration.
 */
public final class SummarySetTests {
    private SummarySetTests() { }

    /** A one-variant score result. No participant-level data is required. */
    public static SetTestResult singleVariant(String id, double score, double information) {
        if (!Double.isFinite(score) || !Double.isFinite(information) || information <= 0)
            throw new IllegalArgumentException("score must be finite and information positive");
        double z = score / Math.sqrt(information);
        double logP = Math.log(2) + Normal.cumulative(-Math.abs(z), 0, 1, true, true);
        return new SetTestResult(id, "single-score", 1, 1, z,
            score / information, 1 / Math.sqrt(information), Double.NaN,
            Math.exp(logP), logP / Math.log(10), "normal-score", null, List.of());
    }

    /** Signed weighted burden, B = G w, with beta and SE per unit of B. */
    public static SetTestResult burden(String id, SetTestScoreState state, double[] weights) {
        SetTestScoreState weighted = weighted(state, weights);
        validate(weighted);
        double u = sum(weighted.scoresView()), v = sum(weighted.informationView());
        if (v == 0 && u == 0) return new SetTestResult(id,"burden-score",state.variants(),state.variants(),
            0,Double.NaN,Double.NaN,Double.NaN,1,0,"degenerate-zero-burden",null,List.of());
        if (!(v > 0)) throw new IllegalArgumentException("burden has no positive information");
        SetTestResult one = singleVariant(id, u, v);
        return new SetTestResult(id, "burden-score", state.variants(), state.variants(),
            one.statistic(), one.beta(), one.standardError(), Double.NaN,
            one.pValue(), one.log10PValue(), one.pValueMethod(), null, List.of());
    }

    /** SKAT Q = sum((w_j U_j)^2), with covariance diag(w) V diag(w). */
    public static SetTestResult skat(String id, SetTestScoreState state, double[] weights) {
        SetTestScoreState weighted = weighted(state, weights);
        validate(weighted);
        try (BackendContext context = BackendContext.select(BackendPolicy.CPU)) {
            return SetTests.kernelResult(id, "skat", state.variants(), state.variants(),
                weighted, List.of(), context.backend());
        }
    }

    /** Correlated minimum-p SKAT-O using the explicitly supplied calibration options. */
    public static SkatOResult skatO(String id, SetTestScoreState state, double[] weights,
            SetTestOptions options) {
        if (options == null) throw new IllegalArgumentException("SKAT-O options required");
        SetTestScoreState weighted = weighted(state, weights);
        validate(weighted);
        // With one variant every rho tests exactly the same hypothesis.
        if (rankOne(weighted)) {
            SetTestResult one = skat(id,state,weights);
            return new SkatOResult(id, state.variants(), state.variants(), List.of(), one.pValue(), one.pValue(),
                one.log10PValue(), 0, 0, List.of());
        }
        try (BackendContext context = BackendContext.select(BackendPolicy.CPU)) {
            return SetTests.skatO(id, state.variants(), state.variants(), List.of(),
                weighted, options, context.backend());
        }
    }

    /** Validate covariance and score consistency, including its numerical null space. */
    public static void validate(SetTestScoreState state) {
        if (state == null) throw new IllegalArgumentException("score state required");
        int n = state.variants();
        double[] u = state.scoresView(), v = state.informationView();
        double scale = 0, norm = 0;
        for (double x : u) { if (!Double.isFinite(x)) throw new IllegalArgumentException("nonfinite score"); norm = Math.hypot(norm,x); }
        for (double x : v) { if (!Double.isFinite(x)) throw new IllegalArgumentException("nonfinite covariance"); scale = Math.max(scale,Math.abs(x)); }
        if (!(scale > 0)) throw new IllegalArgumentException("no positive score information");
        if (n == 1) {
            if (v[0] <= 0) throw new IllegalArgumentException("score information must be positive");
            return;
        }
        for (int i=0;i<n;i++) for (int j=0;j<i;j++)
            if (Math.abs(v[i*n+j]-v[j*n+i]) > 1e-12*scale)
                throw new IllegalArgumentException("score covariance is not symmetric");
        try (BackendContext context = BackendContext.select(BackendPolicy.CPU)) {
            double[] normalized = v.clone();
            for(int i=0;i<normalized.length;i++) normalized[i] /= scale;
            var decomposition = context.backend().dsyev(normalized,n);
            double[] eigen = decomposition.eigenvalues(), vectors = decomposition.eigenvectors();
            for(int j=0;j<n;j++) {
                if(eigen[j] < -1e-10) throw new IllegalArgumentException("score covariance is not positive semidefinite");
                if(Math.abs(eigen[j]) <= 1e-12) {
                    double projection=0;
                    for(int i=0;i<n;i++) projection += vectors[i*n+j]*u[i];
                    if(Math.abs(projection) > 1e-7*Math.max(norm,Math.sqrt(scale)))
                        throw new IllegalArgumentException("score is inconsistent with covariance null space");
                }
            }
        }
    }

    private static SetTestScoreState weighted(SetTestScoreState state, double[] weights) {
        if(state==null || weights==null || weights.length!=state.variants())
            throw new IllegalArgumentException("one weight per variant is required");
        int n=state.variants(); double[] u=state.scores(), v=state.information();
        for(double w:weights) if(!Double.isFinite(w) || w<=0)
            throw new IllegalArgumentException("weights must be finite and positive");
        for(int i=0;i<n;i++) {
            u[i]*=weights[i];
            for(int j=0;j<n;j++) v[i*n+j]*=weights[i]*weights[j];
        }
        return new SetTestScoreState(u,v,n);
    }
    private static double sum(double[] values) {
        double sum=0, correction=0;
        for(double value:values) { double t=sum+value; correction+=Math.abs(sum)>=Math.abs(value)?(sum-t)+value:(value-t)+sum; sum=t; }
        return sum+correction;
    }
    private static boolean rankOne(SetTestScoreState state) {
        int n=state.variants();double[] v=state.informationView();
        for(int i=0;i<n;i++)for(int j=0;j<i;j++) {
            double target=Math.sqrt(v[i*n+i])*Math.sqrt(v[j*n+j]);
            if(Math.abs(Math.abs(v[i*n+j])-target)>1e-12*target)return false;
        }
        return true;
    }
}
