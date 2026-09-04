/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.susie;

import java.util.List;
import org.jlinalg.compute.BackendProvenance;

/** Posterior summaries from SuSiE IBSS. */
public final class SusieResult {
    private final List<String> variableNames;
    private final double[] pip;
    private final double[] posteriorMean;
    private final double[] alpha;
    private final double[] effectPosteriorMean;
    private final double[] logBayesFactors;
    private final List<CredibleSet> credibleSets;
    private final double intercept;
    private final double residualVariance;
    private final int effects;
    private final int iterations;
    private final boolean converged;
    private final double objective;
    private final BackendProvenance backend;

    SusieResult(List<String> variableNames, double[] pip, double[] posteriorMean,
            double[] alpha, double[] effectPosteriorMean,
            double[] logBayesFactors,
            List<CredibleSet> credibleSets, double intercept,
            double residualVariance, int effects, int iterations,
            boolean converged, double objective, BackendProvenance backend) {
        this.variableNames = List.copyOf(variableNames);
        this.pip = pip.clone();
        this.posteriorMean = posteriorMean.clone();
        this.alpha = alpha.clone();
        this.effectPosteriorMean = effectPosteriorMean.clone();
        this.logBayesFactors = logBayesFactors.clone();
        this.credibleSets = List.copyOf(credibleSets);
        this.intercept = intercept;
        this.residualVariance = residualVariance;
        this.effects = effects;
        this.iterations = iterations;
        this.converged = converged;
        this.objective = objective;
        this.backend = backend;
    }
    public List<String> variableNames() { return variableNames; }
    public double[] posteriorInclusionProbabilities() { return pip.clone(); }
    public double[] pip() { return pip.clone(); }
    public double[] posteriorMean() { return posteriorMean.clone(); }
    /** Row-major L-by-P single-effect inclusion probabilities. */
    public double[] alpha() { return alpha.clone(); }
    /** Row-major L-by-P unconditional single-effect means. */
    public double[] effectPosteriorMean() { return effectPosteriorMean.clone(); }
    /** Row-major L-by-P single-effect log Bayes factors. */
    public double[] logBayesFactors() { return logBayesFactors.clone(); }
    public List<CredibleSet> credibleSets() { return credibleSets; }
    public double intercept() { return intercept; }
    public double residualVariance() { return residualVariance; }
    public int effects() { return effects; }
    public int iterations() { return iterations; }
    public boolean converged() { return converged; }
    public double objective() { return objective; }
    public BackendProvenance backend() { return backend; }
}
