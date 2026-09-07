/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mediation;

import org.jlinalg.mixed.SparseLinearMixedModelResult;

/**
 * Frequentist mediation results from sparse Gaussian REML component fits.
 * Component models may contain ordinary or pedigree random-effect terms.
 */
public final class MediationMixedResult {
    private final SparseLinearMixedModelResult mediatorModel;
    private final SparseLinearMixedModelResult outcomeModel;
    private final SparseLinearMixedModelResult totalModel;
    private final MediationEffect aPath;
    private final MediationEffect bPath;
    private final MediationEffect indirectEffect;
    private final MediationEffect directEffect;
    private final MediationEffect totalEffect;
    private final int observations;
    private final int originalObservations;
    private final int[] retainedRows;

    MediationMixedResult(
            SparseLinearMixedModelResult mediatorModel,
            SparseLinearMixedModelResult outcomeModel,
            SparseLinearMixedModelResult totalModel,
            MediationEffect aPath,
            MediationEffect bPath,
            MediationEffect indirectEffect,
            MediationEffect directEffect,
            MediationEffect totalEffect,
            int observations,
            int originalObservations,
            int[] retainedRows) {
        this.mediatorModel = mediatorModel;
        this.outcomeModel = outcomeModel;
        this.totalModel = totalModel;
        this.aPath = aPath;
        this.bPath = bPath;
        this.indirectEffect = indirectEffect;
        this.directEffect = directEffect;
        this.totalEffect = totalEffect;
        this.observations = observations;
        this.originalObservations = originalObservations;
        this.retainedRows = retainedRows.clone();
    }

    /** Sparse REML fit for {@code M ~ X + covariates + random effects}. */
    public SparseLinearMixedModelResult mediatorModel() {
        return mediatorModel;
    }

    /** Sparse REML fit for {@code Y ~ X + M + covariates + random effects}. */
    public SparseLinearMixedModelResult outcomeModel() {
        return outcomeModel;
    }

    /** Sparse REML total-effect fit for {@code Y ~ X + covariates}. */
    public SparseLinearMixedModelResult totalModel() {
        return totalModel;
    }

    /** Treatment-to-mediator path, conventionally called {@code a}. */
    public MediationEffect aPath() { return aPath; }

    /** Mediator-to-outcome path adjusted for treatment, conventionally {@code b}. */
    public MediationEffect bPath() { return bPath; }

    /** Product-of-coefficients indirect effect, {@code a * b}. */
    public MediationEffect indirectEffect() { return indirectEffect; }

    /** Outcome-model treatment coefficient, conventionally {@code c'}. */
    public MediationEffect directEffect() { return directEffect; }

    /** Total-model treatment coefficient, conventionally {@code c}. */
    public MediationEffect totalEffect() { return totalEffect; }

    /** Whether every sparse REML component fit converged. */
    public boolean converged() {
        return mediatorModel.converged()
            && outcomeModel.converged() && totalModel.converged();
    }

    /** Number of observations used by all three component fits. */
    public int observations() { return observations; }

    /** Number of rows supplied before the common complete-row check. */
    public int originalObservations() { return originalObservations; }

    /** Original zero-based row indices retained by the common analysis sample. */
    public int[] retainedRows() { return retainedRows.clone(); }
}
