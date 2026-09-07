/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.mediation;

import org.jlinalg.ols.OlsResult;

/**
 * Frequentist linear mediation results. The three component fits are exposed
 * so that callers can inspect the fitted models and their diagnostics.
 */
public final class MediationResult {
    private final OlsResult mediatorModel;
    private final OlsResult outcomeModel;
    private final OlsResult totalModel;
    private final MediationEffect aPath;
    private final MediationEffect bPath;
    private final MediationEffect indirectEffect;
    private final MediationEffect directEffect;
    private final MediationEffect totalEffect;
    private final int observations;
    private final int originalObservations;
    private final int[] retainedRows;

    MediationResult(
            OlsResult mediatorModel,
            OlsResult outcomeModel,
            OlsResult totalModel,
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

    /** Fit for the mediator response: {@code M ~ X + covariates}. */
    public OlsResult mediatorModel() { return mediatorModel; }

    /** Fit for the outcome response: {@code Y ~ X + M + covariates}. */
    public OlsResult outcomeModel() { return outcomeModel; }

    /** Total-effect fit: {@code Y ~ X + covariates}. */
    public OlsResult totalModel() { return totalModel; }

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

    /** Number of complete observations used by all three component fits. */
    public int observations() { return observations; }

    /** Number of rows supplied before applying the missing-data policy. */
    public int originalObservations() { return originalObservations; }

    /** Original zero-based row indices retained by the common analysis sample. */
    public int[] retainedRows() { return retainedRows.clone(); }
}
