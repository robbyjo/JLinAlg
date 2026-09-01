/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.association;

import java.util.Objects;
import org.jlinalg.compute.BackendPolicy;

/** Parallel execution, backend, chunking, and failure controls. */
public record AssociationEngineOptions(
        int parallelism,
        int chunkSize,
        BackendPolicy backendPolicy,
        AssociationFailurePolicy failurePolicy,
        VariableMissingPolicy predictorMissingPolicy) {
    public AssociationEngineOptions {
        if (parallelism < 1 || chunkSize < 1)
            throw new IllegalArgumentException("parallelism and chunkSize must be positive");
        Objects.requireNonNull(backendPolicy, "backendPolicy");
        Objects.requireNonNull(failurePolicy, "failurePolicy");
        Objects.requireNonNull(predictorMissingPolicy, "predictorMissingPolicy");
    }

    /** Preferred accelerator routing with one submitting thread. */
    public static AssociationEngineOptions defaults() {
        return acceleratedSerial();
    }

    /** Portable-CPU outer parallelism without nested GPU/BLAS workers. */
    public static AssociationEngineOptions cpuParallel() {
        return new AssociationEngineOptions(
            Math.max(1, Runtime.getRuntime().availableProcessors()),
            16, BackendPolicy.CPU, AssociationFailurePolicy.FAIL_FAST,
            VariableMissingPolicy.MEAN_IMPUTE);
    }

    /** Preferred accelerated backend with one submitting thread. */
    public static AssociationEngineOptions acceleratedSerial() {
        return new AssociationEngineOptions(1, 256, BackendPolicy.PREFERRED,
            AssociationFailurePolicy.FAIL_FAST, VariableMissingPolicy.MEAN_IMPUTE);
    }

    public AssociationEngineOptions withParallelism(int value) {
        return new AssociationEngineOptions(value, chunkSize, backendPolicy,
            failurePolicy, predictorMissingPolicy);
    }

    public AssociationEngineOptions withChunkSize(int value) {
        return new AssociationEngineOptions(parallelism, value, backendPolicy,
            failurePolicy, predictorMissingPolicy);
    }

    public AssociationEngineOptions withBackendPolicy(BackendPolicy value) {
        return new AssociationEngineOptions(parallelism, chunkSize, value,
            failurePolicy, predictorMissingPolicy);
    }

    public AssociationEngineOptions withFailurePolicy(
            AssociationFailurePolicy value) {
        return new AssociationEngineOptions(parallelism, chunkSize,
            backendPolicy, value, predictorMissingPolicy);
    }

    public AssociationEngineOptions withPredictorMissingPolicy(
            VariableMissingPolicy value) {
        return new AssociationEngineOptions(parallelism, chunkSize,
            backendPolicy, failurePolicy, value);
    }
}
