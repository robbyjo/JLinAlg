/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.compute;

import java.util.Objects;
import jdistlib.accelerator.Compute;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.ComputeBackends;
import jdistlib.accelerator.ComputeSelection;

/** Owns a JDistlib compute selection and closes its native resources. */
public final class BackendContext implements AutoCloseable {
    private final BackendPolicy requested;
    private final ComputeSelection selection;
    private final ComputeBackend backend;

    private BackendContext(BackendPolicy requested, ComputeSelection selection,
            ComputeBackend backend) {
        this.requested = Objects.requireNonNull(requested, "requested");
        this.selection = Objects.requireNonNull(selection, "selection");
        this.backend = Objects.requireNonNull(backend, "backend");
    }

    private BackendContext(BackendPolicy requested, ComputeSelection selection) {
        this(requested, selection, selection.backend());
    }

    /** Selects the requested policy. Strict policies fail when unavailable. */
    public static BackendContext select(BackendPolicy policy) {
        Objects.requireNonNull(policy, "policy");
        if (policy == BackendPolicy.PREFERRED) {
            return preferred();
        }
        if (policy == BackendPolicy.CHOLMOD) return cholmod(policy);
        return new BackendContext(policy, ComputeBackends.select(toJdistlib(policy)));
    }

    /**
     * Selects CHOLMOD with the best native dense CPU backend when available,
     * then GPU automatic routing, oneMKL, OpenBLAS, and portable Java CPU.
     */
    public static BackendContext preferred() {
        try {
            return cholmod(BackendPolicy.PREFERRED);
        } catch (IllegalStateException | LinkageError unavailable) {
            // Continue through the portable/native dense fallback chain.
        }
        if (available(Compute.GPU)) {
            return new BackendContext(
                BackendPolicy.PREFERRED, ComputeBackends.select(Compute.AUTO));
        }
        if (available(Compute.ONEMKL)) {
            return new BackendContext(
                BackendPolicy.PREFERRED, ComputeBackends.select(Compute.ONEMKL));
        }
        if (available(Compute.OPENBLAS)) {
            return new BackendContext(
                BackendPolicy.PREFERRED, ComputeBackends.select(Compute.OPENBLAS));
        }
        return new BackendContext(
            BackendPolicy.PREFERRED, ComputeBackends.select(Compute.CPU));
    }

    /** Returns the selected JDistlib backend. The context retains ownership. */
    public ComputeBackend backend() {
        return backend;
    }

    /** Returns immutable backend provenance for result objects. */
    public BackendProvenance provenance() {
        return new BackendProvenance(
            requested,
            backend.selectedBackend(),
            backend.deviceInfo().description(),
            backend.capabilities().nativeFactorizations()
                || backend.capabilities().nativeSparseFactorizations(),
            backend.automaticRouting());
    }

    @Override
    public void close() {
        selection.close();
    }

    private static boolean available(Compute compute) {
        try (ComputeSelection probe = ComputeBackends.select(compute)) {
            probe.selectedBackend();
            return true;
        } catch (IllegalStateException | LinkageError exception) {
            return false;
        }
    }

    private static BackendContext cholmod(BackendPolicy requested) {
        ComputeSelection dense = selectDenseCpu();
        try {
            // Load the dense native runtime first. On Windows the CHOLMOD
            // bridge can share its oneMKL runtime dependency.
            CholmodNative.requireAvailable();
            return new BackendContext(requested, dense,
                new CholmodComputeBackend(dense.backend()));
        } catch (RuntimeException | LinkageError failure) {
            dense.close();
            throw failure;
        }
    }

    private static ComputeSelection selectDenseCpu() {
        if (available(Compute.ONEMKL)) return ComputeBackends.select(Compute.ONEMKL);
        if (available(Compute.OPENBLAS)) return ComputeBackends.select(Compute.OPENBLAS);
        return ComputeBackends.select(Compute.CPU);
    }

    private static Compute toJdistlib(BackendPolicy policy) {
        return switch (policy) {
            case AUTO -> Compute.AUTO;
            case GPU -> Compute.GPU;
            case CUDA -> Compute.CUDA;
            case OPENCL -> Compute.OPENCL;
            case VULKAN -> Compute.VULKAN;
            case ONEMKL -> Compute.ONEMKL;
            case OPENBLAS -> Compute.OPENBLAS;
            case CPU -> Compute.CPU;
            case PREFERRED, CHOLMOD -> throw new IllegalArgumentException(
                policy + " must be resolved by the JLinAlg fallback policy");
        };
    }
}
