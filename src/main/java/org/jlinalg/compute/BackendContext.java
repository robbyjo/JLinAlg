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

    private BackendContext(BackendPolicy requested, ComputeSelection selection) {
        this.requested = Objects.requireNonNull(requested, "requested");
        this.selection = Objects.requireNonNull(selection, "selection");
    }

    /** Selects the requested policy. Strict policies fail when unavailable. */
    public static BackendContext select(BackendPolicy policy) {
        Objects.requireNonNull(policy, "policy");
        if (policy == BackendPolicy.PREFERRED) {
            return preferred();
        }
        return new BackendContext(policy, ComputeBackends.select(toJdistlib(policy)));
    }

    /**
     * Selects GPU with JDistlib's automatic size routing when a GPU exists,
     * otherwise oneMKL, OpenBLAS, and finally the portable Java CPU.
     */
    public static BackendContext preferred() {
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
        return selection.backend();
    }

    /** Returns immutable backend provenance for result objects. */
    public BackendProvenance provenance() {
        return new BackendProvenance(
            requested,
            selection.selectedBackend(),
            selection.deviceInfo().description(),
            selection.accelerated(),
            selection.automaticRouting());
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
            case PREFERRED -> throw new IllegalArgumentException(
                "PREFERRED must be resolved by the JLinAlg fallback policy");
        };
    }
}
