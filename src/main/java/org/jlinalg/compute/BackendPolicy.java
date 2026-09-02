/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.compute;

/** Compute backend selection policies supported by JLinAlg. */
public enum BackendPolicy {
    /** CHOLMOD sparse with native dense CPU, then GPU, native CPU, and Java CPU. */
    PREFERRED,
    /** JDistlib automatic GPU/Java-CPU workload routing. */
    AUTO,
    /** Any available GPU provider, selected strictly. */
    GPU,
    /** NVIDIA CUDA, selected strictly. */
    CUDA,
    /** OpenCL, selected strictly. */
    OPENCL,
    /** Vulkan, selected strictly. */
    VULKAN,
    /** Intel oneMKL, selected strictly. */
    ONEMKL,
    /** SuiteSparse CHOLMOD with the best available dense CPU delegate. */
    CHOLMOD,
    /** OpenBLAS, selected strictly. */
    OPENBLAS,
    /** Deterministic portable Java CPU implementation. */
    CPU
}
