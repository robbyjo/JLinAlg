/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.compute;

import java.util.Objects;

/** Reproducibility information for the compute backend used by a fit. */
public record BackendProvenance(
        BackendPolicy requested,
        String selectedBackend,
        String deviceDescription,
        boolean accelerated,
        boolean automaticRouting) {

    /** Creates validated backend provenance. */
    public BackendProvenance {
        Objects.requireNonNull(requested, "requested");
        Objects.requireNonNull(selectedBackend, "selectedBackend");
        Objects.requireNonNull(deviceDescription, "deviceDescription");
    }
}
