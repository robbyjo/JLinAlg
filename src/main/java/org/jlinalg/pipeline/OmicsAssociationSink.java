/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.io.IOException;

/** Streaming destination for generic numeric omics association events. */
public interface OmicsAssociationSink {
    void acceptEstimate(OmicsAssociationEstimate estimate) throws IOException;
    void acceptFailure(AssociationPipelineFailure failure) throws IOException;
}
