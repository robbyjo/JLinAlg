/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.io.IOException;

/** Incremental destination for bounded-memory association and QC output. */
public interface AssociationPipelineSink {
    void acceptEstimate(AssociationPipelineEstimate estimate) throws IOException;
    void acceptExcluded(VariantFilterResult excluded) throws IOException;
    void acceptFailure(AssociationPipelineFailure failure) throws IOException;
}
