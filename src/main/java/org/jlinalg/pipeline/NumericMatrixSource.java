/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.io.IOException;

/** Re-openable source of numeric feature rows aligned by sample ID. */
public interface NumericMatrixSource {
    NumericMatrixMetadata metadata();
    NumericBlockReader open(int[] sampleOrder) throws IOException;
    default NumericBlockReader open() throws IOException { return open(null); }
}
