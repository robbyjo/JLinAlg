/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.io.IOException;

/** Stateful, bounded-memory variant-block reader. */
public interface VariantBlockReader extends AutoCloseable {
    /** Returns the next block, or null at end of input. */
    VariantBlock read(int maximumVariants) throws IOException;

    @Override void close() throws IOException;
}
