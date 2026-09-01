/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pipeline;

import java.io.IOException;

/** Stateful bounded reader for row-feature numeric matrices. */
public interface NumericBlockReader extends AutoCloseable {
    NumericBlock read(int maximumRows) throws IOException;
    @Override void close() throws IOException;
}
