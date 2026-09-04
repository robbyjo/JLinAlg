/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.util.List;

/** Provider boundary for discovery and retrieval of MR instrument candidates. */
interface InstrumentCatalog {
    List<InstrumentStudy> search(String trait, int limit)
        throws IOException, InterruptedException;

    InstrumentStudy study(String accession)
        throws IOException, InterruptedException;

    RemoteSummary openSummaryStatistics(InstrumentStudy study)
        throws IOException, InterruptedException;

    record RemoteSummary(URI source, InputStream input) implements AutoCloseable {
        @Override public void close() throws IOException { input.close(); }
    }
}
