/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.net.URI;
import java.util.List;

/** Search metadata for a GWAS with downloadable summary statistics. */
record InstrumentStudy(
        String accession,
        String trait,
        List<String> ontologyTraits,
        List<String> ancestry,
        String sampleDescription,
        long variantCount,
        URI summaryStatistics,
        URI license) {
    InstrumentStudy {
        ontologyTraits = List.copyOf(ontologyTraits);
        ancestry = List.copyOf(ancestry);
    }
}
