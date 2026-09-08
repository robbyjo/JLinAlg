/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.genetics;

import java.util.List;

/** Final joint tests for leads and conditional tests for other SNPs, plus lead-change history. */
public record ConditionalSignalSelection(List<String> selectedVariantIds,
        List<ConditionalAssociationResult> associations, List<String> changes,
        boolean converged) {
    public ConditionalSignalSelection {
        selectedVariantIds = List.copyOf(selectedVariantIds);
        associations = List.copyOf(associations);
        changes = List.copyOf(changes);
    }
}
