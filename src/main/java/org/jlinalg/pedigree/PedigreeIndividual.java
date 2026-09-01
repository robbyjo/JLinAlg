/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.pedigree;

/** One individual and its optional parents in a pedigree. */
public record PedigreeIndividual(String id, String sireId, String damId) {
    /** Validates an individual. A {@code null} parent denotes an unknown parent. */
    public PedigreeIndividual {
        id = requireId(id, "id");
        sireId = optionalId(sireId, "sireId");
        damId = optionalId(damId, "damId");
        if (id.equals(sireId) || id.equals(damId)) {
            throw new IllegalArgumentException("an individual cannot be its own parent");
        }
        if (sireId != null && sireId.equals(damId)) {
            throw new IllegalArgumentException("sireId and damId must differ");
        }
    }

    /** Creates an individual whose parents are both unknown. */
    public static PedigreeIndividual founder(String id) {
        return new PedigreeIndividual(id, null, null);
    }

    private static String requireId(String value, String name) {
        if (value == null || value.isBlank()) {
            throw new IllegalArgumentException(name + " must not be null or blank");
        }
        return value;
    }

    private static String optionalId(String value, String name) {
        if (value != null && value.isBlank()) {
            throw new IllegalArgumentException(name + " must be null or nonblank");
        }
        return value;
    }
}
