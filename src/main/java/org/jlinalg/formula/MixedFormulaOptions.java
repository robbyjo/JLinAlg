/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.formula;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Objects;
import org.jlinalg.model.MissingDataPolicy;
import org.jlinalg.pedigree.Pedigree;

/** Compilation controls for mixed formulas, including pedigree group mappings. */
public record MixedFormulaOptions(
        FormulaOptions formulaOptions,
        MissingDataPolicy missingDataPolicy,
        Map<String, Pedigree> pedigreeMappings) {

    public MixedFormulaOptions {
        Objects.requireNonNull(formulaOptions, "formulaOptions");
        Objects.requireNonNull(missingDataPolicy, "missingDataPolicy");
        Objects.requireNonNull(pedigreeMappings, "pedigreeMappings");
        LinkedHashMap<String, Pedigree> copy = new LinkedHashMap<>();
        pedigreeMappings.forEach((group, pedigree) -> {
            if (group == null || group.isBlank() || pedigree == null) {
                throw new IllegalArgumentException(
                    "pedigree mappings require nonblank group names and pedigrees");
            }
            copy.put(group, pedigree);
        });
        pedigreeMappings = Map.copyOf(copy);
    }

    public static MixedFormulaOptions defaults() {
        return new MixedFormulaOptions(FormulaOptions.defaults(),
            MissingDataPolicy.ERROR, Map.of());
    }

    public static Builder builder() { return new Builder(); }

    /** Builder for source-readable combinations of alignment and pedigree options. */
    public static final class Builder {
        private FormulaOptions formulaOptions = FormulaOptions.defaults();
        private MissingDataPolicy missingDataPolicy = MissingDataPolicy.ERROR;
        private final Map<String, Pedigree> pedigreeMappings = new LinkedHashMap<>();

        public Builder formulaOptions(FormulaOptions value) {
            formulaOptions = Objects.requireNonNull(value, "formulaOptions");
            return this;
        }

        public Builder missingDataPolicy(MissingDataPolicy value) {
            missingDataPolicy = Objects.requireNonNull(value, "missingDataPolicy");
            return this;
        }

        /** Maps a formula grouping column, such as {@code individual}, to a pedigree. */
        public Builder pedigree(String group, Pedigree value) {
            if (group == null || group.isBlank() || value == null) {
                throw new IllegalArgumentException(
                    "a nonblank group name and pedigree are required");
            }
            pedigreeMappings.put(group, value);
            return this;
        }

        public MixedFormulaOptions build() {
            return new MixedFormulaOptions(formulaOptions,
                missingDataPolicy, pedigreeMappings);
        }
    }
}
