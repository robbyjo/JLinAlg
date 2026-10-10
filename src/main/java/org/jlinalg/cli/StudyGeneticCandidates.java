/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/** Validates and deduplicates automatic irQTL candidate exports. */
final class StudyGeneticCandidates {
    record Candidate(String feature, String id, String chromosome,
            long position, String effect, String other, double beta, double se,
            double log10p, int n) { }

    record Selection(Map<String, List<Candidate>> byFeature,
            Map<String, Integer> exclusions) { }

    private static final Comparator<Candidate> PREFERENCE = Comparator
        .comparingDouble(Candidate::log10p)
        .thenComparing(Candidate::chromosome)
        .thenComparingLong(Candidate::position)
        .thenComparing(Candidate::effect)
        .thenComparing(Candidate::other)
        .thenComparingDouble(Candidate::beta)
        .thenComparingDouble(Candidate::se)
        .thenComparingInt(Candidate::n);

    private StudyGeneticCandidates() { }

    static Selection read(Path path, Double maximumP) throws IOException {
        DelimitedData table = DelimitedData.read(path);
        Map<String, LinkedHashMap<String, Candidate>> retained = new LinkedHashMap<>();
        Map<String, Integer> exclusions = new LinkedHashMap<>();
        double cutoff = maximumP == null ? Double.POSITIVE_INFINITY
            : Math.log10(maximumP);
        for (String[] row : table.rows()) {
            String feature = row[table.column("feature_id")];
            String rawId = row[table.column("variant_id")];
            String[] tokens = rawId.split(":", -1);
            String id = tokens[tokens.length - 1];
            if (!id.matches("rs[1-9][0-9]*")) {
                increment(exclusions, "non_rsid_variant_id");
                continue;
            }
            Double beta = finite(row[table.column("beta")]);
            Double se = finite(row[table.column("se")]);
            Double logp = finite(row[table.column("log10_p")]);
            Double sampleSize = finite(row[table.column("n")]);
            if (beta == null || se == null || logp == null || sampleSize == null) {
                increment(exclusions, "nonfinite_or_invalid_numeric_value");
                continue;
            }
            if (!(se > 0)) {
                increment(exclusions, "nonpositive_standard_error");
                continue;
            }
            if (logp > 0) {
                increment(exclusions, "invalid_log10_p");
                continue;
            }
            if (logp > cutoff) {
                increment(exclusions, "above_instrument_p_threshold");
                continue;
            }
            if (sampleSize < 3 || sampleSize != Math.rint(sampleSize)
                    || sampleSize > Integer.MAX_VALUE) {
                increment(exclusions, "invalid_sample_size");
                continue;
            }
            long position;
            try {
                position = Long.parseLong(row[table.column("position")]);
            } catch (NumberFormatException exception) {
                increment(exclusions, "invalid_position");
                continue;
            }
            String chromosome = row[table.column("chromosome")]
                .replaceFirst("^chr", "");
            if (feature.isBlank() || chromosome.isBlank() || position < 1) {
                increment(exclusions, "invalid_feature_or_coordinate");
                continue;
            }
            String effect = row[table.column("effect_allele")]
                .toUpperCase(Locale.ROOT);
            String other = row[table.column("other_allele")]
                .toUpperCase(Locale.ROOT);
            if (!effect.matches("[ACGT]") || !other.matches("[ACGT]")
                    || effect.equals(other)) {
                increment(exclusions, "unsupported_non_biallelic_snp_alleles");
                continue;
            }
            Candidate candidate = new Candidate(feature, id, chromosome,
                position, effect, other, beta, se, logp,
                sampleSize.intValue());
            LinkedHashMap<String, Candidate> featureRows = retained
                .computeIfAbsent(feature, ignored -> new LinkedHashMap<>());
            Candidate previous = featureRows.get(id);
            if (previous != null) {
                increment(exclusions, "duplicate_variant_id");
                if (PREFERENCE.compare(candidate, previous) < 0)
                    featureRows.put(id, candidate);
            } else featureRows.put(id, candidate);
        }
        Map<String, List<Candidate>> byFeature = new LinkedHashMap<>();
        retained.forEach((feature, rows) -> byFeature.put(feature,
            List.copyOf(rows.values())));
        return new Selection(Map.copyOf(byFeature), Map.copyOf(exclusions));
    }

    static void writeExposure(Path output, List<Candidate> variants)
            throws IOException {
        StringBuilder table = new StringBuilder("feature_id\tvariant_id\tchromosome\tposition\teffect_allele\tother_allele\tbeta\tse\tp_value\n");
        for (Candidate value : variants) table.append(value.feature()).append('\t')
            .append(value.id()).append('\t').append(value.chromosome()).append('\t')
            .append(value.position()).append('\t').append(value.effect()).append('\t')
            .append(value.other()).append('\t').append(value.beta()).append('\t')
            .append(value.se()).append('\t').append(Math.pow(10, value.log10p()))
            .append('\n');
        Files.writeString(output, table);
    }

    private static Double finite(String value) {
        try {
            double parsed = Double.parseDouble(value.trim());
            return Double.isFinite(parsed) ? parsed : null;
        } catch (NumberFormatException exception) {
            return null;
        }
    }

    private static void increment(Map<String, Integer> counts, String reason) {
        counts.merge(reason, 1, Integer::sum);
    }
}
