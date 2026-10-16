/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Audited, P-value-blind QC of generated association tables. */
final class StudyEffectQc {
    private StudyEffectQc() { }
    private static final double MAD_NORMAL = 1.482602218505602;
    private static final double IQR_NORMAL = 1.348979500392164;

    record Result(Map<String, Object> document, List<String> outputs) { }

    static Result apply(Map<String, Path> tables, Map<String, Object> harmonization,
            Path directory) throws IOException {
        Map<String, Object> settings = StudyPublicEvidence.map(harmonization.get("effect_outliers"));
        Object enabledValue = settings.getOrDefault("enabled", true);
        if (!(enabledValue instanceof Boolean)) throw new IllegalArgumentException("effect_outliers.enabled must be boolean");
        boolean enabled = Boolean.TRUE.equals(enabledValue);
        double cutoff = positive(settings.getOrDefault("robust_z", 100), "robust_z");
        double minimum = positive(settings.getOrDefault("min_features", 20), "min_features");
        if (minimum < 3 || minimum != Math.rint(minimum) || minimum > Integer.MAX_VALUE)
            throw new IllegalArgumentException("effect_outliers.min_features must be an integer >= 3");
        double absolute = settings.containsKey("max_abs_effect")
            ? positive(settings.get("max_abs_effect"), "max_abs_effect") : Double.POSITIVE_INFINITY;
        List<Object> summaries = new ArrayList<>();
        StringBuilder summary = new StringBuilder("dataset_id\ttested_rows\tmedian_beta\trobust_scale\tscale_method\trobust_z_cutoff\tmax_abs_effect\tstatus\texcluded_rows\tretained_rows\n");
        StringBuilder exclusions = new StringBuilder("dataset_id\tfeature_id\tbeta\tse\tp_value\trobust_z\treason\n");
        int total = 0;
        for (var entry : tables.entrySet()) {
            DelimitedData table = DelimitedData.read(entry.getValue());
            int beta = table.column("beta"), id = table.column("feature_id");
            double[] values = table.rows().stream().mapToDouble(row -> Double.parseDouble(row[beta])).sorted().toArray();
            double center = quantile(values, .5);
            double[] deviations = Arrays.stream(values).map(value -> Math.abs(value - center)).sorted().toArray();
            double scale = quantile(deviations, .5) * MAD_NORMAL;
            String scaleMethod = "normal_consistent_MAD";
            if (!(scale > 0) || !Double.isFinite(scale)) {
                scale = (quantile(values, .75) - quantile(values, .25)) / IQR_NORMAL;
                scaleMethod = "normal_consistent_IQR";
            }
            boolean relative = enabled && values.length >= minimum && scale > 0 && Double.isFinite(scale);
            String status = !enabled ? "disabled" : values.length < minimum ? "insufficient_features"
                : !relative ? "zero_or_nonfinite_robust_scale" : "ready";
            List<String[]> retained = new ArrayList<>();
            int removed = 0;
            for (String[] row : table.rows()) {
                double effect = Double.parseDouble(row[beta]);
                double score = relative ? Math.abs(effect - center) / scale : Double.NaN;
                List<String> reasons = new ArrayList<>();
                if (enabled && Math.abs(effect) > absolute) reasons.add("absolute_effect_limit");
                if (relative && score > cutoff) reasons.add("extreme_robust_effect");
                if (reasons.isEmpty()) retained.add(row);
                else {
                    removed++;
                    exclusions.append(entry.getKey()).append('\t').append(row[id]).append('\t')
                        .append(row[beta]).append('\t').append(row[table.column("se")]).append('\t')
                        .append(row[table.column("p_value")]).append('\t').append(number(score)).append('\t')
                        .append(String.join(";", reasons)).append('\n');
                }
            }
            if (removed > 0) try (BufferedWriter writer = Files.newBufferedWriter(entry.getValue())) {
                writer.write(String.join("\t", table.header()) + "\n");
                for (String[] row : retained) writer.write(String.join("\t", row) + "\n");
            }
            total += removed;
            Map<String, Object> item = new LinkedHashMap<>();
            item.put("dataset_id", entry.getKey()); item.put("tested_rows", values.length);
            item.put("median_beta", number(center)); item.put("robust_scale", number(scale));
            item.put("scale_method", scaleMethod); item.put("status", status);
            item.put("excluded_rows", removed); item.put("retained_rows", retained.size()); summaries.add(item);
            summary.append(entry.getKey()).append('\t').append(values.length).append('\t').append(number(center))
                .append('\t').append(number(scale)).append('\t').append(scaleMethod).append('\t').append(cutoff)
                .append('\t').append(number(absolute)).append('\t').append(status).append('\t').append(removed)
                .append('\t').append(retained.size()).append('\n');
        }
        Files.writeString(directory.resolve("effect-outlier-summary.tsv"), summary);
        Files.writeString(directory.resolve("effect-outlier-exclusions.tsv"), exclusions);
        Map<String, Object> parameters = new LinkedHashMap<>();
        parameters.put("enabled", enabled); parameters.put("robust_z", cutoff); parameters.put("min_features", (int) minimum);
        parameters.put("max_abs_effect", Double.isFinite(absolute) ? absolute : "not_set");
        parameters.put("comparison", "strictly_greater_than"); parameters.put("scope", "per_dataset_association_estimate");
        Map<String, Object> audit = new LinkedHashMap<>();
        audit.put("method", "conservative robust effect-size QC"); audit.put("parameters", parameters);
        audit.put("performed", enabled && !tables.isEmpty());
        audit.put("software", "JLinAlg " + JLinAlgCli.version());
        audit.put("parameter_justification", "Software QC heuristic; no universal biological cutoff or scientist-specific justification is inferred. Review units, exclusions and an unfiltered sensitivity analysis, especially for a post hoc change.");
        audit.put("datasets", summaries); audit.put("excluded_rows", total);
        audit.put("interpretation", "Scale-aware heuristic, not proof of invalidity or a causal finding. Computed on all numerically eligible rows before FDR, annotation or QTL matching, without using P values, effect significance or cross-cohort agreement. MAD scale=1.482602218505602*median(abs(beta-median(beta))); zero MAD falls back to IQR/1.348979500392164. Insufficient or zero-scale families skip relative screening; an explicit absolute bound still applies. Inspect units, predictor variation and model diagnostics. Raw inputs are unchanged. Excluded estimates are absent from every downstream standardized table; common_universe=intersection also removes their features from all cohorts. Historical sealed outputs are not repaired automatically.");
        StudyPublicEvidence.write(directory.resolve("effect-outlier-qc.yaml"), audit);
        if (summaries.stream().map(StudyPublicEvidence::map).anyMatch(item -> ((Number) item.get("retained_rows")).intValue() == 0))
            throw new IOException("No analyzable association rows after effect-outlier QC; inspect effect-outlier-qc.yaml");
        return new Result(audit, List.of("effect-outlier-summary.tsv", "effect-outlier-exclusions.tsv", "effect-outlier-qc.yaml"));
    }

    private static double positive(Object value, String key) {
        if (!(value instanceof Number number) || !Double.isFinite(number.doubleValue()) || !(number.doubleValue() > 0))
            throw new IllegalArgumentException("effect_outliers." + key + " must be a positive finite number");
        return number.doubleValue();
    }

    private static String number(double value) { return Double.isFinite(value) ? Double.toString(value) : "NA"; }

    private static double quantile(double[] sorted, double probability) {
        if (sorted.length == 0) return Double.NaN;
        double index = (sorted.length - 1) * probability;
        int left = (int) Math.floor(index), right = (int) Math.ceil(index);
        double weight = index - left;
        return sorted[left] * (1 - weight) + sorted[right] * weight;
    }
}
