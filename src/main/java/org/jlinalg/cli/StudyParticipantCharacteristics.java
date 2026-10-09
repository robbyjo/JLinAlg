/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.jlinalg.cli.StudyPublicEvidence.*;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.jlinalg.study.StudyNeed;
import org.jlinalg.study.StudyOperationExecutor;

/** Local-only Table 1 computation. No participant values or identifiers are emitted. */
final class StudyParticipantCharacteristics implements StudyOperationExecutor {
    @Override public Result execute(Context context) throws IOException {
        Files.createDirectories(context.outputDirectory());
        List<Object> tables = new ArrayList<>(), needs = new ArrayList<>();
        for (Map<String, Object> dataset : StudyBioinformatics.datasets(context)) {
            if (!map(dataset.get("files")).containsKey("phenotype")) {
                needs.add(need(dataset, "Supply files.phenotype and columns.characteristics for Table 1."));
                continue;
            }
            try {
                Map<String, Object> summary = summarize(context, dataset);
                tables.add(summary);
                if (!Boolean.TRUE.equals(summary.get("analysis_population_verified")))
                    needs.add(need(dataset, "Confirm the analysis population using an explicit membership file or a documented confirmation that the phenotype file contains exactly the analysis sample. File-level summaries are provisional, not analysis-sample demographics."));
                if (!list(summary.get("missing_requirements")).isEmpty())
                    needs.add(need(dataset, "Resolve Table 1 requirements: " + summary.get("missing_requirements")));
            } catch (IllegalArgumentException | IOException exception) {
                tables.add(Map.of("dataset", dataset.get("id"), "status", "input_required", "reason", exception.getMessage()));
                needs.add(need(dataset, exception.getMessage()));
            }
        }
        Map<String, Object> document = new LinkedHashMap<>();
        document.put("tables", tables); document.put("needs", needs);
        document.put("privacy", "Computed locally. Output contains aggregate summaries and input hashes only; not automatically included in remote LLM payloads.");
        document.put("method", "Per-variable available-case means and sample standard deviations (n-1); missing counts and categorical denominators retained. No imputation. Unverified populations are not publication-ready.");
        write(context.outputDirectory().resolve("participant-characteristics.yaml"), document);
        return new Result(needs.isEmpty() ? State.COMPLETE : State.WAITING_FOR_INPUT,
            "Computed " + tables.size() + " local participant-summary tables; " + needs.size() + " review requirements",
            List.of("participant-characteristics.yaml"), document);
    }

    private static Object need(Map<String, Object> dataset, String reason) {
        return StudyNeed.review("table1:" + dataset.get("id") + ":" + Integer.toHexString(reason.hashCode()),
            "Can you verify the participant-characteristics inputs for " + dataset.get("id") + "?", reason, true).document();
    }

    static Map<String, Object> summarize(Context context, Map<String, Object> dataset) throws IOException {
        Map<String, Object> files = map(dataset.get("files")), columns = map(dataset.get("columns"));
        Map<String, Object> config = map(columns.get("characteristics"));
        Set<String> allowed = Set.of("id", "age", "sex", "male_value", "female_value", "sex_constant", "sex_basis",
            "main_phenotype", "continuous", "categorical", "units", "membership_role", "membership_id_column",
            "analysis_population_confirmed", "population_rationale", "filters", "derived_binary");
        for (String key : config.keySet()) if (!allowed.contains(key)) throw new IllegalArgumentException("Unknown characteristics field: " + key);
        Path source = context.project().resolve(text(files.get("phenotype"))).normalize();
        DelimitedData table = DelimitedData.read(source);
        String beforeHash = FollowupSupport.hash(source);
        List<String[]> rows = new ArrayList<>(table.rows());
        Map<String, Object> filters = map(config.get("filters"));
        for (var filter : filters.entrySet()) {
            int col = table.column(filter.getKey());
            Set<String> values = new HashSet<>();
            for (Object item : list(filter.getValue())) values.add(text(item));
            if (values.isEmpty()) throw new IllegalArgumentException("Each participant filter requires an explicit nonempty list of accepted values");
            rows.removeIf(row -> !values.contains(row[col]));
        }
        Map<String, Object> result = new LinkedHashMap<>();
        result.put("dataset", dataset.get("id")); result.put("cohort", dataset.getOrDefault("cohort", "not_declared"));
        result.put("phenotype_sha256", beforeHash); result.put("input_rows", table.rows().size());
        // A digest of the selection specification allows safe cross-hypothesis deduplication
        // without publishing participant identifiers or filter values.
        result.put("selection_sha256", org.jlinalg.study.StudyReportIdentity.selection(config));
        result.put("selection", filters.isEmpty() ? "All phenotype rows; no implicit association sample reconstruction" : "Explicit column-value filters from study contract");
        String idName = text(config.get("id"));
        boolean verified = Boolean.TRUE.equals(config.get("analysis_population_confirmed"));
        String basis = text(config.get("population_rationale"));
        if (verified && basis.isBlank()) throw new IllegalArgumentException("Analysis-population confirmation needs population_rationale");
        if (config.containsKey("membership_role")) {
            if (idName.isBlank()) throw new IllegalArgumentException("Membership alignment requires characteristics.id");
            String role = text(config.get("membership_role"));
            if (!files.containsKey(role)) throw new IllegalArgumentException("Membership file role is missing");
            Path membershipPath = context.project().resolve(text(files.get(role))).normalize();
            DelimitedData members = DelimitedData.read(membershipPath);
            int mc = members.column(text(config.getOrDefault("membership_id_column", idName))), pc = table.column(idName);
            Set<String> ids = new HashSet<>();
            for (String[] row : members.rows()) if (missing(row[mc]) || !ids.add(row[mc]))
                throw new IllegalArgumentException("Membership identifiers must be nonmissing and unique; identifiers are not printed");
            rows.removeIf(row -> !ids.contains(row[pc]));
            if (rows.stream().map(row -> row[pc]).distinct().count() != ids.size())
                throw new IllegalArgumentException("Some analysis-membership identifiers are absent from the selected phenotype records");
            result.put("membership_sha256", FollowupSupport.hash(membershipPath));
            if (basis.isBlank()) throw new IllegalArgumentException("Membership selection requires a population_rationale describing how the analysis sample was obtained");
            verified = true;
        }
        if (!idName.isBlank()) {
            int id = table.column(idName); Set<String> ids = new HashSet<>();
            for (String[] row : rows) if (missing(row[id]) || !ids.add(row[id]))
                throw new IllegalArgumentException("Participant IDs must be nonmissing and unique after selection; resolve repeated visits before Table 1");
        } else if (verified) throw new IllegalArgumentException("Verified participant counts require an explicit ID column");
        if (rows.isEmpty()) throw new IllegalArgumentException("Participant selection is empty");
        result.put("N", rows.size()); result.put("N_unit", idName.isBlank() ? "records; participant identity unverified" : "unique participants in selected phenotype records");
        result.put("analysis_population_verified", verified); result.put("population_rationale", basis.isBlank() ? "Not supplied" : basis);
        result.put("status", verified ? "analysis_population_declared" : "provisional_phenotype_file_population");
        Map<String, Object> continuous = new LinkedHashMap<>(), categorical = new LinkedHashMap<>();
        List<String> requirements = new ArrayList<>();
        String age = text(config.getOrDefault("age", table.header().contains("Age") ? "Age" : ""));
        if (age.isBlank()) requirements.add("Explicit age column at the relevant study visit");
        else continuous.put("Age", numeric(table, rows, age));
        String outcome = text(StudyBioinformatics.hypothesis(context).get("outcome"));
        String phenotype = text(config.getOrDefault("main_phenotype", map(columns.get("phenotypes")).getOrDefault(outcome,
            table.header().contains(outcome) ? outcome : "")));
        if (phenotype.isBlank()) requirements.add("Main phenotype column and definition (do not treat absent data as zero)");
        else continuous.put("Main phenotype", numeric(table, rows, phenotype));
        result.put("main_phenotype_label", config.containsKey("main_phenotype") ? phenotype : outcome);
        for (var column : map(config.get("continuous")).entrySet()) continuous.put(column.getKey(), numeric(table, rows, text(column.getValue())));
        for (var column : map(config.get("categorical")).entrySet()) categorical.put(column.getKey(), categories(table, rows, text(column.getValue())));
        Map<String, Object> binary = new LinkedHashMap<>();
        for (var entry : map(config.get("derived_binary")).entrySet()) binary.put(entry.getKey(), binary(table, rows, map(entry.getValue())));
        result.put("binary", binary);
        String sex = text(config.getOrDefault("sex", table.header().contains("Sex") ? "Sex" : ""));
        Map<String, Object> sexSummary = new LinkedHashMap<>();
        if (config.containsKey("sex_constant")) {
            String constant = text(config.get("sex_constant"));
            if (!Set.of("female", "male").contains(constant) || text(config.get("sex_basis")).isBlank())
                throw new IllegalArgumentException("sex_constant requires male/female and a documented sex_basis");
            sexSummary.put("male_n", constant.equals("male") ? rows.size() : 0);
            sexSummary.put("male_percent", constant.equals("male") ? 100.0 : 0.0);
            sexSummary.put("denominator", rows.size()); sexSummary.put("missing", 0); sexSummary.put("basis", config.get("sex_basis"));
        } else if (!sex.isBlank() && config.containsKey("male_value") && config.containsKey("female_value")) {
            String male = text(config.get("male_value")), female = text(config.get("female_value"));
            if (male.equals(female)) throw new IllegalArgumentException("Male and female codes must differ");
            int col = table.column(sex), males = 0, known = 0, unknown = 0;
            for (String[] row : rows) {
                if (row[col].equals(male)) { males++; known++; }
                else if (row[col].equals(female)) known++;
                else if (!missing(row[col])) unknown++;
            }
            sexSummary.put("male_n", males); sexSummary.put("denominator", known);
            sexSummary.put("male_percent", known == 0 ? "not_estimable" : 100.0 * males / known);
            sexSummary.put("missing", rows.size() - known - unknown); sexSummary.put("unmapped", unknown); sexSummary.put("column", sex);
            sexSummary.put("male_value", male); sexSummary.put("female_value", female);
            sexSummary.put("basis", config.getOrDefault("sex_basis", "Explicit study-contract coding"));
            if (unknown > 0) requirements.add("Unmapped sex categories require a documented reporting decision");
        } else requirements.add("Sex column with explicit male/female coding, or a documented single-sex cohort declaration");
        result.put("sex", sexSummary); result.put("continuous", continuous); result.put("categorical", categorical);
        result.put("units", map(config.get("units"))); result.put("missing_requirements", requirements);
        if (!beforeHash.equals(FollowupSupport.hash(source))) throw new IOException("Phenotype input changed during summary computation");
        return result;
    }

    static Map<String, Object> numeric(DelimitedData table, List<String[]> rows, String column) {
        int col = table.column(column), n = 0, missingCount = 0; double mean = 0, m2 = 0;
        for (String[] row : rows) {
            if (missing(row[col])) { missingCount++; continue; }
            double value;
            try { value = Double.parseDouble(row[col]); } catch (NumberFormatException e) { throw new IllegalArgumentException("Nonnumeric value in configured continuous column " + column + "; values are not printed"); }
            if (!Double.isFinite(value)) throw new IllegalArgumentException("Nonfinite value in continuous column " + column);
            double delta = value - mean; mean += delta / ++n; m2 += delta * (value - mean);
        }
        return Map.of("column", column, "n", n, "missing", missingCount, "mean", n == 0 ? "not_estimable" : mean,
            "sd", n < 2 ? "not_estimable" : Math.sqrt(Math.max(0, m2 / (n - 1))));
    }

    private static Map<String, Object> binary(DelimitedData table, List<String[]> rows, Map<String, Object> rule) {
        if (!rule.keySet().equals(Set.of("column", "operator", "threshold", "definition"))
                || !Set.of("gt", "ge", "lt", "le").contains(text(rule.get("operator")))
                || !(rule.get("threshold") instanceof Number threshold) || !Double.isFinite(threshold.doubleValue())
                || text(rule.get("definition")).isBlank())
            throw new IllegalArgumentException("derived_binary requires column, operator (gt/ge/lt/le), finite numeric threshold, and definition");
        String column = text(rule.get("column")), operator = text(rule.get("operator"));
        int col = table.column(column), n = 0, cases = 0;
        for (String[] row : rows) {
            if (missing(row[col])) continue;
            double value;
            try { value = Double.parseDouble(row[col]); }
            catch (NumberFormatException e) { throw new IllegalArgumentException("Nonnumeric value in derived-binary source column " + column); }
            if (!Double.isFinite(value)) throw new IllegalArgumentException("Nonfinite value in derived-binary source column " + column);
            boolean positive = switch (operator) {
                case "gt" -> value > threshold.doubleValue(); case "ge" -> value >= threshold.doubleValue();
                case "lt" -> value < threshold.doubleValue(); default -> value <= threshold.doubleValue();
            };
            n++; if (positive) cases++;
        }
        Map<String, Object> result = new LinkedHashMap<>(rule);
        result.put("case_n", cases); result.put("denominator", n); result.put("missing", rows.size() - n);
        result.put("percent", n == 0 ? "not_estimable" : 100.0 * cases / n);
        return result;
    }

    private static Map<String, Object> categories(DelimitedData table, List<String[]> rows, String column) {
        int col = table.column(column), absent = 0; Map<String, Integer> counts = new LinkedHashMap<>();
        for (String[] row : rows) {
            if (missing(row[col])) absent++;
            else counts.merge(row[col], 1, Integer::sum);
        }
        if (counts.size() > 20) throw new IllegalArgumentException("Categorical summary exceeded 20 levels; do not summarize identifier columns");
        return Map.of("column", column, "denominator", rows.size() - absent, "missing", absent, "counts", counts);
    }

    private static boolean missing(String value) { return value == null || Set.of("", "NA", "N/A", "NAN", "NULL", ".").contains(value.trim().toUpperCase(java.util.Locale.ROOT)); }
}
