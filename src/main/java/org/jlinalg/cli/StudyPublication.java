/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.jlinalg.cli.StudyPublicEvidence.*;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jlinalg.study.StudyEvidenceIndex;
import org.jlinalg.study.StudyOperationExecutor.Context;

/** Deterministic, local, source-linked publication scaffold. No LLM or upload needed. */
final class StudyPublication {
    private StudyPublication() { }

    static List<String> write(Context context) throws IOException {
        Map<String, Object> hypothesis = StudyBioinformatics.hypothesis(context);
        List<Map<String, Object>> records = new ArrayList<>();
        List<Object> tables = new ArrayList<>();
        StringBuilder all = new StringBuilder(), caveats = new StringBuilder();
        for (var dependency : context.dependencyDirectories().entrySet()) {
            Path directory = dependency.getValue(), manifestPath = directory.resolve("manifest.yaml");
            if (!Files.isRegularFile(manifestPath)) continue;
            Map<String, Object> manifest = load(manifestPath);
            String type = text(manifest.get("operation_type")), state = text(manifest.get("state"));
            Object findings = StudyEvidenceIndex.keyFindings(directory, type);
            Map<String, Object> record = new LinkedHashMap<>();
            record.put("operation", type); record.put("state", state); record.put("summary", manifest.getOrDefault("summary", ""));
            record.put("details", manifest.getOrDefault("details", Map.of()));
            record.put("manifest", relative(context.outputDirectory(), manifestPath)); record.put("manifest_sha256", FollowupSupport.hash(manifestPath));
            record.put("output_sha256", manifest.getOrDefault("output_sha256", Map.of()));
            record.put("findings", findings == null ? Map.of() : findings);
            records.add(record);
            if (type.equals("participant_characteristics")) tables.addAll(list(map(findings).get("tables")));
            if (type.equals("diagnostic_review") || !state.equals("complete")) {
                caveats.append("<li><strong>").append(esc(label(type))).append(" — ").append(esc(state)).append("</strong>: ")
                    .append(esc(text(manifest.get("summary")))).append(" ").append(link(record)).append("</li>");
                if (type.equals("diagnostic_review")) caveats.append("<li>").append(render(findings)).append("</li>");
            }
            all.append("<section><h2>").append(esc(label(type))).append("</h2><p>").append(esc(state)).append(" — ")
                .append(esc(text(manifest.get("summary")))).append(" ").append(link(record)).append("</p>")
                .append(render(findings)).append("</section>");
        }
        Map<String, Object> ledger = new LinkedHashMap<>();
        ledger.put("hypothesis", hypothesis.getOrDefault("description", "Scientific hypothesis description required"));
        ledger.put("writing_policy", StudyResultsDiscussion.policy());
        ledger.put("impact_ranking", "Not assigned automatically. A scientist or advisory LLM must select distinct findings, link Results paragraphs and justify their impact in an audited revision. Analysis order is not scientific importance.");
        ledger.put("records", records);
        Map<String, Map<String, Object>> references = new LinkedHashMap<>();
        for (Map<String, Object> r : records) {
            Map<String, Object> source = new LinkedHashMap<>();
            source.put("operation_type", r.get("operation")); source.put("key_findings", r.get("findings")); source.put("manifest", r.get("manifest"));
            source.put("state", r.get("state")); source.put("details", r.get("details"));
            source.put("summary", r.get("summary"));
            references.put(text(hypothesis.get("id")) + "." + r.get("operation"), source);
        }
        Map<String, Object> catalog = StudyReferenceCatalog.fromEvidence(references);
        StudyCitations citations = new StudyCitations(catalog);
        Map<String, Object> introduction = StudyIntroduction.template("Publication first draft", List.of(text(ledger.get("hypothesis"))));
        Map<String, Object> introductionReview = new LinkedHashMap<>(StudyIntroduction.assess(list(introduction.get("sections"))));
        var scopedStudy = new LinkedHashMap<>(context.study()); scopedStudy.put("datasets", StudyBioinformatics.datasets(context));
        var methodsInventory = StudyMethods.inventory(scopedStudy, references);
        var methods = StudyMethods.template("Publication first draft", methodsInventory, catalog);
        var methodsReview = StudyMethods.assess(list(methods.get("sections")), methodsInventory);
        var narrative = StudyResultsDiscussion.template("Publication first draft", references);
        var narrativeReview = StudyResultsDiscussion.assess(narrative, references);
        StudyPublicEvidence.write(context.outputDirectory().resolve("results-discussion-inventory.yaml"), StudyResultsDiscussion.inventory(references));
        StudyPublicEvidence.write(context.outputDirectory().resolve("results-discussion-template.yaml"), narrative);
        StudyPublicEvidence.write(context.outputDirectory().resolve("results-discussion-review.yaml"), narrativeReview);
        ledger.put("results_discussion_review", narrativeReview);
        StudyPublicEvidence.write(context.outputDirectory().resolve("methods-inventory.yaml"), methodsInventory);
        StudyPublicEvidence.write(context.outputDirectory().resolve("methods-template.yaml"), methods);
        StudyPublicEvidence.write(context.outputDirectory().resolve("methods-review.yaml"), methodsReview);
        ledger.put("methods_review", methodsReview);
        introductionReview.put("literature_search_coverage", catalog.get("search_coverage"));
        ledger.put("introduction_review", introductionReview);
        StudyPublicEvidence.write(context.outputDirectory().resolve("introduction-template.yaml"), introduction);
        StudyPublicEvidence.write(context.outputDirectory().resolve("introduction-review.yaml"), introductionReview);
        StudyPublicEvidence.write(context.outputDirectory().resolve("reference-catalog.yaml"), catalog);
        Files.writeString(context.outputDirectory().resolve("reference-catalog.html"), StudyReferenceCatalog.html(catalog));
        StudyPublicEvidence.write(context.outputDirectory().resolve("publication-claims.yaml"), ledger);
        StringBuilder html = new StringBuilder(start("Publication first draft"));
        html.append("<p><strong>Unreviewed first draft generated locally by JLinAlg.</strong> No language model was called. Missing information remains explicit; this is not ready for submission.</p>")
            .append("<p><a href=\"introduction-template.yaml\">Introduction writing template</a> · <a href=\"introduction-review.yaml\">Missing introduction inputs</a> · <a href=\"reference-catalog.html\">Reference catalog and search coverage</a></p>")
            .append("<h2>Introduction</h2>").append(StudyIntroduction.scaffold(list(introduction.get("sections"))))
            .append("<h2>Methods</h2>").append(StudyMethods.scaffold(methods, citations))
            .append("<p><a href=\"methods-template.yaml\">Methods writing template</a> · <a href=\"methods-review.yaml\">Methods review needs</a> · <a href=\"methods-inventory.yaml\">Method and parameter provenance</a></p>")
            .append("<details><summary>Declared cohort analysis descriptions and execution inventory</summary>").append(render(StudyMethods.displayInventory(methodsInventory))).append("</details>");
        html.append("<h2>Results</h2><h3>Table 1. Participant characteristics</h3>")
            .append("<p>SD denotes sample standard deviation. Provisional phenotype-file populations are not verified analysis samples.</p>")
            .append(table1(tables));
        for (Object raw : list(narrative.get("sections"))) {
            var p = map(raw); if (!"results".equals(p.get("section"))) continue;
            var source = references.get(list(p.get("sources")).get(0).toString());
            html.append("<section id=\"").append(esc(text(p.get("id")))).append("\"><h3>").append(esc(text(p.get("title")))).append("</h3><p>")
                .append(esc(text(source.get("summary")))).append(" ").append(link(source)).append("</p>")
                .append(StudyIntroduction.scaffold(List.of(p))).append(render(source.get("key_findings"))).append("</section>");
        }
        html.append("<h2>Discussion</h2>");
        for (Object raw : list(narrative.get("sections"))) {
            var p = map(raw); if (!"discussion".equals(p.get("section"))) continue;
            html.append(StudyIntroduction.scaffold(List.of(p)));
            if ("limitations".equals(p.get("discussion_role"))) html.append("<details><summary>Unresolved evidence and review needs</summary><ul>").append(caveats)
                .append("</ul><p>Unavailable evidence is not negative evidence. Exploratory partial-data colocalization posteriors are conditional on observed selected variants; incomplete coverage, significance filtering, unequal coverage and winner's curse can bias them, so a null result is not evidence against a shared causal variant. Association, gene-level annotations and literature matches do not establish isoform-specific causality or drug efficacy.</p></details>");
        }
        html.append("<h2>Conclusion</h2>").append(StudyIntroduction.scaffold(list(narrative.get("sections")).stream().filter(p -> "conclusion".equals(map(p).get("section"))).toList()))
            .append("<p><a href=\"results-discussion-template.yaml\">Results and Discussion writing template</a> · <a href=\"results-discussion-review.yaml\">Narrative review needs</a> · <a href=\"results-discussion-inventory.yaml\">Findings and bioinformatics context inventory</a></p>")
            .append("<p><a href=\"publication-supplement.html\">All evidence layers and supplemental findings</a> · <a href=\"publication-claims.yaml\">Source and writing-policy ledger</a></p>").append(citations.bibliography()).append("</body></html>");
        Files.writeString(context.outputDirectory().resolve("publication-draft.html"), expandAbbreviations(html.toString(),
            map(map(hypothesis.get("metadata")).get("abbreviations"))));
        Files.writeString(context.outputDirectory().resolve("publication-supplement.html"), start("Supplemental evidence")
            + "<h2>Supplemental Table 1. Detailed participant characteristics</h2>" + detailedTable1(tables)
            + "<p>All evidence layers are retained, including null, contradictory, blocked and unreviewed results. Top-row previews are not exhaustive; use the linked manifests for complete artifacts.</p>" + all + "</body></html>");
        Files.writeString(context.outputDirectory().resolve("participant-characteristics.csv"), table1Csv(tables));
        return List.of("publication-draft.html", "publication-supplement.html", "publication-claims.yaml", "participant-characteristics.csv",
            "introduction-template.yaml", "introduction-review.yaml", "reference-catalog.yaml", "reference-catalog.html",
            "methods-template.yaml", "methods-review.yaml", "methods-inventory.yaml",
            "results-discussion-template.yaml", "results-discussion-review.yaml", "results-discussion-inventory.yaml");
    }
    static String table1(List<Object> tables) {
        if (tables.isEmpty()) return "<p>[SCIENTIST INPUT REQUIRED: phenotype files and verified analysis-sample membership]</p>";
        List<Map<String, Object>> columns = new ArrayList<>();
        for (Object raw : tables) {
            Map<String, Object> t = new LinkedHashMap<>(map(raw));
            Map<String, Object> continuous = new LinkedHashMap<>(map(t.get("continuous")));
            if (continuous.containsKey("Main phenotype")) {
                String name = text(t.get("main_phenotype_label"));
                if (name.isBlank()) name = text(map(continuous.get("Main phenotype")).get("column"));
                continuous.put(name.isBlank() ? "Main phenotype" : name, continuous.remove("Main phenotype"));
            }
            t.put("continuous", continuous); t.put("binary", new LinkedHashMap<>(map(t.get("binary"))));
            Map<String, Object> match = columns.stream().filter(c -> compatible(c, t)).findFirst().orElse(null);
            if (match == null) columns.add(t);
            else for (String field : List.of("continuous", "binary")) {
                Map<String, Object> merged = new LinkedHashMap<>(map(match.get(field))); merged.putAll(map(t.get(field))); match.put(field, merged);
            }
        }
        StringBuilder html = new StringBuilder("<table id=\"table-1\"><thead><tr><th scope=\"col\">Characteristic</th>");
        for (Map<String, Object> t : columns) {
            String name = text(t.getOrDefault("cohort", t.get("dataset")));
            boolean repeated = columns.stream().filter(c -> java.util.Objects.equals(c.get("cohort"), t.get("cohort"))).count() > 1;
            html.append("<th scope=\"col\">").append(esc(name.toUpperCase(java.util.Locale.ROOT)))
                .append(repeated ? " (" + esc(text(t.get("dataset"))) + ")" : "").append("</th>");
        }
        html.append("</tr></thead><tbody>");
        row(html, "N", columns.stream().map(t -> text(t.getOrDefault("N", "Input required"))).toList());
        row(html, "Male, n (%)", columns.stream().map(t -> proportion(map(t.get("sex")), "male_n", "male_percent", t.get("N"))).toList());
        java.util.Set<String> continuous = new java.util.LinkedHashSet<>(), binary = new java.util.LinkedHashSet<>();
        columns.forEach(t -> { continuous.addAll(map(t.get("continuous")).keySet()); binary.addAll(map(t.get("binary")).keySet()); });
        // A declared binary phenotype is a proportion, not the mean and SD of a 0/1 variable.
        continuous.removeAll(binary);
        for (String name : continuous) {
            java.util.Set<String> units = new java.util.LinkedHashSet<>();
            columns.forEach(t -> { String unit = text(map(t.get("units")).get(name)); if (!unit.isBlank()) units.add(unit); });
            String unit = units.size() == 1 ? " (" + units.iterator().next() + ")" : "";
            row(html, name + unit + ", mean ± SD", columns.stream().map(t -> mean(map(t.get("continuous")).get(name), t.get("N"))).toList());
        }
        for (String name : binary) row(html, name + ", n (%)", columns.stream().map(t -> proportion(map(map(t.get("binary")).get(name)), "case_n", "percent", t.get("N"))).toList());
        html.append("</tbody></table><p>Available-case summaries; when a variable has fewer observations than N, its reduced n is shown. Detailed denominators, missingness, coding and source hashes are in Supplemental Table 1. Display values are rounded to two decimal places; CSV and evidence records retain full precision.</p>");
        java.util.Set<String> definitions = new java.util.LinkedHashSet<>();
        columns.forEach(t -> map(t.get("binary")).values().forEach(v -> definitions.add(text(map(v).get("definition")))));
        for (String definition : definitions) html.append("<p>").append(esc(definition)).append("</p>");
        if (columns.stream().anyMatch(t -> !Boolean.TRUE.equals(t.get("analysis_population_verified"))))
            html.append("<p>Provisional phenotype-file populations: analysis-sample membership remains unverified.</p>");
        return html.toString();
    }
    private static boolean compatible(Map<String, Object> a, Map<String, Object> b) {
        if (!a.containsKey("selection_sha256") || !b.containsKey("selection_sha256")) return false;
        for (String key : List.of("cohort", "phenotype_sha256", "selection_sha256", "membership_sha256", "N", "N_unit", "analysis_population_verified", "sex", "units"))
            if (!java.util.Objects.equals(a.get(key), b.get(key))) return false;
        for (String field : List.of("continuous", "binary")) for (var entry : map(a.get(field)).entrySet())
            if (map(b.get(field)).containsKey(entry.getKey()) && !entry.getValue().equals(map(b.get(field)).get(entry.getKey()))) return false;
        return true;
    }
    private static void row(StringBuilder html, String name, List<String> cells) {
        html.append("<tr><th scope=\"row\">").append(esc(name)).append("</th>");
        cells.forEach(cell -> html.append("<td>").append(esc(cell)).append("</td>")); html.append("</tr>");
    }
    private static String denominator(Object n, Object total) { return java.util.Objects.equals(n, total) ? "" : " [n=" + text(n) + "]"; }
    private static String proportion(Map<String, Object> m, String count, String percent, Object total) {
        return m.isEmpty() ? "Input required" : text(m.get(count)) + " (" + StudyReportNumbers.format(m.get(percent)) + "%)" + denominator(m.get("denominator"), total);
    }
    private static String mean(Object value, Object total) { Map<String, Object> m = map(value); return m.isEmpty() ? "Input required" : StudyReportNumbers.format(m.get("mean")) + " ± " + StudyReportNumbers.format(m.get("sd")) + denominator(m.get("n"), total); }
    static String detailedTable1(List<Object> tables) {
        StringBuilder h = new StringBuilder("<p><a href=\"participant-characteristics.csv\">Full-precision participant-summary CSV</a></p><table><thead><tr><th>Dataset</th><th>Characteristic</th><th>Available n</th><th>Missing</th><th>Mean or n (%)</th><th>SD</th></tr></thead><tbody>");
        for (Object raw : tables) {
            Map<String, Object> t = map(raw);
            for (var entry : map(t.get("continuous")).entrySet()) {
                Map<String, Object> v = map(entry.getValue());
                detailRow(h, List.of(text(t.get("dataset")), entry.getKey().equals("Main phenotype") ? text(t.getOrDefault("main_phenotype_label", "Main phenotype")) : entry.getKey(), text(v.get("n")), text(v.get("missing")), StudyReportNumbers.format(v.get("mean")), StudyReportNumbers.format(v.get("sd"))));
            }
            Map<String, Object> sex = map(t.get("sex"));
            if (!sex.isEmpty()) detailRow(h, List.of(text(t.get("dataset")), "Male", text(sex.get("denominator")), text(sex.get("missing")), proportion(sex, "male_n", "male_percent", sex.get("denominator")), "—"));
            for (var entry : map(t.get("binary")).entrySet()) {
                Map<String, Object> v = map(entry.getValue());
                detailRow(h, List.of(text(t.get("dataset")), entry.getKey(), text(v.get("denominator")), text(v.get("missing")), proportion(v, "case_n", "percent", v.get("denominator")), "—"));
            }
        }
        h.append("</tbody></table>");
        for (Object raw : tables) { Map<String, Object> t = map(raw); h.append("<details><summary>").append(esc(text(t.get("dataset")))).append(" coding, covariates and provenance</summary>").append(render(t)).append("</details>"); }
        return h.toString();
    }
    private static void detailRow(StringBuilder h, List<String> values) { h.append("<tr>"); values.forEach(v -> h.append("<td>").append(esc(v)).append("</td>")); h.append("</tr>"); }
    static String table1Csv(List<Object> tables) {
        StringBuilder csv = new StringBuilder("dataset,cohort,field,value\n");
        for (Object raw : tables) { Map<String, Object> t = map(raw); csvFields(csv, text(t.get("dataset")), text(t.get("cohort")), "", t); }
        return csv.toString();
    }
    private static void csvFields(StringBuilder csv, String dataset, String cohort, String path, Object value) {
        if (value instanceof Map<?, ?> m) { m.forEach((key, v) -> csvFields(csv, dataset, cohort, path + "/" + key.toString().replace("~", "~0").replace("/", "~1"), v)); return; }
        if (value instanceof List<?> l) { for (int i = 0; i < l.size(); i++) csvFields(csv, dataset, cohort, path + "/" + i, l.get(i)); return; }
        csv.append(csvCell(dataset)).append(',').append(csvCell(cohort)).append(',').append(csvCell(path)).append(',').append(csvCell(text(value))).append('\n');
    }
    private static String csvCell(String s) { return "\"" + s.replace("\"", "\"\"") + "\""; }
    private static String start(String title) { return "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>" + esc(title) + "</title><style>body{font:16px/1.5 system-ui;max-width:1150px;margin:2rem auto;padding:0 1rem;color:#17212b}table{border-collapse:collapse;display:block;overflow:auto}td,th{border-bottom:1px solid #ccd3da;padding:.45rem;text-align:left;vertical-align:top}dt{font-weight:600}dd{margin-bottom:.7rem}details{margin:1rem 0}a{color:#175787}li{margin:.3rem 0}</style><body><h1>" + esc(title) + "</h1>"; }
    private static String link(Map<String, Object> record) { return "<a href=\"" + esc(text(record.get("manifest"))) + "\">Audited source manifest</a>"; }
    private static String relative(Path from, Path to) { return from.toAbsolutePath().normalize().relativize(to.toAbsolutePath().normalize()).toString().replace('\\', '/'); }
    private static String label(String value) { return value.replace('_', ' '); }
    private static String esc(String value) { return value.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("\"", "&quot;").replace("'", "&#39;"); }
    private static Map<String, Object> load(Path path) throws IOException {
        return StudyYaml.loadArtifact(path);
    }
    static String render(Object value) {
        if (value instanceof Map<?, ?> m) { StringBuilder h = new StringBuilder("<dl>"); m.forEach((k, v) -> h.append("<dt>").append(esc(label(text(k)))).append("</dt><dd>").append(v instanceof Map<?, ?> || v instanceof List<?> ? render(v) : esc(StudyReportNumbers.field(text(k), v))).append("</dd>")); return h.append("</dl>").toString(); }
        if (value instanceof List<?> l) { StringBuilder h = new StringBuilder("<ul>"); l.forEach(v -> h.append("<li>").append(render(v)).append("</li>")); return h.append("</ul>").toString(); }
        return esc(value instanceof Number ? StudyReportNumbers.format(value) : text(value));
    }

    static String expandAbbreviations(String html, Map<String, Object> supplied) {
        Map<String, Object> dictionary = new LinkedHashMap<>(Map.ofEntries(
            Map.entry("BMI", "body mass index"), Map.entry("MR", "Mendelian randomization"),
            Map.entry("LD", "linkage disequilibrium"), Map.entry("GWAS", "genome-wide association study"),
            Map.entry("irQTL", "isoform ratio quantitative trait locus"), Map.entry("QTL", "quantitative trait locus"),
            Map.entry("FDR", "false discovery rate"), Map.entry("BH", "Benjamini–Hochberg"),
            Map.entry("SD", "standard deviation"), Map.entry("SE", "standard error"),
            Map.entry("GO", "Gene Ontology"), Map.entry("LLM", "large language model"),
            Map.entry("PBMC", "peripheral blood mononuclear cells"), Map.entry("RNA", "ribonucleic acid"),
            Map.entry("DNA", "deoxyribonucleic acid")));
        dictionary.putAll(supplied);
        java.util.Set<String> seen = new java.util.HashSet<>();
        // Operate only on text nodes after <body>, never on paths, URLs or markup.
        StringBuilder result = new StringBuilder();
        boolean body = false;
        for (String part : html.split("(?=<)|(?<=>)", -1)) {
            if (part.equals("<body>")) body = true;
            if (!body || part.startsWith("<")) { result.append(part); continue; }
            String text = part;
            for (var entry : dictionary.entrySet()) {
                String shortName = entry.getKey(), longName = esc(StudyPublicEvidence.text(entry.getValue()));
                if (seen.contains(shortName) || !shortName.matches("[A-Za-z][A-Za-z0-9]{1,15}")) continue;
                var matcher = java.util.regex.Pattern.compile("(?<![A-Za-z0-9_.-])" + java.util.regex.Pattern.quote(shortName) + "(?![A-Za-z0-9_.-])").matcher(text);
                if (matcher.find()) {
                    seen.add(shortName);
                    String before = text.substring(0, matcher.start());
                    if (!before.toLowerCase(java.util.Locale.ROOT).endsWith(longName.toLowerCase(java.util.Locale.ROOT) + " ("))
                        text = before + longName + " (" + shortName + ")" + text.substring(matcher.end());
                }
            }
            result.append(text);
        }
        return result.toString();
    }
}
