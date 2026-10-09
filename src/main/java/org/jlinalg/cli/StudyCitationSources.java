/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.jlinalg.cli.StudyPublicEvidence.*;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Instant;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.jlinalg.study.StudyOperationExecutor.Context;
import org.jlinalg.study.StudyProject;
import org.yaml.snakeyaml.Yaml;
import org.yaml.snakeyaml.LoaderOptions;
import org.yaml.snakeyaml.constructor.SafeConstructor;

/** Explicit public identifier retrieval, separate from offline manuscript rendering. */
final class StudyCitationSources {
    private StudyCitationSources() { }
    static Path fetch(Path project, Path request, Path output, StudyProject.Actor actor, boolean offline, StudyPublicEvidence client) throws IOException {
        Path out = output.toAbsolutePath().normalize(), contract = project.toRealPath().resolve("resolved-study.yaml");
        if (Files.exists(out)) throw new IOException("Reference bundle exists; choose a new output");
        var study = load(contract); var spec = load(request);
        if (!Integer.valueOf(1).equals(spec.get("schema_version")) || !(spec.get("references") instanceof List<?>) || list(spec.get("references")).isEmpty()
                || list(spec.get("references")).size() > 50 || !Set.of("schema_version", "references").containsAll(spec.keySet()))
            throw new IOException("Reference request needs schema_version: 1 and 1 to 50 references");
        for (Object raw : list(spec.get("references"))) {
            var r = map(raw); String id = StudyReferenceCatalog.normalize(text(r.get("id")));
            if (!r.keySet().equals(Set.of("id", "rationale")) || text(r.get("rationale")).isBlank()
                || !(id.matches("pmid:[1-9][0-9]*") || id.matches("10\\.[0-9]{4,9}/[^\\s\\\"\\[\\]{}]+")))
                throw new IOException("Each public reference needs a DOI or pmid:ID and a local scientific rationale");
        }
        String contractHash = FollowupSupport.hash(contract), requestHash = FollowupSupport.hash(request), started = Instant.now().toString();
        if (offline) { var resources = map(study.get("resources")); resources.put("network", "offline"); study.put("resources", resources); }
        Files.createDirectories(out);
        for (String folder : List.of("metadata", "artifacts", "scripts", "logs")) Files.createDirectory(out.resolve(folder));
        Files.copy(request, out.resolve("metadata/request.yaml")); Files.copy(contract, out.resolve("metadata/resolved-study.yaml"));
        var context = new Context(project, out, Map.of("type", "reference_acquisition", "hypothesis", "manuscript"), study, Map.of(), Map.of());
        List<Object> queries = new ArrayList<>(); boolean complete = true;
        for (Object raw : list(spec.get("references"))) {
            var r = map(raw); String id = StudyReferenceCatalog.normalize(text(r.get("id"))), rationale = text(r.get("rationale"));
            String query = id.startsWith("pmid:") ? "EXT_ID:" + id.substring(5) + " AND SRC:MED" : "DOI:\"" + id + "\"";
            String url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=" + encode(query) + "&format=json&resultType=core&pageSize=5";
            var response = client.get(context, "Europe PMC", url, "https://europepmc.org/RestfulWebService", rationale);
            List<Object> articles = new ArrayList<>();
            if (response.ok()) for (Object item : list(map(map(response.data()).get("resultList")).get("result"))) {
                var article = map(item);
                if (matches(id, article)) articles.add(StudyLiterature.abstractEvidence(article, "", "", response.audit()));
            }
            // Crossref supplies publisher-deposited bibliographic metadata for methods not indexed by PubMed.
            if (articles.isEmpty() && response.ok() && id.startsWith("10.")) {
                var crossref = client.get(context, "Crossref", "https://api.crossref.org/works/" + encode(id), "https://www.crossref.org/documentation/retrieve-metadata/rest-api/", rationale);
                if (crossref.ok()) {
                    var record = crossrefArticle(map(map(crossref.data()).get("message")));
                    if (matches(id, record)) articles.add(record);
                }
            }
            complete &= !articles.isEmpty();
            queries.add(Map.of("requested_reference", id, "query", query, "scientific_rationale", rationale,
                "status", articles.isEmpty() ? "unresolved" : "candidate_matches", "articles", articles, "audit", response.audit(),
                "search_strategy", "exact_identifier", "source_manifest", "metadata/manifest.yaml"));
        }
        write(out.resolve("artifacts/queries.yaml"), Map.of("queries", queries));
        write(out.resolve("artifacts/reference-catalog.yaml"), StudyReferenceCatalog.build(queries));
        Files.writeString(out.resolve("artifacts/reference-catalog.html"), StudyReferenceCatalog.html(StudyReferenceCatalog.build(queries)));
        write(out.resolve("logs/requests.yaml"), Map.of("requests", client.requests()));
        write(out.resolve("metadata/software-versions.yaml"), Map.of("jlinalg", JLinAlgCli.version(), "java", System.getProperty("java.runtime.version"), "actor", actor.asMap()));
        Files.writeString(out.resolve("metadata/README.md"), "# Native public reference acquisition\n\nExact public DOI/PMID queries only; no manuscript text, local rationale, study files or patient rows transmitted. Europe PMC core records and publisher-deposited Crossref fallback metadata are archived with requests, timing and response hashes. Retrieval is not verification of claim support. Unknown identifiers and blocked queries remain unresolved. This command respects the source study's transmission/download policy and --offline. Import the manifest with its SHA-256 as reference_bundles in a new offline manuscript revision. Source contract: " + contract + "\n");
        Files.writeString(out.resolve("SESSION_HISTORY.md"), started + " - JLinAlg public reference acquisition; actor " + actor.asMap() + "; outcome " + (complete ? "complete" : "needs_input") + ". No analyses or source contract modified.\n");
        if (!requestHash.equals(FollowupSupport.hash(request)) || !contractHash.equals(FollowupSupport.hash(contract))) throw new IOException("Reference request or study policy changed during acquisition");
        Map<String, Object> hashes = new LinkedHashMap<>();
        try (var paths = Files.walk(out)) { for (Path p : paths.filter(Files::isRegularFile).sorted().toList()) hashes.put(out.relativize(p).toString().replace('\\', '/'), FollowupSupport.hash(p)); }
        write(out.resolve("metadata/manifest.yaml"), Map.of("schema_version", 1, "kind", "reference_acquisition", "state", complete ? "complete" : "needs_input",
            "started_at", started, "completed_at", Instant.now().toString(), "request_sha256", requestHash, "study_sha256", contractHash, "files_sha256", hashes));
        return out;
    }
    private static boolean matches(String id, Map<String, Object> a) {
        return StudyReferenceCatalog.identifiers(a).contains(id.startsWith("10.") ? "doi:" + id : id);
    }
    static Map<String, Object> crossrefArticle(Map<String, Object> message) {
        Map<String, Object> a = new LinkedHashMap<>();
        a.put("doi", text(message.get("DOI"))); a.put("title", first(message.get("title")));
        a.put("journalInfo", Map.of("volume", text(message.get("volume")), "journal", Map.of("title", first(message.get("container-title")))));
        a.put("pageInfo", message.getOrDefault("page", message.getOrDefault("article-number", "")));
        List<?> dates = list(map(message.getOrDefault("published-print", message.get("published"))).get("date-parts"));
        a.put("pubYear", dates.isEmpty() || list(dates.get(0)).isEmpty() ? "" : text(list(dates.get(0)).get(0)));
        List<Object> authors = new ArrayList<>();
        for (Object raw : list(message.get("author"))) {
            var author = map(raw); String given = text(author.get("given"));
            String initials = java.util.Arrays.stream(given.split("[\\s.-]+" )).filter(s -> !s.isBlank()).map(s -> s.substring(0, 1)).collect(java.util.stream.Collectors.joining());
            authors.add(Map.of("lastName", text(author.get("family")), "initials", initials, "fullName", text(author.getOrDefault("name", author.get("family")))));
        }
        a.put("authorList", Map.of("author", authors)); a.put("retrieval_provider", "Crossref publisher-deposited metadata"); return a;
    }
    private static String first(Object value) { return list(value).isEmpty() ? "" : text(list(value).get(0)); }
    static List<Object> imports(Map<String, Object> revision, Path request, List<Path> manifests) throws IOException {
        if (revision.containsKey("reference_bundles") && !(revision.get("reference_bundles") instanceof List<?>)) throw new IOException("reference_bundles must be a list");
        List<Object> queries = new ArrayList<>(); List<Object> normalized = new ArrayList<>();
        for (Object raw : list(revision.get("reference_bundles"))) {
            var ref = map(raw);
            if (!ref.keySet().equals(Set.of("manifest", "sha256"))) throw new IOException("reference_bundles entries require manifest and sha256");
            Path manifest = request.toAbsolutePath().getParent().resolve(text(ref.get("manifest"))).toRealPath();
            if (!text(ref.get("sha256")).equals(FollowupSupport.hash(manifest))) throw new IOException("Reference bundle manifest hash mismatch");
            verify(manifest); manifests.add(manifest);
            normalized.add(Map.of("manifest", manifest.toString(), "sha256", ref.get("sha256")));
            var loaded = load(manifest.getParent().getParent().resolve("artifacts/queries.yaml"));
            for (Object q : list(loaded.get("queries"))) { var item = map(q); item.put("source_manifest", manifest.toString()); queries.add(item); }
        }
        if (revision.containsKey("reference_bundles")) revision.put("reference_bundles", normalized);
        return queries;
    }
    static void verify(Path manifest) throws IOException {
        var m = load(manifest); Path root = manifest.getParent().getParent().toRealPath();
        if (!"reference_acquisition".equals(m.get("kind"))) throw new IOException("Not a native reference acquisition bundle");
        var hashes = map(m.get("files_sha256"));
        if (!hashes.containsKey("artifacts/queries.yaml") || !hashes.containsKey("logs/requests.yaml")) throw new IOException("Incomplete reference bundle inventory");
        for (var e : hashes.entrySet()) {
            Path file = root.resolve(e.getKey()).toRealPath();
            if (!file.startsWith(root) || !e.getValue().equals(FollowupSupport.hash(file))) throw new IOException("Reference bundle artifact hash mismatch: " + e.getKey());
        }
    }
    private static Map<String, Object> load(Path file) throws IOException {
        var options = new LoaderOptions(); options.setAllowDuplicateKeys(false); options.setMaxAliasesForCollections(0);
        return map(new Yaml(new SafeConstructor(options)).load(Files.readString(file)));
    }
}
