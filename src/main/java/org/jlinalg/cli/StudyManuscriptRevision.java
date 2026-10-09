/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.time.Instant;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.jlinalg.study.StudyEvidenceIndex;
import org.jlinalg.study.StudyProject;
import org.yaml.snakeyaml.LoaderOptions;
import org.yaml.snakeyaml.Yaml;
import org.yaml.snakeyaml.constructor.SafeConstructor;

/** Offline editorial revisions derived from verified evidence, never new analyses. */
final class StudyManuscriptRevision {
    private StudyManuscriptRevision() { }
    private static final Set<String> SECTIONS = Set.of("introduction", "methods", "results", "discussion", "conclusion");
    private static final Pattern TOKEN = Pattern.compile("\\{\\{([A-Za-z0-9._-]+)#(/[^{}]*)}}" );
    private static final String LIMIT = "Unreviewed manuscript revision. Evidence integrity is checked; scientific interpretation, authorship, ethics, novelty and publication readiness are not certified. No network or external LLM call is made by this command.";

    static Path create(Path project, Path request, Path output, StudyProject.Actor actor) throws IOException {
        Path root = project.toRealPath(), out = output.toAbsolutePath().normalize();
        if (Files.exists(out, LinkOption.NOFOLLOW_LINKS)) throw new IOException("Revision output exists; choose a new bundle");
        if (out.startsWith(root)) throw new IOException("Revision bundle must be outside the frozen source project");
        String requestHash = FollowupSupport.hash(request);
        Map<String, Object> edit = load(request);
        keys(edit, Set.of("schema_version", "id", "title", "rationale", "previous_revision", "sections", "remove_sections", "comments", "abbreviations", "participant_project", "discussion_selection_note", "paragraph_order", "citation_policy", "reference_bundles"), "revision");
        if (!Integer.valueOf(1).equals(edit.get("schema_version"))) throw new IOException("Revision schema_version must be 1");
        identifier(required(edit, "id")); required(edit, "title"); required(edit, "rationale");
        String start = Instant.now().toString();
        Path ledger = root.resolve("ledger/events.jsonl");
        if (!Files.isRegularFile(ledger)) throw new IOException("A resolved, executed and ledger-anchored source project is required");
        String ledgerHash = FollowupSupport.hash(ledger);
        StudyProject.verify(root);
        Map<String, Object> index = StudyEvidenceIndex.snapshot(root);
        Map<String, Object> study = load(root.resolve("resolved-study.yaml"));
        Map<String, Map<String, Object>> evidence = new LinkedHashMap<>();
        for (Object value : list(index.get("entries"))) {
            Map<String, Object> e = map(value);
            if (!Boolean.TRUE.equals(e.get("current"))) continue;
            String key = e.get("hypothesis") + "." + e.get("operation_type");
            if (evidence.putIfAbsent(key, e) != null) throw new IOException("Ambiguous current operation: " + key);
        }
        if (evidence.isEmpty()) throw new IOException("No current executed evidence available");
        Map<String, Object> previous = new LinkedHashMap<>();
        Map<String, Object> attribution = new LinkedHashMap<>();
        Path predecessor = null;
        if (edit.containsKey("previous_revision")) {
            Map<String, Object> ref = map(edit.get("previous_revision"));
            keys(ref, Set.of("manifest", "sha256"), "previous_revision");
            predecessor = request.toAbsolutePath().getParent().resolve(required(ref, "manifest")).normalize();
            if (!FollowupSupport.hash(predecessor).equals(required(ref, "sha256"))) throw new IOException("Predecessor manifest hash mismatch");
            verifyBundle(predecessor.getParent().getParent());
            previous = load(predecessor.getParent().resolve("resolved-revision.yaml"));
            attribution.putAll(load(predecessor.getParent().resolve("section-attribution.yaml")));
        }
        Map<String, Object> merged = new LinkedHashMap<>(previous);
        merged.putAll(edit);
        // Explicit, hash-pinned native Table 1 refresh, without replacing association evidence.
        Path participantRoot = null;
        String participantLedgerHash = null;
        Set<String> participantSources = new LinkedHashSet<>();
        if (merged.containsKey("participant_project")) {
            Map<String, Object> ref = map(merged.get("participant_project"));
            keys(ref, Set.of("path", "ledger_sha256", "rationale"), "participant_project");
            required(ref, "rationale");
            participantRoot = request.toAbsolutePath().getParent().resolve(required(ref, "path")).toRealPath();
            participantLedgerHash = required(ref, "ledger_sha256");
            if (!FollowupSupport.hash(participantRoot.resolve("ledger/events.jsonl")).equals(participantLedgerHash)) throw new IOException("Participant project ledger hash mismatch");
            StudyProject.verify(participantRoot);
            for (Object raw : list(StudyEvidenceIndex.snapshot(participantRoot).get("entries"))) {
                Map<String, Object> entry = map(raw);
                if (!Boolean.TRUE.equals(entry.get("current")) || !"participant_characteristics".equals(entry.get("operation_type"))) continue;
                String key = entry.get("hypothesis") + ".participant_characteristics";
                if (!evidence.containsKey(key)) throw new IOException("Participant refresh hypothesis does not match the source project: " + key);
                if (!participantSources.add(key)) throw new IOException("Ambiguous participant refresh source: " + key);
                if (!participantIdentities(entry).equals(participantIdentities(evidence.get(key))))
                    throw new IOException("Participant refresh datasets/cohorts differ; amend the source study first: " + key);
                entry.put("source_project", participantRoot.toString());
                evidence.put(key, entry);
            }
            if (participantSources.isEmpty()) throw new IOException("Participant project contains no current participant summaries");
        }
        List<Object> changes = new ArrayList<>();
        Map<String, Object> paragraphs = new LinkedHashMap<>();
        for (Object p : list(previous.get("sections"))) paragraphs.put(required(map(p), "id"), p);
        Set<String> submitted = new LinkedHashSet<>();
        for (String field : List.of("sections", "remove_sections", "comments"))
            if (edit.containsKey(field) && !(edit.get(field) instanceof List<?>)) throw new IOException(field + " must be a list");
        for (Object value : list(edit.get("sections"))) {
            Map<String, Object> p = map(value); String id = required(p, "id"); identifier(id);
            if (!submitted.add(id)) throw new IOException("Duplicate section ID: " + id);
            keys(p, Set.of("id", "section", "title", "text", "sources", "rationale", "references", "introduction_role", "methods_role", "cohort_ids", "parameters", "parameter_justifications", "uncommon_parameters", "results_role", "supporting_evidence", "discussion_role", "impact_rank", "impact_reason", "result_ids", "citation_basis"), "paragraph");
            if (!SECTIONS.contains(required(p, "section"))) throw new IOException("Unknown manuscript section");
            required(p, "text"); required(p, "rationale");
            Object old = paragraphs.put(id, p);
            attribution.put(id, Map.of("actor", actor.asMap(), "revision", edit.get("id"), "recorded_at", start));
            changes.add(Map.of("id", id, "action", old == null ? "added" : "replaced", "before", old == null ? Map.of() : old, "after", p, "actor", actor.asMap()));
        }
        for (Object value : list(edit.get("remove_sections"))) {
            Map<String, Object> removal = map(value); keys(removal, Set.of("id", "rationale"), "removal");
            String id = required(removal, "id"); required(removal, "rationale");
            if (submitted.contains(id)) throw new IOException("Cannot replace and remove the same section: " + id);
            Object old = paragraphs.remove(id);
            attribution.remove(id);
            if (old == null) throw new IOException("Unknown section to remove: " + id);
            changes.add(Map.of("id", id, "action", "removed", "before", old, "rationale", removal.get("rationale"), "actor", actor.asMap()));
        }
        // Explicit reordering changes placement, not authorship. Old order requests are not replayed on later additions.
        merged.remove("paragraph_order");
        if (edit.containsKey("paragraph_order")) {
            if (!(edit.get("paragraph_order") instanceof List<?>)) throw new IOException("paragraph_order must list every retained paragraph ID exactly once");
            Map<String, Object> ordered = new LinkedHashMap<>();
            for (Object id : list(edit.get("paragraph_order"))) {
                if (!(id instanceof String) || !paragraphs.containsKey(id) || ordered.putIfAbsent(id.toString(), paragraphs.get(id)) != null)
                    throw new IOException("paragraph_order contains an unknown or repeated paragraph ID");
            }
            if (ordered.size() != paragraphs.size()) throw new IOException("paragraph_order must include every retained paragraph ID");
            changes.add(Map.of("action", "paragraphs_reordered", "before", new ArrayList<>(paragraphs.keySet()), "after", new ArrayList<>(ordered.keySet()), "actor", actor.asMap()));
            paragraphs.clear(); paragraphs.putAll(ordered); merged.put("paragraph_order", edit.get("paragraph_order"));
        }
        merged.put("sections", new ArrayList<>(paragraphs.values()));
        Map<String, Object> comments = new LinkedHashMap<>();
        for (Object c : list(previous.get("comments"))) comments.put(required(map(c), "id"), c);
        Set<String> commentIds = new LinkedHashSet<>();
        for (Object c : list(edit.get("comments"))) {
            String id = required(map(c), "id");
            if (!commentIds.add(id)) throw new IOException("Duplicate comment ID: " + id);
            Object before = comments.put(id, c);
            changes.add(Map.of("id", id, "action", "comment_updated", "before", before == null ? Map.of() : before, "after", c, "actor", actor.asMap()));
        }
        merged.put("comments", new ArrayList<>(comments.values()));
        Map<String, Object> abbreviations = new LinkedHashMap<>(map(previous.get("abbreviations")));
        abbreviations.putAll(map(edit.get("abbreviations"))); merged.put("abbreviations", abbreviations);
        List<Object> claims = new ArrayList<>();
        List<Path> referenceManifests = new ArrayList<>();
        Map<String, Object> catalog = StudyReferenceCatalog.fromEvidence(evidence, StudyCitationSources.imports(merged, request, referenceManifests));
        Map<String, Object> citationReview = StudyCitations.assess(merged, catalog);
        Map<String, Object> introductionReview = new LinkedHashMap<>(StudyIntroduction.assess(list(merged.get("sections"))));
        introductionReview.put("literature_search_coverage", catalog.get("search_coverage"));
        var methodsInventory = StudyMethods.inventory(study, evidence);
        var methodsReview = StudyMethods.assess(list(merged.get("sections")), methodsInventory);
        var narrativeReview = StudyResultsDiscussion.assess(merged, evidence);
        for (Object value : paragraphs.values()) {
            Map<String, Object> p = map(value);
            if (!(p.get("sources") instanceof List<?>)) throw new IOException("Paragraph sources must be a list (empty allowed for editorial context)");
            for (Object source : list(p.get("sources"))) if (!evidence.containsKey(source.toString())) throw new IOException("Unknown evidence source: " + source);
            if (p.containsKey("references") && !(p.get("references") instanceof List<?>)) throw new IOException("references must be a list of cataloged DOIs, identifiers or reference IDs");
            for (Object reference : list(p.get("references"))) if (StudyReferenceCatalog.lookup(catalog, reference.toString()).isEmpty()) throw new IOException("Reference is not present in current retrieved evidence: " + reference);
            if ("results".equals(p.get("section"))) {
                if (list(p.get("sources")).stream().noneMatch(s -> "complete".equals(evidence.get(s.toString()).get("state"))))
                    throw new IOException("A main Results point requires at least one completed evidence source");
            }
            String rendered = expand(required(p, "text"), evidence, claims, required(p, "id"));
            claims.add(Map.of("section_id", p.get("id"), "text", rendered, "sources", list(p.get("sources")),
                "assessment", "Editorial text requires scientific review; only template bindings are mechanically checked", "attribution", attribution.get(p.get("id"))));
            for (Object raw : list(p.get("supporting_evidence"))) {
                var support = map(raw);
                String relationship = expand(required(support, "relationship"), evidence, claims, required(p, "id"));
                claims.add(Map.of("section_id", p.get("id"), "supporting_source", support.get("source"), "relationship", relationship,
                    "assessment", "Attributed interpretation; source identity is checked, biological relevance and independence require scientific review", "attribution", attribution.get(p.get("id"))));
            }
            if (p.containsKey("impact_reason")) claims.add(Map.of("section_id", p.get("id"), "impact_rank", p.get("impact_rank"),
                "impact_reason", expand(required(p, "impact_reason"), evidence, claims, required(p, "id")), "result_ids", p.get("result_ids"), "attribution", attribution.get(p.get("id"))));
        }
        for (Object value : list(merged.get("comments"))) {
            Map<String, Object> c = map(value); keys(c, Set.of("id", "comment", "response", "section", "status"), "comment");
            required(c, "id"); required(c, "comment");
            if (!Set.of("addressed", "pending", "deferred").contains(required(c, "status"))) throw new IOException("Invalid comment status");
            if (c.containsKey("section") && !paragraphs.containsKey(c.get("section"))) throw new IOException("Comment references unknown section");
            if ("addressed".equals(c.get("status"))) required(c, "response");
        }
        // All semantic validation precedes publication. The unique output is never overwritten.
        String editorialPayload = StudyAssistance.json(StudyEditorialReview.payload(merged, evidence, methodsInventory, catalog));
        Files.createDirectories(out.getParent()); Files.createDirectory(out);
        for (String directory : List.of("artifacts", "scripts", "metadata", "logs", "previews")) Files.createDirectory(out.resolve(directory));
        Files.copy(request, out.resolve("metadata/revision-request.yaml"));
        writeYaml(out.resolve("metadata/resolved-revision.yaml"), merged);
        writeYaml(out.resolve("metadata/section-attribution.yaml"), attribution);
        writeYaml(out.resolve("metadata/resolved-study.yaml"), study);
        writeYaml(out.resolve("metadata/change-log.yaml"), Map.of("changes", changes, "comments", list(merged.get("comments")),
            "note", "Omitted paragraphs are retained from the predecessor; replacements and removals require explicit IDs and reasons. Evidence-bound values are refreshed from the current source snapshot."));
        Files.copy(ledger, out.resolve("metadata/source-ledger.jsonl"));
        if (predecessor != null) Files.copy(predecessor, out.resolve("metadata/predecessor-manifest.yaml"));
        List<Object> inputs = new ArrayList<>();
        int referenceSequence = 0;
        for (Path manifest : referenceManifests) {
            inputs.add(input(manifest, FollowupSupport.hash(manifest), "audited public reference acquisition"));
            Path copy = out.resolve("metadata/reference-bundles/" + (++referenceSequence)); Files.createDirectories(copy);
            Path referenceRoot = manifest.getParent().getParent();
            try (var paths = Files.walk(referenceRoot)) {
                for (Path file : paths.filter(Files::isRegularFile).toList()) {
                    if (file.toRealPath().startsWith(referenceRoot.toRealPath())) {
                        Path target = copy.resolve(referenceRoot.relativize(file)); Files.createDirectories(target.getParent()); Files.copy(file, target);
                    } else throw new IOException("Reference bundle path escapes its root");
                }
            }
        }
        if (participantRoot != null) {
            inputs.add(input(participantRoot.resolve("ledger/events.jsonl"), participantLedgerHash, "explicit participant-summary refresh ledger"));
            Path participantStudy = participantRoot.resolve("resolved-study.yaml");
            inputs.add(input(participantStudy, FollowupSupport.hash(participantStudy), "participant-summary contract"));
            Files.copy(participantStudy, out.resolve("metadata/participant-study.yaml"));
            Files.copy(participantRoot.resolve("ledger/events.jsonl"), out.resolve("metadata/participant-ledger.jsonl"));
        }
        inputs.add(input(request.toAbsolutePath(), requestHash, "revision request"));
        inputs.add(input(ledger, ledgerHash, "source ledger"));
        inputs.add(input(root.resolve("resolved-study.yaml"), FollowupSupport.hash(root.resolve("resolved-study.yaml")), "resolved study"));
        for (Map.Entry<String, Map<String, Object>> item : evidence.entrySet()) {
            Path manifest = contained(participantSources.contains(item.getKey()) ? participantRoot : root, item.getValue().get("manifest").toString());
            String hash = FollowupSupport.hash(manifest);
            item.getValue().put("manifest_sha256", hash);
            Path copy = out.resolve("metadata/sources/" + item.getKey() + ".yaml");
            Files.createDirectories(copy.getParent()); Files.copy(manifest, copy);
            inputs.add(input(manifest, hash, "verified source operation"));
            Path oldDraft = manifest.getParent().resolve("publication-draft.html");
            if ("report".equals(item.getValue().get("operation_type")) && Files.isRegularFile(oldDraft)) {
                String oldHash = FollowupSupport.hash(oldDraft);
                if (!oldHash.equals(map(load(manifest).get("output_sha256")).get("publication-draft.html"))) throw new IOException("Previous draft is not ledger-anchored");
                Files.copy(oldDraft, out.resolve("metadata/sources/" + item.getKey() + "-previous-draft.txt"));
                inputs.add(input(oldDraft, oldHash, "previous generated draft preserved as inert text"));
            }
        }
        writeYaml(out.resolve("artifacts/evidence.yaml"), evidence);
        writeYaml(out.resolve("artifacts/publication-claims.yaml"), Map.of("claims", claims, "limitations", LIMIT));
        write(out.resolve("artifacts/editorial-payload.json"), editorialPayload);
        writeYaml(out.resolve("artifacts/citation-review.yaml"), citationReview);
        writeYaml(out.resolve("artifacts/methods-inventory.yaml"), methodsInventory);
        writeYaml(out.resolve("artifacts/study-context-review.yaml"), StudyContextReview.inventory(study, null));
        writeYaml(out.resolve("artifacts/methods-review.yaml"), methodsReview);
        writeYaml(out.resolve("artifacts/methods-template.yaml"), StudyMethods.template(required(edit, "title"), methodsInventory, catalog));
        writeYaml(out.resolve("artifacts/results-discussion-inventory.yaml"), StudyResultsDiscussion.inventory(evidence));
        writeYaml(out.resolve("artifacts/results-discussion-review.yaml"), narrativeReview);
        writeYaml(out.resolve("artifacts/results-discussion-template.yaml"), StudyResultsDiscussion.template(required(edit, "title"), evidence));
        String manuscript = manuscript(merged, evidence, study, catalog);
        write(out.resolve("artifacts/manuscript.html"), StudyPublication.expandAbbreviations(manuscript, abbreviations));
        write(out.resolve("artifacts/supplement.html"), supplement(evidence));
        write(out.resolve("artifacts/participant-characteristics.csv"), StudyPublication.table1Csv(participantTables(evidence)));
        writeYaml(out.resolve("metadata/reporting-policy.yaml"), Map.of("decimal_places", 2, "scientific_notation", "lowercase e; two mantissa decimals; nonzero magnitudes below 0.01 or at least 1000000",
            "integer_counts", "unchanged", "machine_readable_precision", "unchanged; raw values and bindings retained", "participant_refresh", merged.getOrDefault("participant_project", Map.of())));
        write(out.resolve("artifacts/literature-review.html"), literature(evidence));
        writeYaml(out.resolve("artifacts/reference-catalog.yaml"), catalog);
        write(out.resolve("artifacts/reference-catalog.html"), StudyReferenceCatalog.html(catalog));
        writeYaml(out.resolve("artifacts/introduction-review.yaml"), introductionReview);
        writeYaml(out.resolve("artifacts/introduction-template.yaml"), StudyIntroduction.template(required(edit, "title"),
            list(study.get("hypotheses")).stream().map(StudyManuscriptRevision::map).filter(h -> !"retired".equals(h.get("status"))).map(h -> String.valueOf(h.get("description"))).toList()));
        write(out.resolve("artifacts/index.html"), start("Manuscript revision evidence index")
            + "<p><a href=\"results-discussion-template.yaml\">Results and Discussion template</a> · <a href=\"results-discussion-review.yaml\">Narrative review needs</a> · <a href=\"results-discussion-inventory.yaml\">Findings and contextual evidence</a></p>"
            + "<p><a href=\"methods-template.yaml\">Methods writing template</a> · <a href=\"methods-review.yaml\">Methods review needs</a> · <a href=\"methods-inventory.yaml\">Methods provenance</a></p>"
            + "<p>" + esc(LIMIT) + "</p><ul><li><a href=\"manuscript.html\">Revised manuscript</a></li><li><a href=\"supplement.html\">All evidence layers and limitations</a></li><li><a href=\"literature-review.html\">Literature review candidates</a></li><li><a href=\"reference-catalog.html\">Reference catalog and search coverage</a></li><li><a href=\"introduction-review.yaml\">Introduction structure and research needs</a></li><li><a href=\"introduction-template.yaml\">Introduction paragraph template</a></li><li><a href=\"publication-claims.yaml\">Claim and value provenance</a></li><li><a href=\"../metadata/change-log.yaml\">Editorial changes and responses</a></li><li><a href=\"../metadata/manifest.yaml\">Integrity manifest</a></li></ul>" + end());
        Map<String, Object> software = new LinkedHashMap<>();
        software.put("jlinalg", JLinAlgCli.version()); software.put("java", System.getProperty("java.runtime.version"));
        software.put("os", System.getProperty("os.name") + " " + System.getProperty("os.version")); software.put("actor", actor.asMap());
        software.put("network_requests", 0); software.put("external_llm_calls", 0);
        try {
            Path runtime = Path.of(StudyManuscriptRevision.class.getProtectionDomain().getCodeSource().getLocation().toURI());
            if (Files.isRegularFile(runtime) && runtime.toString().endsWith(".jar")) {
                Files.copy(runtime, out.resolve("scripts/jlinalg.jar")); software.put("jar_sha256", FollowupSupport.hash(runtime));
            } else software.put("jar_sha256", "unknown: class-directory execution, no packaged JAR");
        } catch (java.net.URISyntaxException ex) { throw new IOException(ex); }
        writeYaml(out.resolve("metadata/software-versions.yaml"), software);
        writeYaml(out.resolve("metadata/inputs.yaml"), Map.of("inputs", inputs));
        write(out.resolve("metadata/README.md"), "# Native JLinAlg manuscript revision\n\n" + LIMIT + "\n\n"
            + "Created: " + start + "\n\nSource project: " + root + "\n\nPurpose: " + required(edit, "rationale")
            + "\n\nThis is a reporting-only evidence snapshot, not a study amendment or new analysis. Source metadata and participant-population uncertainties cannot be fixed by editing prose. Use an approved study amendment and an analysis rerun for scientific changes."
            + "\n\nIntroduction development: artifacts/introduction-template.yaml supplies the ordered paragraph roles; artifacts/introduction-review.yaml records missing structure/citations and search coverage. artifacts/reference-catalog.html and .yaml retain deduplicated retrieved sources and query provenance. Catalog inclusion and structural completeness do not verify source support, sufficient literature coverage or novelty. Existing untagged Introduction prose is preserved for explicit review, not silently rewritten."
            + "\n\nMethods development: artifacts/methods-template.yaml follows cohort descriptions with the performed analyses. methods-inventory.yaml retains declared external models and recorded execution parameters; methods-review.yaml lists missing citations, parameters, justifications and cohort coverage. Historical gaps are not filled using current software defaults. Editorial parameters explicitly supplied in a revision are part of the hash-approved disclosure payload and require sensitivity review."
            + "\n\nResults and Discussion: results-discussion-template.yaml follows Methods order after native Table 1. supporting_evidence links bioinformatics context beside its finding with an attributed relationship. Discussion prioritizes three distinct findings with impact reasons and Results paragraph links, then future potential, next steps, limitations and the overall conclusion. Fewer supported findings require a discussion_selection_note, not invented findings. results-discussion-review.yaml flags unclassified older prose and incomplete structure; no impact ranking is inferred from analysis type or P-values. These editorial fields enter the hash-approved payload; automatic native inventories do not."
            + "\n\nReproduction: java -jar scripts/jlinalg.jar study manuscript metadata/revision-request.yaml --project \"" + root + "\" --out NEW_BUNDLE --actor-type " + actor.type() + " --actor-id \"" + actor.id() + "\""
            + (actor.model() == null ? "" : " --model \"" + actor.model() + "\"")
            + "\n\nSee software-versions.yaml for exact versions and input/output manifests for SHA-256. All local evidence layers are retained in artifacts/evidence.yaml; original large result tables and public response caches remain in the predecessor project. Keep both bundles together. No patient-level files were copied or read by this revision command. Text fields may still contain sensitive information and must be reviewed before sharing."
            + "\n\nHTML editors must not alter this frozen bundle. Revise section IDs in a new revision request and name this manifest plus its SHA-256 as previous_revision to preserve and track prose changes. Manual edits outside that workflow require a separate documented audit trail.\n");
        write(out.resolve("SESSION_HISTORY.md"), "# Manuscript revision notebook\n\n" + start + " — Started reporting-only revision " + edit.get("id")
            + ". Actor: " + actor.asMap() + ". Native JLinAlg " + JLinAlgCli.version() + " verified the source ledger and generated the manuscript, supplement, literature candidates and claim map. Editorial inputs and decisions are preserved in metadata. No analysis rerun or network request.\n");
        // Fail if either source authority or request changed while building the report.
        StudyProject.verify(root);
        if (participantRoot != null) StudyProject.verify(participantRoot);
        if (!ledgerHash.equals(FollowupSupport.hash(ledger)) || !requestHash.equals(FollowupSupport.hash(request))) throw new IOException("Source or revision request changed during generation; output remains incomplete");
        for (Path manifest : referenceManifests) StudyCitationSources.verify(manifest);
        for (Object raw : inputs) {
            Map<String, Object> i = map(raw);
            if (!i.get("sha256").equals(FollowupSupport.hash(Path.of(i.get("path").toString())))) throw new IOException("Source changed during generation");
        }
        writeYaml(out.resolve("logs/validation.yaml"), Map.of("source_project_integrity", "passed before and after generation",
            "source_unchanged", true, "numeric_bindings", "resolved against verified native aggregate evidence", "scientific_review", "required", "visual_qa", "not performed by this command", "completed_at", Instant.now().toString()));
        Map<String, Object> hashes = new LinkedHashMap<>();
        try (var paths = Files.walk(out)) { for (Path p : paths.filter(Files::isRegularFile).sorted().toList()) hashes.put(relative(out, p), FollowupSupport.hash(p)); }
        writeYaml(out.resolve("metadata/manifest.yaml"), Map.of("schema_version", 1, "id", edit.get("id"), "state", "complete",
            "started_at", start, "completed_at", Instant.now().toString(), "actor", actor.asMap(), "source_ledger_sha256", ledgerHash, "files_sha256", hashes));
        verifyBundle(out);
        return out;
    }

    static void verifyBundle(Path bundle) throws IOException {
        Map<String, Object> manifest = load(bundle.resolve("metadata/manifest.yaml"));
        if (!"complete".equals(manifest.get("state"))) throw new IOException("Predecessor revision is incomplete");
        Map<String, Object> hashes = map(manifest.get("files_sha256"));
        if ("editorial_review".equals(manifest.get("bundle_type"))) {
            if (!hashes.containsKey("metadata/approved-payload.json") || !hashes.containsKey("artifacts/decision.yaml")) throw new IOException("Incomplete editorial review inventory");
        } else if (!hashes.containsKey("metadata/resolved-revision.yaml") || !hashes.containsKey("artifacts/manuscript.html")) throw new IOException("Incomplete revision inventory");
        for (var entry : hashes.entrySet()) {
            Path p = contained(bundle, entry.getKey());
            if (!entry.getValue().equals(FollowupSupport.hash(p))) throw new IOException("Revision artifact hash mismatch: " + entry.getKey());
        }
    }

    private static String manuscript(Map<String, Object> revision, Map<String, Map<String, Object>> evidence,
            Map<String, Object> study, Map<String, Object> catalog) throws IOException {
        StringBuilder h = new StringBuilder(start(required(revision, "title")));
        StudyCitations citations = new StudyCitations(catalog);
        h.append("<nav><a href=\"#references\">References</a> · <a href=\"citation-review.yaml\">Citation review</a> · <a href=\"study-context-review.yaml\">Study context and candidate explanations</a></nav>");
        if ("needs_input".equals(StudyCitations.assess(revision, catalog).get("status"))) h.append("<p class=\"notice\">Citation placement or classification needs review; see citation-review.yaml. Journal-derived statements require citations throughout the paper.</p>");
        h.append("<p class=\"notice\">").append(esc(LIMIT)).append("</p><nav><a href=\"supplement.html\">Supplement</a> · <a href=\"reference-catalog.html\">Reference catalog</a> · <a href=\"introduction-review.yaml\">Introduction review</a> · <a href=\"publication-claims.yaml\">Claim audit</a> · <a href=\"../metadata/change-log.yaml\">Changes</a></nav>");
        if ("needs_input".equals(StudyIntroduction.assess(list(revision.get("sections"))).get("status")))
            h.append("<p class=\"notice\">Introduction development remains incomplete. See the introduction review for missing paragraphs or literature inputs; existing prose has been preserved.</p>");
        h.append("<p><a href=\"results-discussion-review.yaml\">Results and Discussion review</a> · <a href=\"results-discussion-template.yaml\">Writing template</a> · <a href=\"results-discussion-inventory.yaml\">Findings and context inventory</a></p>");
        if ("needs_input".equals(StudyResultsDiscussion.assess(revision, evidence).get("status")))
            h.append("<p class=\"notice\">Results and Discussion need structural or scientific inputs; see the narrative review. Older prose is preserved, not silently rewritten.</p>");
        Set<String> previewed = new LinkedHashSet<>();
        for (String section : List.of("introduction", "methods", "results", "discussion", "conclusion")) {
            if (section.equals("conclusion") && list(revision.get("sections")).stream().map(StudyManuscriptRevision::map).noneMatch(p -> "conclusion".equals(p.get("section")))) continue;
            h.append("<h2>").append(label(section)).append("</h2>");
            if (section.equals("methods")) {
                if (list(revision.get("sections")).stream().map(StudyManuscriptRevision::map).noneMatch(p -> "cohorts".equals(p.get("methods_role")))) {
                    var template = StudyMethods.template("Methods", StudyMethods.inventory(study, evidence), catalog);
                    h.append(StudyMethods.scaffold(Map.of("sections", list(template.get("sections")).stream().filter(p -> "cohorts".equals(map(p).get("methods_role"))).toList()), citations));
                }
            }
            if (section.equals("results")) {
                h.append("<h3>Table 1. Participant characteristics</h3><p>SD means sample standard deviation. <a href=\"supplement.html#supplemental-table-1\">Supplemental Table 1</a> retains the detailed summaries.</p>");
                h.append(StudyPublication.table1(participantTables(evidence)));
            }
            for (Object raw : list(revision.get("sections"))) {
                Map<String, Object> p = map(raw); if (!section.equals(p.get("section"))) continue;
                boolean intro = section.equals("introduction");
                h.append("<section id=\"").append(esc(p.get("id"))).append("\">");
                if (!intro && !p.containsKey("discussion_role")) h.append("<h3>").append(esc(p.getOrDefault("title", label(section)))).append("</h3>");
                String text = expand(required(p, "text"), evidence, new ArrayList<>(), required(p, "id"), true);
                String[] proseParagraphs = text.split("\\R\\s*\\R");
                for (int i = 0; i < proseParagraphs.length; i++) h.append("<p>").append(citations.render(proseParagraphs[i], i == proseParagraphs.length - 1 ? list(p.get("references")).stream().filter(id -> !StudyCitations.identifiers(Map.of("text", text), "text").contains(id.toString())).toList() : List.of())).append("</p>");
                if (!intro && !section.equals("methods")) {
                    h.append("<details><summary>Editorial rationale and evidence</summary><p class=\"source\">").append(esc(p.get("rationale"))).append("</p>");
                    if (p.containsKey("impact_reason")) {
                        h.append("<p>Impact rank ").append(esc(p.get("impact_rank"))).append(": ")
                            .append(esc(expand(required(p, "impact_reason"), evidence, new ArrayList<>(), required(p, "id"), true))).append("</p>");
                        for (Object id : list(p.get("result_ids"))) h.append(" <a href=\"#").append(esc(id)).append("\">Discussed result</a>");
                    }
                    for (Object source : list(p.get("sources"))) h.append(sourceLink(source.toString()));
                    h.append("</details>");
                }
                for (Object item : list(p.get("supporting_evidence"))) {
                    var support = map(item); String id = support.get("source").toString();
                    h.append("<aside><h4>Related bioinformatics and contextual evidence</h4><p>")
                        .append(citations.render(expand(required(support, "relationship"), evidence, new ArrayList<>(), required(p, "id"), true), list(support.get("references")))).append("</p>")
                        .append(sourceLink(id)).append("<p>Source status: ").append(esc(evidence.get(id).get("state")))
                        .append(". This relationship is an editorial interpretation, not automatically verified support.</p></aside>");
                }
                if (section.equals("results")) for (Object id : list(p.get("sources"))) {
                    var source = evidence.get(id.toString());
                    if ("meta_analysis".equals(source.get("operation_type")) && StudyResultsDiscussion.completed(source) && previewed.add(id.toString())) {
                        h.append("<details><summary>Selected molecular findings for this analysis</summary><p>Source-selected previews, not independently validated discoveries; inspect effect scale and heterogeneity.</p>")
                            .append(table(list(map(source.get("key_findings")).get("top_features")).stream().limit(10).toList(), List.of("gene_symbol", "feature_id", "beta", "se", "p_value", "BH_FDR", "i_squared", "direction"))).append("</details>");
                    }
                }
                h.append("</section>");
                if ("limitations".equals(p.get("discussion_role"))) appendLimitations(h, evidence);
                if (section.equals("methods") && (p.containsKey("parameters") || p.containsKey("parameter_justifications")))
                    h.append("<details><summary>Parameters and recorded justifications</summary>").append(StudyPublication.render(Map.of("parameters", map(p.get("parameters")), "parameter_justifications", map(p.get("parameter_justifications"))))).append("</details>");
            }
            if (section.equals("methods")) {
                h.append("<p><a href=\"methods-template.yaml\">Methods template</a> · <a href=\"methods-review.yaml\">Methods review</a> · <a href=\"methods-inventory.yaml\">Execution and parameter inventory</a></p>");
                h.append("<details><summary>Declared cohort analyses, recorded parameters and unresolved inputs</summary>").append(StudyPublication.render(StudyMethods.displayInventory(StudyMethods.inventory(study, evidence)))).append("</details>");
                if (revision.containsKey("participant_project")) h.append("<p>Participant summaries were refreshed separately; see <a href=\"../metadata/participant-study.yaml\">the descriptive-summary contract</a>. Supplied association models were not refitted.</p>");
            }
            if (section.equals("introduction") && list(revision.get("sections")).stream().map(StudyManuscriptRevision::map).noneMatch(p -> "hypotheses_and_motivation".equals(p.get("introduction_role")))) {
                h.append("<p>The declared study hypotheses are: ");
                for (Object raw : list(study.get("hypotheses"))) if (!"retired".equals(map(raw).get("status"))) h.append(esc(map(raw).get("description"))).append(" ");
                h.append("</p>");
            }
            if (section.equals("discussion") && list(revision.get("sections")).stream().map(StudyManuscriptRevision::map).noneMatch(p -> p.containsKey("discussion_role"))) appendLimitations(h, evidence);
        }
        h.append("<h2>Reproducibility and use of artificial intelligence</h2><p>JLinAlg verified prior results and generated this revision locally. The revision request contains attributed editorial text; it may have been supplied by a scientist or an LLM advisor. The command made no remote model calls and read no patient-level phenotype rows. Checksums, comments, numeric bindings and before/after paragraph records are available in the audit bundle. This is integrity checking, not a guarantee that an editorial interpretation is correct.</p>");
        h.append("<p>Software: JLinAlg ").append(esc(JLinAlgCli.version())).append(". See <a href=\"../metadata/software-versions.yaml\">runtime and actor provenance</a>.</p>");
        h.append(citations.bibliography());
        return h.append(end()).toString();
    }

    private static void appendLimitations(StringBuilder h, Map<String, Map<String, Object>> evidence) {
        h.append("<details><summary>Outstanding evidence and required inputs</summary><ul>");
        for (var entry : evidence.entrySet()) if (!Set.of("complete", "not_applicable").contains(entry.getValue().get("state")))
            h.append("<li>").append(esc(entry.getKey())).append(": ").append(esc(entry.getValue().get("summary"))).append(sourceLink(entry.getKey())).append("</li>");
        h.append("</ul><p>Missing evidence is not negative evidence. When JLinAlg reports exploratory partial-data colocalization, its posterior probabilities are conditional on observed selected variants rather than conventional complete-region evidence. Missing variants, significance filtering, unequal coverage and winner's curse can materially bias the result; a null colocalization result is not evidence against a shared causal variant. Candidate counts are not independent instruments or causal estimates. Association and gene-level literature do not establish isoform-specific causality or drug efficacy.</p></details>");
    }

    private static String supplement(Map<String, Map<String, Object>> evidence) {
        StringBuilder h = new StringBuilder(start("Supplemental evidence and limitations"));
        h.append("<p><a href=\"manuscript.html\">Manuscript</a> · <a href=\"evidence.yaml\">Complete local evidence snapshot</a></p><p>All current operations are retained, including failed, waiting, null and contradictory evidence. Previews below are not exhaustive result tables. Original manifests record full result paths and hashes.</p><label>Find evidence <input id=\"search\" type=\"search\"></label>");
        h.append("<section class=\"evidence\" id=\"supplemental-table-1\"><h2>Supplemental Table 1. Detailed participant characteristics</h2>").append(StudyPublication.detailedTable1(participantTables(evidence))).append("</section>");
        for (var entry : evidence.entrySet()) {
            Map<String, Object> e = entry.getValue();
            h.append("<section class=\"evidence\" id=\"").append(esc(entry.getKey())).append("\"><h2>").append(esc(entry.getKey())).append("</h2><p>").append(esc(e.get("state"))).append(": ").append(esc(e.get("summary"))).append("</p>")
                .append("<p><a href=\"../metadata/sources/").append(esc(entry.getKey())).append(".yaml\">Verified source manifest</a>; SHA-256 ").append(esc(e.get("manifest_sha256"))).append("</p>");
            Object findings = e.get("key_findings");
            if (findings instanceof List<?> rows && !rows.isEmpty()) h.append(table(new ArrayList<>(rows), new ArrayList<>(map(rows.get(0)).keySet())));
            else h.append("<details><summary>Native findings and prerequisites</summary>").append(StudyPublication.render(findings)).append("</details>");
            h.append("</section>");
        }
        return h.append("<script>document.getElementById('search').addEventListener('input',function(){const q=this.value.toLowerCase();document.querySelectorAll('.evidence').forEach(e=>e.hidden=!e.textContent.toLowerCase().includes(q))});</script>").append(end()).toString();
    }

    private static String literature(Map<String, Map<String, Object>> evidence) {
        StringBuilder h = new StringBuilder(start("Literature candidates and relevance review"));
        h.append("<p><a href=\"manuscript.html\">Manuscript</a></p><p>Deduplicated by DOI, then source and article ID. One short previously extracted abstract fragment per article is shown. These are lexical retrieval candidates, not a systematic review or verified hypothesis support. Read full sources to assess species, design, tissue, direction, retractions and isoform relevance. Absence from this table does not mean absence of literature.</p>");
        Map<String, Map<String, Object>> articles = new LinkedHashMap<>();
        for (var source : evidence.entrySet()) if (source.getKey().endsWith(".literature_search")) {
            for (Object raw : list(map(source.getValue().get("key_findings")).get("evidence"))) {
                Map<String, Object> query = map(raw);
                for (Object value : list(query.get("articles"))) {
                    Map<String, Object> a = map(value);
                    String key = a.getOrDefault("doi", a.getOrDefault("source", "") + ":" + a.getOrDefault("id", "")).toString().toLowerCase(java.util.Locale.ROOT);
                    if (key.isBlank() || key.equals(":")) continue;
                    Map<String, Object> record = articles.computeIfAbsent(key, ignored -> new LinkedHashMap<>(a));
                    if (!record.containsKey("excerpt") && a.containsKey("excerpt")) record.put("excerpt", a.get("excerpt"));
                    record.computeIfAbsent("retrieval_context", ignored -> new ArrayList<String>());
                    @SuppressWarnings("unchecked") List<String> context = (List<String>) record.get("retrieval_context");
                    String contextText = source.getKey() + "; gene=" + query.get("gene_symbol") + "; " + query.get("search_strategy");
                    if (!context.contains(contextText)) context.add(contextText);
                    String title = String.valueOf(a.get("title")).toLowerCase(java.util.Locale.ROOT);
                    if ("FGR".equals(query.get("gene_symbol")) && (title.contains("fetal growth restriction") || title.contains("foetal growth restriction")))
                        record.put("relevance_warning", "Possible acronym collision: FGR denotes fetal growth restriction in this title; not established evidence about the FGR gene.");
                }
            }
        }
        h.append("<p>Unique retrieved article identities: ").append(articles.size()).append(".</p>");
        for (var a : articles.values()) {
            h.append("<article id=\"").append(referenceId(String.valueOf(a.getOrDefault("doi", a.get("id"))))).append("\"><h2>").append(esc(a.get("title"))).append("</h2><p>").append(esc(a.getOrDefault("authorString", ""))).append(" ").append(esc(a.getOrDefault("pubYear", ""))).append("</p>");
            String pmid = String.valueOf(a.getOrDefault("pmid", "")), doi = String.valueOf(a.getOrDefault("doi", ""));
            if (pmid.matches("[0-9]+")) h.append("<a href=\"https://pubmed.ncbi.nlm.nih.gov/").append(pmid).append("/\">PubMed source</a> ");
            if (!doi.isBlank()) h.append("<a href=\"https://doi.org/").append(esc(doi)).append("\">DOI ").append(esc(doi)).append("</a>");
            if (a.containsKey("excerpt")) {
                String excerpt = a.get("excerpt").toString();
                if (excerpt.strip().split("\\s+").length <= 12) h.append("<blockquote>").append(esc(excerpt)).append("</blockquote>");
            }
            h.append("<p>").append(esc(a.getOrDefault("relevance_warning", "Unreviewed candidate; gene-level or other-tissue findings do not validate the measured isoform."))).append("</p><p class=\"source\">").append(esc(a.get("retrieval_context"))).append("</p></article>");
        }
        return h.append(end()).toString();
    }

    static String expand(String template, Map<String, Map<String, Object>> evidence, List<Object> claims, String section) throws IOException {
        return expand(template, evidence, claims, section, false);
    }
    private static List<Object> participantTables(Map<String, Map<String, Object>> evidence) {
        List<Object> tables = new ArrayList<>();
        for (var item : evidence.entrySet()) if (item.getKey().endsWith(".participant_characteristics")) tables.addAll(list(map(item.getValue().get("key_findings")).get("tables")));
        return tables;
    }
    private static Set<String> participantIdentities(Map<String, Object> entry) {
        Set<String> identities = new LinkedHashSet<>();
        for (Object raw : list(map(entry.get("key_findings")).get("tables"))) {
            Map<String, Object> t = map(raw); identities.add(t.get("dataset") + "\n" + t.get("cohort"));
        }
        return identities;
    }
    private static String expand(String template, Map<String, Map<String, Object>> evidence, List<Object> claims, String section, boolean display) throws IOException {
        Matcher m = TOKEN.matcher(template); StringBuilder result = new StringBuilder();
        while (m.find()) {
            Map<String, Object> e = evidence.get(m.group(1));
            if (e == null) throw new IOException("Unknown bound source: " + m.group(1));
            Object value = e.get("key_findings");
            for (String part : m.group(2).substring(1).split("/", -1)) {
                String key = part.replace("~1", "/").replace("~0", "~");
                if (value instanceof Map<?, ?> map) value = map.get(key);
                else if (value instanceof List<?> list && key.matches("0|[1-9][0-9]*")) {
                    try { value = list.get(Integer.parseInt(key)); } catch (IndexOutOfBoundsException | NumberFormatException ex) { throw new IOException("Invalid evidence pointer: " + m.group(), ex); }
                } else throw new IOException("Invalid evidence pointer: " + m.group());
            }
            if (!(value instanceof String || value instanceof Number || value instanceof Boolean)) throw new IOException("Evidence binding requires an existing scalar: " + m.group());
            claims.add(Map.of("section_id", section, "binding", m.group(), "value", value, "source", m.group(1), "source_manifest", e.get("manifest")));
            String key = m.group(2).substring(m.group(2).lastIndexOf('/') + 1);
            m.appendReplacement(result, Matcher.quoteReplacement(display ? StudyReportNumbers.field(key, value) : value.toString()));
        }
        m.appendTail(result);
        if (result.toString().contains("{{") || result.toString().contains("}}")) throw new IOException("Malformed or unresolved evidence binding");
        return result.toString();
    }
    private static String table(List<Object> rows, List<String> columns) {
        if (rows.isEmpty()) return "<p>No preview rows available.</p>";
        StringBuilder h = new StringBuilder("<div class=\"table-wrap\"><table><thead><tr>");
        for (String column : columns) h.append("<th>").append(esc(label(column))).append("</th>");
        h.append("</tr></thead><tbody>");
        for (Object row : rows) { h.append("<tr>"); for (String column : columns) h.append("<td>").append(esc(StudyReportNumbers.field(column, map(row).getOrDefault(column, "—")))).append("</td>"); h.append("</tr>"); }
        return h.append("</tbody></table></div>").toString();
    }
    private static String sourceLink(String id) { return " <a class=\"source\" href=\"supplement.html#" + esc(id) + "\">[" + esc(id) + "]</a> "; }
    private static String referenceId(String doi) { return "ref-" + java.util.Base64.getUrlEncoder().withoutPadding().encodeToString(doi.toLowerCase(java.util.Locale.ROOT).getBytes(StandardCharsets.UTF_8)); }
    private static String start(String title) { return "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><meta http-equiv=\"Content-Security-Policy\" content=\"default-src 'none'; style-src 'unsafe-inline'; script-src 'unsafe-inline'; img-src 'none'; connect-src 'none'\"><title>" + esc(title) + "</title><style>body{font:17px/1.65 Georgia,serif;max-width:1080px;margin:2.5rem auto;padding:0 1.3rem;color:#182835;background:#fff}h1,h2,h3,h4,nav,.source,th{font-family:system-ui,sans-serif}h1{line-height:1.2;font-size:2.25rem}h2{border-top:2px solid #dbe4e9;padding-top:1rem;margin-top:2.5rem}h3{line-height:1.4}a{color:#125680;overflow-wrap:anywhere}.notice{padding:1rem;background:#fff7e5;border-left:4px solid #a66b00}.source{font-size:.82rem;color:#435767}table{border-collapse:collapse;font:14px/1.5 system-ui;display:block;overflow-x:auto;max-width:100%}td,th{padding:.5rem;border-bottom:1px solid #dbe4e9;text-align:left;vertical-align:top}th{background:#eff4f7}pre{font:13px/1.5 monospace;white-space:pre-wrap;overflow-wrap:anywhere}details{margin:1rem 0}dt{font-weight:bold}blockquote{border-left:3px solid #91afc1;padding-left:1rem}input{padding:.5rem}nav{margin:1.2rem 0}@media print{body{max-width:none;font-size:11pt}details> *{display:block}nav{display:none}}</style></head><body><h1>" + esc(title) + "</h1>"; }
    private static String end() { return "</body></html>"; }
    private static String label(String s) { return s.isEmpty() ? s : Character.toUpperCase(s.charAt(0)) + s.substring(1).replace('_', ' '); }
    private static String esc(Object value) { return String.valueOf(value == null ? "Not reported" : value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("\"", "&quot;").replace("'", "&#39;"); }
    private static void write(Path path, String value) throws IOException { Files.writeString(path, value, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW); }
    private static void writeYaml(Path path, Object value) throws IOException { write(path, new Yaml().dump(value)); }
    private static Map<String, Object> load(Path path) throws IOException { LoaderOptions o = new LoaderOptions(); o.setAllowDuplicateKeys(false); o.setCodePointLimit(32 * 1024 * 1024); return map(new Yaml(new SafeConstructor(o)).load(Files.readString(path))); }
    private static Map<String, Object> map(Object o) { Map<String, Object> m = new LinkedHashMap<>(); if (o instanceof Map<?, ?> raw) raw.forEach((k, v) -> m.put(k.toString(), v)); return m; }
    private static List<Object> list(Object o) { return o instanceof List<?> l ? new ArrayList<>(l) : List.of(); }
    private static String required(Map<String, Object> m, String key) throws IOException { if (!(m.get(key) instanceof String s) || s.isBlank()) throw new IOException("Required nonempty text: " + key); return (String) m.get(key); }
    private static void identifier(String id) throws IOException { if (!id.matches("[A-Za-z][A-Za-z0-9._-]{0,127}")) throw new IOException("Invalid revision or paragraph ID"); }
    private static void keys(Map<String, Object> m, Set<String> allowed, String name) throws IOException { for (String key : m.keySet()) if (!allowed.contains(key)) throw new IOException("Unknown " + name + " field: " + key); }
    private static Map<String, Object> input(Path path, String hash, String role) { return Map.of("path", path.toString(), "sha256", hash, "role", role); }
    private static String relative(Path from, Path to) { return from.relativize(to).toString().replace('\\', '/'); }
    private static Path contained(Path base, String name) throws IOException { Path root = base.toRealPath(), p = root.resolve(name).normalize(); if (!p.startsWith(root) || !p.toRealPath().startsWith(root) || !Files.isRegularFile(p)) throw new IOException("Evidence path escapes bundle or is absent"); return p; }
}
