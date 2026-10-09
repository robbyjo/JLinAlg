/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.net.URI;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.time.Instant;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.jlinalg.study.StudyOperationExecutor;
import org.jlinalg.study.StudyProject;
import org.yaml.snakeyaml.LoaderOptions;
import org.yaml.snakeyaml.Yaml;
import org.yaml.snakeyaml.constructor.SafeConstructor;

/** Offline, browser-independent exchange. Browser events are human attestations, not telemetry. */
final class StudyBrowserHandoff {
    private StudyBrowserHandoff() { }
    static void verify(Path directory) throws IOException { verifyPacket(directory, load(directory.resolve("packet.yaml"))); }
    static final String GAP = "JLinAlg did not observe the browser submission or authenticate the provider, model, timestamps, answer or citations. These are human-attested records. Hashes detect changes, not scientific truth or provider authenticity.";

    static Map<String, Object> prepare(StudyOperationExecutor.Context context, String purpose,
            String prompt, Map<String, Object> policy) throws IOException {
        if (prompt.length() > 150000) throw new IOException("Browser prompt exceeds 150000 characters");
        Path dir = context.outputDirectory().resolve("handoff");
        Files.createDirectory(dir);
        Map<String, Object> packet = new LinkedHashMap<>();
        packet.put("schema_version", 1); packet.put("purpose", purpose);
        packet.put("hypothesis", context.operation().getOrDefault("hypothesis", "editorial"));
        packet.put("project_id", map(context.study().get("project")).getOrDefault("id", "not_declared"));
        packet.put("resolved_study_hash", map(context.study().get("_resolution")).getOrDefault("resolved_study_hash", "not_available"));
        Map<String, Object> sources = new LinkedHashMap<>();
        for (var entry : context.dependencyDirectories().entrySet()) {
            Path manifest = entry.getValue().resolve("manifest.yaml");
            if (Files.isRegularFile(manifest)) {
                Map<String, Object> m = load(manifest);
                sources.put(entry.getKey(), Map.of("state", m.getOrDefault("state", "unknown"),
                    "operation_type", m.getOrDefault("operation_type", "unknown"), "output_sha256", m.getOrDefault("output_sha256", Map.of())));
            }
        }
        packet.put("source_identity_sha256", hash(StudyAssistance.json(sources)));
        packet.put("prompt_sha256", hash(prompt)); packet.put("data_policy", policy);
        Map<String, Object> llm = map(policy.get("llm"));
        boolean allowed = !StudyAssistance.denied(context.study())
            && !llm.containsKey("minimum_cell_count") && !llm.containsKey("eligible_remote_artifacts");
        packet.put("disclosure_permitted_after_human_review", allowed);
        packet.put("disclosure_note", allowed
            ? "Inspect the exact prompt before copying. Review patient privacy AND unpublished-research confidentiality against the destination's terms. Only explicit human approval permits manual submission."
            : "DO NOT SUBMIT: project policy denies disclosure or requires unsupported artifact/cell-count checks. An approved policy amendment is needed, not a browser workaround.");
        packet.put("request_id", hash(StudyAssistance.json(packet)));
        packet.put("created_at", Instant.now().toString());
        packet.put("jlinalg_network_requests", 0); packet.put("audit_gap", GAP);
        Files.writeString(dir.resolve("prompt.txt"), prompt, StandardOpenOption.CREATE_NEW);
        Map<String, Object> template = new LinkedHashMap<>();
        template.put("schema_version", 1); template.put("request_id", packet.get("request_id"));
        template.put("prompt_sha256", packet.get("prompt_sha256"));
        template.put("provider", ""); template.put("model", "not_disclosed"); template.put("url", "");
        template.put("submitted_at", ""); template.put("received_at", "");
        template.put("exact_prompt_submitted", false); template.put("disclosure_approved", false);
        template.put("no_patient_data", false); template.put("response", ""); template.put("references", List.of());
        template.put("decision", "pending_review"); template.put("rationale", "");
        write(dir.resolve("response-template.yaml"), template);
        String html = page("Audited browser handoff", "<p>" + esc(packet.get("disclosure_note"))
            + "</p><p>No browser opened and no data transmitted by JLinAlg. Do not upload the project, phenotype files, or this entire folder.</p>"
            + "<ol><li>Read <a href=\"prompt.txt\">the exact prompt</a> and review disclosure policy.</li>"
            + "<li>If permitted, manually copy only that prompt into your chosen AI website in a normal browser.</li>"
            + "<li>Copy <a href=\"response-template.yaml\">the response template</a> OUTSIDE this frozen packet. Fill in the complete answer, references, provider/model label, UTC-offset timestamps and attestations. Use not_disclosed when the model is unknown.</li>"
            + "<li>Import it with: <code>java -jar jlinalg.jar study handoff-import RESPONSE.yaml --handoff PACKET_DIRECTORY --out NEW_RECEIPT_DIRECTORY --actor-id SCIENTIST</code>.</li>"
            + "<li>Inspect the receipt. Importing does not apply manuscript edits or authorize analyses. Accept selected changes through a new manuscript revision or study amendment.</li></ol>"
            + "<p>Request ID: <code>" + packet.get("request_id") + "</code><br>Prompt SHA-256: <code>" + packet.get("prompt_sha256")
            + "</code></p><p>" + esc(GAP) + "</p><p>This local page needs no JavaScript, network, account or browser extension. Linux Links can read it; use a modern browser on an authorized workstation for AI websites. No private endpoint or cookie replay is used.</p>"
            + "<h2>Exact prompt</h2><pre>" + esc(prompt) + "</pre>");
        Files.writeString(dir.resolve("index.html"), html, StandardOpenOption.CREATE_NEW);
        packet.put("files_sha256", Map.of("prompt.txt", FollowupSupport.hash(dir.resolve("prompt.txt")),
            "index.html", FollowupSupport.hash(dir.resolve("index.html")), "response-template.yaml", FollowupSupport.hash(dir.resolve("response-template.yaml"))));
        write(dir.resolve("packet.yaml"), packet);
        return packet;
    }

    static Path importReply(Path handoff, Path replyFile, Path output, StudyProject.Actor actor) throws IOException {
        if (!"human".equals(actor.type()) || "unspecified".equals(actor.id()) || actor.id().isBlank())
            throw new IOException("Browser import requires a named human --actor-id; an LLM cannot attest browser disclosure");
        Path dir = handoff.toRealPath(), out = output.toAbsolutePath().normalize();
        if (replyFile.toRealPath().startsWith(dir)) throw new IOException("Copy and fill the response template outside the frozen handoff packet");
        if (Files.exists(out) || out.startsWith(dir) || out.equals(replyFile.toAbsolutePath().getParent()))
            throw new IOException("Choose a new receipt directory outside the frozen packet");
        Map<String, Object> packet = load(dir.resolve("packet.yaml"));
        verifyPacket(dir, packet);
        String packetHash = FollowupSupport.hash(dir.resolve("packet.yaml"));
        String originalReply = readBounded(replyFile), replyHash = hash(originalReply);
        Map<String, Object> reply = parse(originalReply);
        validateReply(packet, reply);
        String prompt = readBounded(dir.resolve("prompt.txt"));
        Files.createDirectories(out.getParent()); Files.createDirectory(out);
        Files.writeString(out.resolve("response.yaml"), originalReply, StandardOpenOption.CREATE_NEW);
        Files.copy(dir.resolve("packet.yaml"), out.resolve("packet.yaml"));
        Files.writeString(out.resolve("prompt.txt"), prompt, StandardOpenOption.CREATE_NEW);
        Map<String, Object> receipt = new LinkedHashMap<>();
        receipt.put("schema_version", 1); receipt.put("bundle_type", "browser_handoff_receipt");
        receipt.put("request_id", packet.get("request_id")); receipt.put("prompt_sha256", packet.get("prompt_sha256"));
        receipt.put("packet_sha256", packetHash); receipt.put("response_sha256", replyHash);
        receipt.put("imported_at", Instant.now().toString()); receipt.put("actor", actor.asMap());
        receipt.put("purpose", packet.get("purpose")); receipt.put("hypothesis", packet.get("hypothesis"));
        receipt.put("audit_gap", GAP); receipt.put("disposition", "recorded_only_not_applied");
        receipt.put("jlinalg_network_requests", 0);
        String response = required(reply, "response");
        Files.writeString(out.resolve("index.html"), page("Browser response audit",
            "<p>Recorded only. No analyses executed or manuscript edits applied.</p><p>" + esc(GAP)
            + "</p><p>Provider: " + esc(reply.get("provider")) + "; model label: " + esc(reply.get("model"))
            + "; decision: " + esc(reply.get("decision")) + "</p><p>Rationale: " + esc(reply.get("rationale"))
            + "</p><p><a href=\"receipt.yaml\">Receipt and hashes</a> | <a href=\"prompt.txt\">Exact prompt</a> | <a href=\"response.yaml\">Attested response and references</a></p><pre>" + esc(response) + "</pre>"), StandardOpenOption.CREATE_NEW);
        Files.writeString(out.resolve("SESSION_HISTORY.md"), "# Browser handoff notebook\n\n" + receipt.get("imported_at")
            + " — JLinAlg " + JLinAlgCli.version() + " imported a browser transcript attested by " + actor.asMap()
            + ". Provider/model and external event timestamps are in response.yaml; they were not observed by JLinAlg. No network calls, automatic edits or analysis execution.\n\n" + GAP + "\n", StandardOpenOption.CREATE_NEW);
        receipt.put("software", Map.of("jlinalg", JLinAlgCli.version(), "java", System.getProperty("java.runtime.version")));
        try {
            Path runtime = Path.of(StudyBrowserHandoff.class.getProtectionDomain().getCodeSource().getLocation().toURI());
            receipt.put("runtime_sha256", Files.isRegularFile(runtime) ? FollowupSupport.hash(runtime) : "unknown: class-directory execution");
        } catch (java.net.URISyntaxException e) { throw new IOException(e); }
        receipt.put("files_sha256", Map.of("packet.yaml", packetHash, "prompt.txt", hash(prompt), "response.yaml", replyHash,
            "index.html", FollowupSupport.hash(out.resolve("index.html")), "SESSION_HISTORY.md", FollowupSupport.hash(out.resolve("SESSION_HISTORY.md"))));
        verifyPacket(dir, packet);
        if (!packetHash.equals(FollowupSupport.hash(dir.resolve("packet.yaml"))) || !replyHash.equals(FollowupSupport.hash(replyFile)))
            throw new IOException("Handoff source changed during import; receipt remains incomplete");
        write(out.resolve("receipt.yaml"), receipt);
        verifyReceipt(out.resolve("receipt.yaml"));
        return out;
    }

    static Map<String, Object> verifyReceipt(Path file) throws IOException {
        Map<String, Object> receipt = load(file);
        if (!"browser_handoff_receipt".equals(receipt.get("bundle_type"))) throw new IOException("Not a browser receipt");
        if (!"receipt.yaml".equals(file.getFileName().toString())) throw new IOException("Browser receipt must be named receipt.yaml");
        Path dir = file.toAbsolutePath().getParent();
        verifyFiles(dir, map(receipt.get("files_sha256")), Set.of("packet.yaml", "prompt.txt", "response.yaml", "index.html", "SESSION_HISTORY.md"));
        Map<String, Object> packet = load(dir.resolve("packet.yaml")), reply = load(dir.resolve("response.yaml"));
        verifyIdentity(packet);
        Map<String, Object> actor = map(receipt.get("actor"));
        if (!"human".equals(actor.get("type")) || required(actor, "id").equals("unspecified")
                || !packet.get("purpose").equals(receipt.get("purpose")) || !packet.get("hypothesis").equals(receipt.get("hypothesis"))) throw new IOException("Invalid browser receipt attribution");
        if (!receipt.get("request_id").equals(packet.get("request_id")) || !receipt.get("prompt_sha256").equals(hash(readBounded(dir.resolve("prompt.txt"))))
                || !receipt.get("packet_sha256").equals(FollowupSupport.hash(dir.resolve("packet.yaml")))
                || !receipt.get("response_sha256").equals(FollowupSupport.hash(dir.resolve("response.yaml")))) throw new IOException("Browser receipt identity mismatch");
        validateReply(packet, reply);
        return reply;
    }

    static Map<String, Object> matchingReply(Path project, Map<String, Object> config, Map<String, Object> packet, Path output) throws IOException {
        Path match = null; Map<String, Object> answer = null;
        for (Object value : list(config.get("response_files"))) {
            Path file = project.resolve(value.toString()).normalize();
            if (!file.getFileName().toString().equals("receipt.yaml")) throw new IOException("response_files entries must end in receipt.yaml");
            if (!Files.exists(file)) continue;
            Map<String, Object> receipt = load(file);
            if (!packet.get("request_id").equals(receipt.get("request_id"))) continue;
            if (match != null) throw new IOException("Multiple browser receipts match the same request; choose one explicitly");
            answer = verifyReceipt(file); match = file;
        }
        if (match != null) {
            Path target = output.resolve("browser-response"); Files.createDirectory(target);
            for (String name : List.of("receipt.yaml", "packet.yaml", "prompt.txt", "response.yaml", "index.html", "SESSION_HISTORY.md"))
                Files.copy(match.getParent().resolve(name), target.resolve(name));
            verifyReceipt(target.resolve("receipt.yaml"));
            answer.put("_import_actor", load(target.resolve("receipt.yaml")).get("actor"));
        }
        return answer;
    }

    private static void validateReply(Map<String, Object> packet, Map<String, Object> reply) throws IOException {
        Set<String> fields = Set.of("schema_version", "request_id", "prompt_sha256", "provider", "model", "url", "submitted_at", "received_at",
            "exact_prompt_submitted", "disclosure_approved", "no_patient_data", "response", "references", "decision", "rationale");
        if (!fields.equals(reply.keySet()) || !Integer.valueOf(1).equals(reply.get("schema_version"))) throw new IOException("Invalid browser response fields/schema");
        if (!Boolean.TRUE.equals(packet.get("disclosure_permitted_after_human_review"))) throw new IOException("Policy does not permit this handoff disclosure");
        if (!packet.get("request_id").equals(reply.get("request_id")) || !packet.get("prompt_sha256").equals(reply.get("prompt_sha256"))) throw new IOException("Browser response does not match the exact request and prompt");
        for (String field : List.of("exact_prompt_submitted", "disclosure_approved", "no_patient_data"))
            if (!Boolean.TRUE.equals(reply.get(field))) throw new IOException("Human attestation required: " + field);
        for (String field : List.of("provider", "model", "rationale", "response")) required(reply, field);
        if (required(reply, "provider").length() > 200 || required(reply, "model").length() > 200 || required(reply, "response").length() > 250000) throw new IOException("Browser response exceeds bounds");
        safeUrl(required(reply, "url"));
        try {
            Instant sent = Instant.parse(required(reply, "submitted_at")), received = Instant.parse(required(reply, "received_at"));
            if (sent.isAfter(received) || received.isAfter(Instant.now().plusSeconds(300)) || sent.isBefore(Instant.parse(required(packet, "created_at")))) throw new IOException("Browser timestamps must follow preparation, be ordered and not in the future");
        } catch (java.time.DateTimeException invalid) { throw new IOException("Use ISO-8601 timestamps with UTC offset", invalid); }
        if (!Set.of("pending_review", "accepted_for_revision", "rejected").contains(reply.get("decision"))) throw new IOException("Unknown browser review decision");
        if (!(reply.get("references") instanceof List<?> refs) || refs.size() > 100) throw new IOException("References must be a list of at most 100 entries");
        for (Object raw : list(reply.get("references"))) {
            Map<String, Object> ref = map(raw);
            if (!ref.keySet().equals(Set.of("title", "url"))) throw new IOException("References require title and url only");
            required(ref, "title"); safeUrl(required(ref, "url"));
        }
    }

    private static void verifyPacket(Path dir, Map<String, Object> packet) throws IOException {
        verifyIdentity(packet);
        verifyFiles(dir, map(packet.get("files_sha256")), Set.of("prompt.txt", "index.html", "response-template.yaml"));
        if (!packet.get("prompt_sha256").equals(FollowupSupport.hash(dir.resolve("prompt.txt")))) throw new IOException("Prompt hash mismatch");
    }
    private static void verifyIdentity(Map<String, Object> packet) throws IOException {
        Map<String, Object> identity = new LinkedHashMap<>(packet);
        Map<String, Object> policy = map(packet.get("data_policy")), llm = map(policy.get("llm"));
        boolean permitted = !"deny_all".equals(policy.get("remote_transmission")) && !"deny_all".equals(llm.get("remote_transmission"))
            && !llm.containsKey("minimum_cell_count") && !llm.containsKey("eligible_remote_artifacts");
        if (!Boolean.valueOf(permitted).equals(packet.get("disclosure_permitted_after_human_review"))) throw new IOException("Handoff disclosure policy mismatch");
        for (String field : List.of("request_id", "created_at", "jlinalg_network_requests", "audit_gap", "files_sha256")) identity.remove(field);
        if (!hash(StudyAssistance.json(identity)).equals(packet.get("request_id"))) throw new IOException("Handoff request identity mismatch");
    }
    private static void verifyFiles(Path dir, Map<String, Object> hashes, Set<String> required) throws IOException {
        if (!hashes.keySet().equals(required)) throw new IOException("Incomplete or unexpected browser inventory");
        for (String name : required) {
            Path p = dir.resolve(name);
            if (Files.isSymbolicLink(p) || !hashes.get(name).equals(FollowupSupport.hash(p))) throw new IOException("Browser artifact hash mismatch: " + name);
        }
    }
    private static void safeUrl(String url) throws IOException {
        try {
            URI u = URI.create(url);
            if (!"https".equals(u.getScheme()) || u.getHost() == null || u.getUserInfo() != null || u.getFragment() != null || u.getQuery() != null)
                throw new IllegalArgumentException();
        } catch (IllegalArgumentException bad) { throw new IOException("Use an HTTPS URL without credentials, query parameters or fragments; do not record session tokens"); }
    }
    static String hash(String text) { try { return HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256").digest(text.getBytes(StandardCharsets.UTF_8))); } catch (NoSuchAlgorithmException e) { throw new IllegalStateException(e); } }
    private static String readBounded(Path file) throws IOException { if (Files.size(file) > 2000000) throw new IOException("Browser artifact exceeds 2 MB"); return Files.readString(file); }
    private static Map<String, Object> load(Path file) throws IOException { return parse(readBounded(file)); }
    private static Map<String, Object> parse(String text) throws IOException { LoaderOptions o = new LoaderOptions(); o.setAllowDuplicateKeys(false); o.setMaxAliasesForCollections(0); try { Object value = new Yaml(new SafeConstructor(o)).load(text); if (!(value instanceof Map<?, ?>)) throw new IOException("Expected a YAML mapping"); return map(value); } catch (org.yaml.snakeyaml.error.YAMLException e) { throw new IOException("Invalid browser YAML", e); } }
    private static String required(Map<String, Object> m, String k) throws IOException { if (!(m.get(k) instanceof String s) || s.isBlank()) throw new IOException("Required browser field: " + k); return s; }
    private static void write(Path p, Object value) throws IOException { Files.writeString(p, new Yaml().dump(value), StandardOpenOption.CREATE_NEW); }
    private static Map<String, Object> map(Object o) { Map<String, Object> r = new LinkedHashMap<>(); if (o instanceof Map<?, ?> m) m.forEach((k,v) -> r.put(k.toString(),v)); return r; }
    private static List<?> list(Object o) { return o instanceof List<?> l ? l : List.of(); }
    private static String esc(Object o) { return String.valueOf(o).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("\"", "&quot;"); }
    private static String page(String title, String body) { return "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"referrer\" content=\"no-referrer\"><meta http-equiv=\"Content-Security-Policy\" content=\"default-src 'none'; style-src 'unsafe-inline'; form-action 'none'; base-uri 'none'\"><title>" + esc(title) + "</title><style>body{max-width:70rem;margin:2rem auto;padding:1rem;font-family:system-ui}pre{white-space:pre-wrap;overflow-wrap:anywhere}</style></head><body><h1>" + esc(title) + "</h1>" + body + "</body></html>"; }
}
