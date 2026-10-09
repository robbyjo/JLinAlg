/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.net.URI;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Instant;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.jlinalg.study.StudyOperationExecutor.Context;
import org.yaml.snakeyaml.Yaml;

/** Bounded, audited public evidence requests. No input file is ever uploaded. */
final class StudyPublicEvidence {
    static final Set<String> HOSTS = Set.of("www.ebi.ac.uk", "ftp.ebi.ac.uk",
        "rest.uniprot.org", "alphafold.ebi.ac.uk", "api.genome.ucsc.edu", "gtexportal.org", "api.crossref.org",
        "api.ncbi.nlm.nih.gov", "pan-ukb-us-east-1.s3.amazonaws.com",
        "ftp.1000genomes.ebi.ac.uk", "hgdownload.soe.ucsc.edu");
    record Reply(int status, String body) { }
    @FunctionalInterface interface Transport {
        Reply get(URI uri) throws IOException, InterruptedException;
    }
    record Response(String status, Object data, String audit, String reason) {
        boolean ok() { return status.equals("retrieved"); }
    }
    private final Transport transport;
    private final Map<String, Response> cache = new LinkedHashMap<>();
    private final List<Object> requests = new ArrayList<>();
    private static long lastUcscRequest;

    StudyPublicEvidence() { this(StudyPublicEvidence::send); }
    StudyPublicEvidence(Transport transport) { this.transport = transport; }
    List<Object> requests() { return List.copyOf(requests); }

    Response get(Context context, String provider, String url, String citation, String rationale) throws IOException {
        if (rationale == null || rationale.isBlank() || citation == null || citation.isBlank())
            throw new IllegalArgumentException("Every evidence request requires a scientific rationale and source citation");
        URI uri = URI.create(url);
        validate(uri);
        String cacheKey = context.outputDirectory().toAbsolutePath().normalize() + "\n" + url;
        Path folder = context.outputDirectory().resolve("queries").resolve(String.format("%04d", requests.size() + 1));
        Files.createDirectories(folder);
        String relative = context.outputDirectory().relativize(folder).toString().replace('\\', '/');
        Map<String, Object> record = new LinkedHashMap<>();
        record.put("provider", provider);
        record.put("url", url);
        record.put("citation", citation);
        record.put("scientific_rationale", rationale);
        record.put("hypothesis_id", context.operation().get("hypothesis"));
        record.put("operation", context.operation().get("type"));
        record.put("started_at", Instant.now().toString());
        record.put("request_attempted", false);
        record.put("payload", "URL query containing public feature identifiers, phenotype labels or reference intervals only; no files or participant records");
        String status = blocked(context.study());
        if (status.isEmpty() && requests.size() >= 200) status = "query_budget_exhausted";
        String reason = status.isEmpty() ? "" : "Request not sent: " + status;
        Object data = Map.of();
        // Each use retains its own rationale. A cached response must not bypass a
        // subsequently tightened disclosure policy or refer into another operation.
        if (status.isEmpty() && cache.containsKey(cacheKey)) {
            Response cached = cache.get(cacheKey);
            status = cached.status(); data = cached.data(); reason = cached.reason();
            record.put("cached_response_audit", cached.audit());
        }
        if (status.isEmpty()) {
            record.put("request_attempted", true);
            write(folder.resolve("request.yaml"), record);
            try {
                Reply reply = transport.get(uri);
                record.put("http_status", reply.status());
                // Provider responses are untrusted data, never instructions or executable code.
                if (reply.body().length() > 8 * 1024 * 1024) throw new IOException("Response exceeded 8 MiB limit");
                Files.writeString(folder.resolve("response.txt"), reply.body());
                record.put("response_sha256", FollowupSupport.hash(folder.resolve("response.txt")));
                if (reply.status() >= 200 && reply.status() < 300) {
                    data = SimpleJson.parse(reply.body());
                    status = "retrieved";
                } else {
                    status = switch (reply.status()) {
                        case 401, 403 -> "access_denied";
                        case 404 -> "not_found";
                        case 410 -> "endpoint_retired";
                        case 429 -> "rate_limited";
                        default -> "http_error";
                    };
                    reason = "HTTP " + reply.status() + "; not evidence of no scientific matches";
                }
            } catch (InterruptedException exception) {
                Thread.currentThread().interrupt();
                throw new IOException("Evidence acquisition interrupted", exception);
            } catch (IOException exception) {
                status = "unavailable";
                reason = String.valueOf(exception.getMessage());
            }
        }
        record.put("finished_at", Instant.now().toString());
        record.put("status", status);
        record.put("reason", reason);
        record.put("audit", relative);
        write(folder.resolve("request.yaml"), record);
        requests.add(record);
        Response response = new Response(status, data, relative, reason);
        cache.putIfAbsent(cacheKey, response);
        return response;
    }

    static String blocked(Map<String, Object> study) {
        Map<String, Object> policy = map(map(study.get("project")).get("data_policy"));
        Map<String, Object> llm = map(policy.get("llm"));
        if ("deny_all".equals(policy.get("remote_transmission")) || "deny_all".equals(llm.get("remote_transmission")))
            return "policy_denied";
        if (Boolean.TRUE.equals(llm.get("require_human_approval")) || llm.containsKey("minimum_cell_count")
                || llm.containsKey("eligible_remote_artifacts")) return "disclosure_review_required";
        if ("offline".equals(map(study.get("resources")).get("network"))) return "offline_policy";
        if (!"automatic".equals(map(study.get("resources")).getOrDefault("downloads", "automatic")))
            return "automatic_downloads_disabled";
        return "";
    }

    /** Generic public dataset installation never transmits a study payload. */
    static String resourceInstallBlocked(Map<String, Object> study) {
        Map<String, Object> resources = map(study.get("resources"));
        if ("offline".equals(resources.get("network")))
            return "offline_policy";
        if (!"automatic".equals(resources.getOrDefault("downloads",
                "automatic")))
            return "automatic_downloads_disabled";
        return "";
    }

    static void validate(URI uri) {
        if (!"https".equals(uri.getScheme()) || !HOSTS.contains(uri.getHost())
                || uri.getUserInfo() != null || uri.getFragment() != null
                || (uri.getPort() != -1 && uri.getPort() != 443))
            throw new IllegalArgumentException("Evidence URL must use a supported public HTTPS host without credentials or redirects");
    }

    private static Reply send(URI uri) throws IOException, InterruptedException {
        validate(uri);
        if (uri.getHost().equals("api.genome.ucsc.edu")) synchronized (StudyPublicEvidence.class) {
            long remaining = 1000000000L - (System.nanoTime() - lastUcscRequest);
            if (remaining > 0) java.util.concurrent.TimeUnit.NANOSECONDS.sleep(remaining);
            lastUcscRequest = System.nanoTime();
        }
        var response = PublicResourceDownloader.responseStrict(uri, 8 * 1024 * 1024);
        return new Reply(response.status(), response.body());
    }

    static String encode(String value) { return URLEncoder.encode(value, StandardCharsets.UTF_8); }
    static String text(Object value) { return value == null ? "" : String.valueOf(value); }
    static Map<String, Object> map(Object value) {
        Map<String, Object> result = new LinkedHashMap<>();
        if (value instanceof Map<?, ?> source) source.forEach((key, item) -> result.put(String.valueOf(key), item));
        return result;
    }
    static List<?> list(Object value) { return value instanceof List<?> items ? items : List.of(); }
    static void write(Path path, Object value) throws IOException {
        // Cached provider responses may be referenced by several isoforms. Emit
        // independent containers so the strict study reader never sees aliases.
        Files.writeString(path, new Yaml().dump(detach(value)));
    }

    private static Object detach(Object value) {
        if (value instanceof Map<?, ?> source) {
            Map<String, Object> result = new LinkedHashMap<>();
            source.forEach((key, item) -> result.put(String.valueOf(key), detach(item)));
            return result;
        }
        if (value instanceof List<?> source) {
            List<Object> result = new ArrayList<>(source.size());
            source.forEach(item -> result.add(detach(item)));
            return result;
        }
        return value;
    }
}
