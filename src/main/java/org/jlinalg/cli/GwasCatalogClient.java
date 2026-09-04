/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.net.URLEncoder;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.charset.StandardCharsets;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/** Public NHGRI-EBI GWAS Catalog client used by the MR instrument CLI. */
final class GwasCatalogClient implements InstrumentCatalog {
    private static final URI API = URI.create(
        "https://www.ebi.ac.uk/gwas/rest/api/v2/");
    private static final Duration REQUEST_TIMEOUT = Duration.ofSeconds(60);
    private static final Pattern DOWNLOAD = Pattern.compile(
        "(?i)href=\"([^\"/?#]+\\.(?:tsv|txt)(?:\\.gz)?)\"");
    private final HttpClient http;

    GwasCatalogClient() {
        http = HttpClient.newBuilder()
            .connectTimeout(Duration.ofSeconds(20))
            .version(HttpClient.Version.HTTP_1_1)
            .followRedirects(HttpClient.Redirect.NORMAL)
            .build();
    }

    @Override public List<InstrumentStudy> search(String trait, int limit)
            throws IOException, InterruptedException {
        String query = "studies?efo_trait=" + encode(trait)
            + "&full_pvalue_set=true&size=" + limit;
        Map<String, Object> root = object(json(API.resolve(query)));
        Object embeddedValue = root.get("_embedded");
        if (!(embeddedValue instanceof Map<?, ?> embedded)) return List.of();
        Object studiesValue = embedded.get("studies");
        if (!(studiesValue instanceof List<?> studies)) return List.of();
        List<InstrumentStudy> result = new ArrayList<>();
        for (Object value : studies) {
            if (value instanceof Map<?, ?> item)
                result.add(studyFrom(item));
        }
        return List.copyOf(result);
    }

    @Override public InstrumentStudy study(String accession)
            throws IOException, InterruptedException {
        return studyFrom(object(json(API.resolve(
            "studies/" + encode(accession)))));
    }

    @Override public RemoteSummary openSummaryStatistics(InstrumentStudy study)
            throws IOException, InterruptedException {
        URI location = secureEbiUri(study.summaryStatistics());
        String path = location.getPath().toLowerCase(java.util.Locale.ROOT);
        URI source = path.endsWith(".gz") || path.endsWith(".tsv")
            || path.endsWith(".txt") ? location : discover(location, study.accession());
        HttpResponse<InputStream> response = send(source,
            HttpResponse.BodyHandlers.ofInputStream());
        if (response.statusCode() / 100 != 2) {
            response.body().close();
            throw new IOException("GWAS summary-statistics server returned HTTP "
                + response.statusCode() + " for " + source);
        }
        return new RemoteSummary(source, response.body());
    }

    private URI discover(URI directory, String accession)
            throws IOException, InterruptedException {
        URI base = URI.create(directory.toString().replaceAll("/*$", "/"));
        String listing = text(base);
        List<String> files = new ArrayList<>();
        Matcher matcher = DOWNLOAD.matcher(listing);
        while (matcher.find()) files.add(matcher.group(1));
        if (files.isEmpty()) throw new IOException(
            "no TSV summary-statistics file was listed at " + base);
        files.sort(Comparator.comparingInt(file -> rank(file, accession)));
        return secureEbiUri(base.resolve(files.get(0)));
    }

    private static int rank(String file, String accession) {
        String lower = file.toLowerCase(java.util.Locale.ROOT);
        String id = accession.toLowerCase(java.util.Locale.ROOT);
        if (lower.equals(id + ".h.tsv.gz")) return 0;
        if (lower.equals(id + ".tsv.gz")) return 1;
        if (lower.endsWith(".h.tsv.gz")) return 2;
        if (lower.endsWith(".tsv.gz")) return 3;
        if (lower.endsWith(".tsv")) return 4;
        return 5;
    }

    private InstrumentStudy studyFrom(Map<?, ?> item) throws IOException {
        String accession = requiredString(item, "accession_id");
        String trait = requiredString(item, "disease_trait");
        String summary = requiredString(item, "full_summary_stats");
        return new InstrumentStudy(accession, trait,
            nestedStrings(item.get("efo_traits"), "efo_trait"),
            strings(item.get("discovery_ancestry")),
            string(item.get("initial_sample_size")),
            number(item.get("snp_count")), secureEbiUri(URI.create(summary)),
            optionalUri(item.get("terms_of_license")));
    }

    private String json(URI uri) throws IOException, InterruptedException {
        return text(uri);
    }

    private String text(URI uri) throws IOException, InterruptedException {
        HttpResponse<String> response = send(uri,
            HttpResponse.BodyHandlers.ofString(StandardCharsets.UTF_8));
        if (response.statusCode() / 100 != 2)
            throw new IOException("GWAS Catalog returned HTTP "
                + response.statusCode() + " for " + uri);
        return response.body();
    }

    private <T> HttpResponse<T> send(URI uri,
            HttpResponse.BodyHandler<T> handler)
            throws IOException, InterruptedException {
        HttpRequest request = HttpRequest.newBuilder(uri)
            .timeout(REQUEST_TIMEOUT)
            .header("Accept", "application/json, text/html;q=0.8")
            .header("User-Agent", "JLinAlg/0.1 MR-instrument-client")
            .GET().build();
        return http.send(request, handler);
    }

    private static URI secureEbiUri(URI input) throws IOException {
        String host = input.getHost();
        if (host == null || !(host.equalsIgnoreCase("ftp.ebi.ac.uk")
                || host.equalsIgnoreCase("www.ebi.ac.uk")))
            throw new IOException("GWAS Catalog returned an untrusted download host: "
                + input);
        try {
            return new URI("https", input.getUserInfo(), host, input.getPort(),
                input.getPath(), input.getQuery(), input.getFragment());
        } catch (java.net.URISyntaxException exception) {
            throw new IOException("invalid GWAS Catalog download URI", exception);
        }
    }

    private static Map<String, Object> object(String json) throws IOException {
        Object parsed = SimpleJson.parse(json);
        if (!(parsed instanceof Map<?, ?> values))
            throw new IOException("GWAS Catalog response is not a JSON object");
        Map<String, Object> result = new java.util.LinkedHashMap<>();
        for (Map.Entry<?, ?> entry : values.entrySet())
            if (entry.getKey() instanceof String key)
                result.put(key, entry.getValue());
        return result;
    }

    private static String requiredString(Map<?, ?> values, String key)
            throws IOException {
        String result = string(values.get(key));
        if (result.isBlank()) throw new IOException(
            "GWAS Catalog response is missing " + key);
        return result;
    }

    private static String string(Object value) {
        return value instanceof String text ? text : "";
    }

    private static long number(Object value) {
        return value instanceof Number number ? number.longValue() : -1L;
    }

    private static List<String> strings(Object value) {
        if (!(value instanceof List<?> values)) return List.of();
        return values.stream().filter(String.class::isInstance)
            .map(String.class::cast).toList();
    }

    private static List<String> nestedStrings(Object value, String key) {
        if (!(value instanceof List<?> values)) return List.of();
        List<String> result = new ArrayList<>();
        for (Object item : values) if (item instanceof Map<?, ?> map) {
            String text = string(map.get(key));
            if (!text.isBlank()) result.add(text);
        }
        return List.copyOf(result);
    }

    private static URI optionalUri(Object value) {
        String text = string(value);
        return text.isBlank() ? null : URI.create(text);
    }

    private static String encode(String value) {
        return URLEncoder.encode(value, StandardCharsets.UTF_8);
    }
}
