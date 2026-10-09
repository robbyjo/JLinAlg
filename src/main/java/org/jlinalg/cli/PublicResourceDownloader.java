/*
 * Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package org.jlinalg.cli;

import java.io.BufferedInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InterruptedIOException;
import java.io.OutputStream;
import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.charset.StandardCharsets;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.time.Duration;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import org.apache.hc.client5.http.classic.methods.HttpGet;
import org.apache.hc.client5.http.config.ConnectionConfig;
import org.apache.hc.client5.http.config.RequestConfig;
import org.apache.hc.client5.http.impl.classic.CloseableHttpClient;
import org.apache.hc.client5.http.impl.classic.HttpClients;
import org.apache.hc.client5.http.impl.io.PoolingHttpClientConnectionManagerBuilder;
import org.apache.hc.core5.http.HttpEntity;
import org.apache.hc.core5.util.Timeout;

/** Bounded HTTPS downloads used only for public reference resources. */
final class PublicResourceDownloader {
    record Result(String source, String lastModified, long bytes) { }
    record Response(int status, String body) { }

    private PublicResourceDownloader() { }

    static Result download(URI source, Path output, long maximumBytes,
            String userAgent, java.io.PrintStream progress,
            long progressInterval) throws IOException, InterruptedException {
        return download(source, output, maximumBytes, userAgent, progress, progressInterval, true);
    }

    static Result downloadStrict(URI source, Path output, long maximumBytes) throws IOException, InterruptedException {
        StudyPublicEvidence.validate(source);
        return download(source, output, maximumBytes, "JLinAlg-study/1", null, 0, false);
    }

    /** Stream one complete immutable provider file with no regional Range request. */
    static Result downloadComplete(URI source, Path output,
            long maximumBytes) throws IOException, InterruptedException {
        StudyPublicEvidence.validate(source);
        validate(source, maximumBytes);
        HttpClient client = HttpClient.newBuilder()
            .connectTimeout(Duration.ofSeconds(30))
            .followRedirects(HttpClient.Redirect.NEVER)
            .version(HttpClient.Version.HTTP_1_1).build();
        HttpRequest head = HttpRequest.newBuilder(source)
            .timeout(Duration.ofMinutes(2))
            .header("User-Agent", "JLinAlg-public-resource/1")
            .header("Accept-Encoding", "identity")
            .method("HEAD", HttpRequest.BodyPublishers.noBody()).build();
        HttpResponse<Void> metadata = client.send(head,
            HttpResponse.BodyHandlers.discarding());
        if (metadata.statusCode() < 200 || metadata.statusCode() >= 300)
            throw new IOException("HTTP " + metadata.statusCode()
                + " inspecting complete public resource " + source);
        long declared = metadata.headers().firstValueAsLong(
            "Content-Length").orElse(-1);
        if (declared < 1)
            throw new IOException("Public resource did not declare a positive Content-Length: " + source);
        if (declared > maximumBytes)
            throw new IOException("Download exceeds safety limit: " + source);
        try (FileChannel file = FileChannel.open(output,
                StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING,
                StandardOpenOption.WRITE)) {
            file.write(ByteBuffer.wrap(new byte[] {0}), declared - 1);
        }
        final long chunkBytes = 32L * 1024 * 1024;
        ExecutorService workers = Executors.newFixedThreadPool(4);
        List<Future<Long>> tasks = new ArrayList<>();
        try {
            for (long start = 0; start < declared; start += chunkBytes) {
                long first = start;
                long last = Math.min(declared - 1, first + chunkBytes - 1);
                tasks.add(workers.submit(() -> downloadRange(client, source,
                    output, first, last, declared)));
            }
            long copied = 0;
            for (Future<Long> task : tasks) try {
                copied += task.get();
            } catch (ExecutionException failure) {
                Throwable cause = failure.getCause();
                if (cause instanceof IOException io) throw io;
                if (cause instanceof InterruptedException interrupted)
                    throw interrupted;
                throw new IOException("Complete public-resource chunk failed",
                    cause);
            }
            if (copied != declared || Files.size(output) != declared)
                throw new IOException("Empty or incomplete download: " + source);
        } finally {
            tasks.forEach(task -> task.cancel(true));
            workers.shutdownNow();
        }
        if (Files.size(output) != declared)
            throw new IOException("Empty or incomplete download: " + source);
        return new Result(source.toString(), metadata.headers()
            .firstValue("Last-Modified").orElse(null), declared);
    }

    private static long downloadRange(HttpClient client, URI source,
            Path output, long first, long last, long total)
            throws IOException, InterruptedException {
        HttpRequest request = HttpRequest.newBuilder(source)
            .timeout(Duration.ofHours(24))
            .header("User-Agent", "JLinAlg-public-resource/1")
            .header("Accept-Encoding", "identity")
            .header("Range", "bytes=" + first + "-" + last).GET().build();
        HttpResponse<InputStream> response = client.send(request,
            HttpResponse.BodyHandlers.ofInputStream());
        if (response.statusCode() != 206) {
            response.body().close();
            throw new IOException("Server did not honor complete-file byte range "
                + first + "-" + last + " for " + source + ": HTTP "
                + response.statusCode());
        }
        String expectedRange = "bytes " + first + "-" + last + "/" + total;
        String actualRange = response.headers().firstValue("Content-Range")
            .orElse("");
        if (!expectedRange.equals(actualRange)) {
            response.body().close();
            throw new IOException("Unexpected Content-Range for " + source
                + ": " + actualRange);
        }
        long expected = last - first + 1, copied = 0;
        try (InputStream input = new BufferedInputStream(response.body());
                FileChannel file = FileChannel.open(output,
                    StandardOpenOption.WRITE)) {
            byte[] bytes = new byte[1024 * 1024];
            for (int count; (count = input.read(bytes)) >= 0;) {
                if (Thread.currentThread().isInterrupted())
                    throw new InterruptedIOException(
                        "complete public-resource download interrupted");
                if (count == 0) continue;
                if (copied + count > expected)
                    throw new IOException("Byte range exceeded declared size for "
                        + source);
                ByteBuffer buffer = ByteBuffer.wrap(bytes, 0, count);
                long position = first + copied;
                while (buffer.hasRemaining())
                    position += file.write(buffer, position);
                copied += count;
            }
        }
        if (copied != expected)
            throw new IOException("Incomplete byte range " + first + "-"
                + last + " for " + source);
        return copied;
    }

    /** Preserve non-success bodies for the evidence audit without following redirects. */
    static Response responseStrict(URI source, long maximumBytes) throws IOException, InterruptedException {
        StudyPublicEvidence.validate(source);
        validate(source, maximumBytes);
        HttpGet request = request(source, "JLinAlg-study/1");
        request.addHeader("Accept", "application/json");
        try (CloseableHttpClient client = client(false)) {
            try {
                return client.execute(request, response -> {
                    HttpEntity entity = response.getEntity();
                    if (entity == null) return new Response(response.getCode(), "");
                    if (entity.getContentLength() > maximumBytes)
                        throw new IOException("Response exceeds safety limit: " + source);
                    ByteArrayOutputStream output = new ByteArrayOutputStream();
                    try (InputStream input = entity.getContent()) {
                        copy(input, output, maximumBytes, source, null, 0);
                    }
                    return new Response(response.getCode(), output.toString(StandardCharsets.UTF_8));
                });
            } catch (InterruptedIOException interrupted) {
                // A socket timeout is an unavailable service, not a cancelled study run.
                if (Thread.currentThread().isInterrupted())
                    throw new InterruptedException("public evidence request interrupted");
                throw interrupted;
            }
        }
    }

    private static Result download(URI source, Path output, long maximumBytes,
            String userAgent, java.io.PrintStream progress, long progressInterval,
            boolean redirects) throws IOException, InterruptedException {
        validate(source, maximumBytes);
        HttpGet request = request(source, userAgent);
        try (CloseableHttpClient client = client(redirects)) {
            try {
                return client.execute(request, response -> {
                    int code = response.getCode();
                    if (code < 200 || code >= 300)
                        throw new IOException("HTTP " + code
                            + " downloading " + source);
                    HttpEntity entity = response.getEntity();
                    if (entity == null)
                        throw new IOException("Empty HTTP response downloading "
                            + source);
                    long declared = entity.getContentLength();
                    if (declared > maximumBytes)
                        throw new IOException("Download exceeds safety limit: "
                            + source);
                    long copied;
                    try (InputStream input = new BufferedInputStream(
                                entity.getContent());
                            OutputStream file = Files.newOutputStream(output,
                                StandardOpenOption.CREATE,
                                StandardOpenOption.TRUNCATE_EXISTING)) {
                        copied = copy(input, file, maximumBytes, source,
                            progress, progressInterval);
                    }
                    if (copied == 0 || declared >= 0 && copied != declared)
                        throw new IOException("Empty or incomplete download: "
                            + source);
                    var modified = response.getFirstHeader("Last-Modified");
                    return new Result(source.toString(),
                        modified == null ? null : modified.getValue(), copied);
                });
            } catch (InterruptedIOException interrupted) {
                Thread.currentThread().interrupt();
                throw new InterruptedException(
                    "public reference download interrupted");
            }
        }
    }

    static String text(URI source, long maximumBytes, String userAgent)
            throws IOException, InterruptedException {
        validate(source, maximumBytes);
        HttpGet request = request(source, userAgent);
        try (CloseableHttpClient client = client()) {
            try {
                return client.execute(request, response -> {
                    int code = response.getCode();
                    if (code < 200 || code >= 300)
                        throw new IOException("HTTP " + code + " from "
                            + source);
                    HttpEntity entity = response.getEntity();
                    if (entity == null || entity.getContentLength() > maximumBytes)
                        throw new IOException("Empty or oversized response from "
                            + source);
                    ByteArrayOutputStream output = new ByteArrayOutputStream();
                    try (InputStream input = entity.getContent()) {
                        copy(input, output, maximumBytes, source, null, 0);
                    }
                    return output.toString(StandardCharsets.UTF_8);
                });
            } catch (InterruptedIOException interrupted) {
                Thread.currentThread().interrupt();
                throw new InterruptedException(
                    "public reference download interrupted");
            }
        }
    }

    private static HttpGet request(URI source, String userAgent) {
        HttpGet request = new HttpGet(source);
        request.addHeader("User-Agent", userAgent);
        request.addHeader("Accept-Encoding", "identity");
        return request;
    }

    private static CloseableHttpClient client() {
        return client(true);
    }

    private static CloseableHttpClient client(boolean redirects) {
        RequestConfig configuration = RequestConfig.custom()
            .setConnectionRequestTimeout(Timeout.ofSeconds(30))
            .setResponseTimeout(Timeout.ofMinutes(2))
            .setRedirectsEnabled(redirects)
            .build();
        var connections = PoolingHttpClientConnectionManagerBuilder.create()
            .setDefaultConnectionConfig(ConnectionConfig.custom()
                .setConnectTimeout(Timeout.ofSeconds(30))
                .setSocketTimeout(Timeout.ofMinutes(2))
                .build())
            .build();
        var builder = HttpClients.custom()
            .setConnectionManager(connections)
            .setDefaultRequestConfig(configuration)
            .disableContentCompression();
        // Each strict acquisition is one auditable attempt; retries belong to resume.
        if (!redirects) builder.disableAutomaticRetries();
        return builder.build();
    }

    private static long copy(InputStream input, OutputStream output,
            long maximumBytes, URI source, java.io.PrintStream progress,
            long progressInterval) throws IOException {
        long copied = 0;
        long nextProgress = progressInterval;
        byte[] buffer = new byte[1024 * 1024];
        for (int count; (count = input.read(buffer)) >= 0;) {
            if (Thread.currentThread().isInterrupted())
                throw new InterruptedIOException(
                    "public reference download interrupted");
            if (count == 0) continue;
            copied += count;
            if (copied > maximumBytes)
                throw new IOException("Download exceeds safety limit: "
                    + source);
            output.write(buffer, 0, count);
            if (progress != null && progressInterval > 0
                    && copied >= nextProgress) {
                progress.println("  downloaded " + humanBytes(copied));
                nextProgress += progressInterval;
            }
        }
        return copied;
    }

    private static void validate(URI source, long maximumBytes)
            throws IOException {
        if (!"https".equalsIgnoreCase(source.getScheme())
                || source.getHost() == null || source.getUserInfo() != null)
            throw new IOException(
                "public reference source must be HTTPS without credentials");
        if (maximumBytes < 1)
            throw new IOException("download safety limit must be positive");
    }

    private static String humanBytes(long bytes) {
        if (bytes >= 1024L * 1024 * 1024)
            return String.format(java.util.Locale.ROOT, "%.1f GiB",
                bytes / (1024.0 * 1024 * 1024));
        if (bytes >= 1024L * 1024)
            return String.format(java.util.Locale.ROOT, "%.1f MiB",
                bytes / (1024.0 * 1024));
        if (bytes >= 1024L)
            return String.format(java.util.Locale.ROOT, "%.1f KiB",
                bytes / 1024.0);
        return bytes + " B";
    }
}
