/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.*;
import java.net.*;
import java.net.http.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.security.*;
import java.time.Duration;
import java.util.*;
import java.util.concurrent.TimeUnit;

/** Strict options, atomic result directories and audited external operations. */
final class FollowupSupport {
    private FollowupSupport() { }
    static Map<String,String> options(String[] args,Set<String> allowed) {
        Map<String,String> map=new LinkedHashMap<>();
        for(int i=0;i<args.length;i++) {
            String key=args[i];
            if(!key.startsWith("--") || !allowed.contains(key.substring(2))) throw new IllegalArgumentException("Unknown option: "+key);
            if(map.containsKey(key.substring(2))) throw new IllegalArgumentException("Duplicate option: "+key);
            if(i+1==args.length || args[i+1].startsWith("--")) throw new IllegalArgumentException(key+" requires a value");
            map.put(key.substring(2),args[++i]);
        }
        return map;
    }
    static String required(Map<String,String> o,String key) {
        String value=o.get(key); if(value==null || value.isBlank()) throw new IllegalArgumentException("--"+key+" is required"); return value;
    }
    static double number(Map<String,String> o,String key,double fallback) {
        double value=Double.parseDouble(o.getOrDefault(key,Double.toString(fallback)));
        if(!Double.isFinite(value)) throw new IllegalArgumentException(key+" must be finite"); return value;
    }
    static int integer(Map<String,String> o,String key,int fallback,int min,int max) {
        int value=Integer.parseInt(o.getOrDefault(key,Integer.toString(fallback)));
        if(value<min || value>max) throw new IllegalArgumentException(key+" must be in ["+min+", "+max+"]"); return value;
    }
    static boolean bool(Map<String,String> o,String key,boolean fallback) {
        String value=o.getOrDefault(key,Boolean.toString(fallback));
        if(!value.equals("true") && !value.equals("false")) throw new IllegalArgumentException(key+" must be true or false"); return Boolean.parseBoolean(value);
    }
    static String hash(Path path) throws IOException {
        try {
            MessageDigest digest=MessageDigest.getInstance("SHA-256");
            try(InputStream in=Files.newInputStream(path)) { byte[] buffer=new byte[65536]; int n; while((n=in.read(buffer))>=0)digest.update(buffer,0,n); }
            return HexFormat.of().formatHex(digest.digest());
        } catch(NoSuchAlgorithmException impossible) { throw new AssertionError(impossible); }
    }
    interface Job { void execute(Path directory, Map<String,Object> manifest) throws Exception; }
    static String executable(String value) {
        return value.contains("/")||value.contains("\\")||Files.exists(Path.of(value))?Path.of(value).toAbsolutePath().normalize().toString():value;
    }
    static void transaction(Map<String,String> options,String operation,Job job) throws Exception {
        Path destination=Path.of(required(options,"out")).toAbsolutePath().normalize();
        if(Files.exists(destination)) throw new IOException("Output already exists: "+destination);
        Files.createDirectories(destination.getParent());
        Path staging=Files.createTempDirectory(destination.getParent(),".jlinalg-followup-");
        try {
            Map<String,Object> manifest=new LinkedHashMap<>();
            manifest.put("schema_version",1); manifest.put("operation",operation); manifest.put("jlinalg_version",JLinAlgCli.version());
            manifest.put("started_at",java.time.Instant.now().toString()); manifest.put("options",options);
            Map<String,String> hashes=new LinkedHashMap<>();
            for(String key:List.of("input","counts","cells","samples","features","feature-list","gene-sets","edges","hits","matrix","matrix-b","pheno","regulators","activities","source-file","script","fasta","annotations","components","reference-matrix","executable","rscript","perl","solver-path")) {
                if(options.containsKey(key) && Files.isRegularFile(Path.of(options.get(key)))) hashes.put(key,hash(Path.of(options.get(key))));
            }
            manifest.put("input_sha256",hashes);
            job.execute(staging,manifest);
            Map<String,String> outputs=new TreeMap<>();
            try(var paths=Files.walk(staging)) { for(Path p:paths.filter(Files::isRegularFile).toList()) outputs.put(staging.relativize(p).toString(),hash(p)); }
            manifest.put("output_sha256",outputs); manifest.put("completed_at",java.time.Instant.now().toString());
            Files.writeString(staging.resolve("manifest.yaml"),new org.yaml.snakeyaml.Yaml().dump(manifest));
            Files.move(staging,destination);
        } finally {
            if(Files.exists(staging)) try(var paths=Files.walk(staging)) { for(Path p:paths.sorted(Comparator.reverseOrder()).toList())Files.deleteIfExists(p); }
        }
    }
    static void table(Path file,List<String> header,List<String[]> rows) throws IOException {
        try(var out=Files.newBufferedWriter(file)) {
            row(out,header.toArray(String[]::new)); for(String[] row:rows) { if(row.length!=header.size())throw new IOException("Internal table width mismatch"); row(out,row); }
        }
    }
    static void row(Writer out,String[] fields) throws IOException {
        for(int i=0;i<fields.length;i++) {
            if(i>0)out.write('\t'); String value=fields[i];
            if(value.indexOf('\n')>=0 || value.indexOf('\r')>=0) throw new IOException("Multiline table field is unsupported");
            out.write(value.indexOf('\t')>=0 || value.indexOf('"')>=0 ? '"'+value.replace("\"","\"\"")+'"' : value);
        }
        out.write('\n');
    }
    static double finite(String value,String name) {
        double result=Double.parseDouble(value); if(!Double.isFinite(result))throw new IllegalArgumentException("Nonfinite "+name); return result;
    }
    static void request(String url,byte[] body,String type,Path output,Map<String,String> options) throws Exception {
        URI uri=URI.create(url);
        if(!Set.of("https","http").contains(uri.getScheme()) || uri.getUserInfo()!=null) throw new IllegalArgumentException("Require HTTP(S) URL without embedded credentials");
        int timeout=integer(options,"timeout",120,1,86400);
        HttpRequest.Builder builder=HttpRequest.newBuilder(uri).timeout(Duration.ofSeconds(timeout)).header("Accept",type);
        if(options.containsKey("token-env")) {
            String token=System.getenv(required(options,"token-env"));
            if(token==null || token.isBlank()) throw new IllegalArgumentException("Authentication environment variable is unset");
            builder.header("Authorization","Bearer "+token);
        }
        if(body!=null)builder.header("Content-Type",type).POST(HttpRequest.BodyPublishers.ofByteArray(body)); else builder.GET();
        HttpClient client=HttpClient.newBuilder().connectTimeout(Duration.ofSeconds(timeout)).followRedirects(HttpClient.Redirect.NEVER).build();
        HttpResponse<InputStream> response=client.send(builder.build(),HttpResponse.BodyHandlers.ofInputStream());
        try(InputStream in=response.body()) {
            if(response.statusCode()!=200) throw new IOException("HTTP request failed: status="+response.statusCode()+" (redirects require an explicit final URL)");
            long max=(long)integer(options,"max-download-mb",1024,1,1000000)*1024*1024;
            try(OutputStream out=Files.newOutputStream(output)) { byte[] buffer=new byte[65536]; long total=0; int n; while((n=in.read(buffer))>=0) { total+=n; if(total>max)throw new IOException("Response exceeds max-download-mb"); out.write(buffer,0,n); } }
        }
    }
    static void process(List<String> command,Path dir,int timeout) throws Exception {
        Process process=new ProcessBuilder(command).directory(dir.toFile()).redirectErrorStream(true).redirectOutput(dir.resolve("external.log").toFile()).start();
        try {
            if(!process.waitFor(timeout,TimeUnit.SECONDS)) { process.descendants().forEach(ProcessHandle::destroyForcibly); process.destroyForcibly(); throw new IOException("External tool timed out after "+timeout+" seconds"); }
            if(process.exitValue()!=0) {
                String log=Files.readString(dir.resolve("external.log"));
                throw new IOException("External tool failed ("+process.exitValue()+"): "+log.substring(Math.max(0,log.length()-5000)));
            }
        } catch(InterruptedException e) { process.destroyForcibly(); Thread.currentThread().interrupt(); throw e; }
    }
}
