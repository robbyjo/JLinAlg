/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.PrintStream;
import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HexFormat;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.regex.Pattern;

/** Explicit, versioned downloads; incomplete packages are never published as installed. */
final class EnrichmentDownloads {
    private final HttpClient client=HttpClient.newBuilder().connectTimeout(Duration.ofSeconds(30))
        .version(HttpClient.Version.HTTP_1_1).followRedirects(HttpClient.Redirect.NORMAL).build();
    private final PrintStream console;
    private final long maximumBytes;
    private final Properties manifest=new Properties();
    private Path stage;
    EnrichmentDownloads(PrintStream console,long maximumBytes) {
        if(maximumBytes<1) throw new IllegalArgumentException("maximum download size must be positive");
        this.console=console; this.maximumBytes=maximumBytes;
    }
    static String choices() {
        return "GO[:BP|MF|CC], Reactome, HPO, GWASCatalog, Monarch:Mondo, OpenTargets, KEGG, MSigDB:<collection>, Custom\n"
            + "MSigDB: use --source-file registered-download.gmt (or an accessible --source-url) and --gene-id-type.\n"
            + "KEGG: academic REST download requires --kegg-academic true; otherwise import local data.\n"
            + "All providers accept local GMT/term2gene files via --db-format and --source-file.\n";
    }
    void install(String type,Path target,Map<String,String> options) throws IOException,InterruptedException {
        String provider=EnrichmentDatabase.provider(type);
        Path destination=target.toAbsolutePath().normalize();
        if(Files.exists(destination)) throw new IOException("download destination already exists; choose a new version directory: "+destination);
        Files.createDirectories(destination.getParent());
        stage=Files.createTempDirectory(destination.getParent(),".enrichment-download-").toAbsolutePath().normalize();
        boolean complete=false;
        try {
            manifest.setProperty("format_version","1"); manifest.setProperty("provider",provider);
            manifest.setProperty("collection",type); manifest.setProperty("downloaded_at",Instant.now().toString());
            String species=options.getOrDefault("--species","human"), release=options.getOrDefault("--db-release","latest");
            if (Set.of("HPO", "GWASCatalog", "OpenTargets").contains(provider) && !species.equalsIgnoreCase("human"))
                throw new IOException(provider+" automatic source supports human only");
            if(!release.matches("[A-Za-z0-9._-]+")) throw new IOException("invalid database release");
            manifest.setProperty("species",species); manifest.setProperty("release",release);
            for(String key:List.of("--catalog-gene-source","--catalog-p-threshold","--trait-id-type","--association-min-score","--monarch-evidence","--evidence","--propagate"))
                if(options.containsKey(key)) manifest.setProperty("option."+key,options.get(key));
            if(options.containsKey("--source-file") || options.containsKey("--source-url")) {
                if(options.containsKey("--source-file") && options.containsKey("--source-url")) throw new IOException("choose --source-file or --source-url");
                String source=options.getOrDefault("--source-file",options.get("--source-url"));
                String sourcePath=options.containsKey("--source-url") ? URI.create(source).getPath() : source;
                String format=options.getOrDefault("--db-format",sourcePath.matches("(?i).*\\.gmt(?:\\.gz)?$") ? "gmt" : "term2gene");
                String extension=format.equals("gmt") ? "gmt" : sourcePath.matches("(?i).*\\.csv(?:\\.gz)?$") ? "csv" : "tsv";
                String name="data."+extension+(sourcePath.toLowerCase(java.util.Locale.ROOT).endsWith(".gz")?".gz":"");
                if(options.containsKey("--source-file")) {
                    Path input=Path.of(source); if(Files.size(input)>maximumBytes) throw new IOException("source exceeds size limit");
                    Files.copy(input,stage.resolve(name)); manifest.setProperty("source."+name,input.toAbsolutePath().toString());
                    manifest.setProperty("sha256."+name,sha256(stage.resolve(name)));
                } else fetch(source,name);
                configure(format,name,options.getOrDefault("--gene-id-type","unspecified"));
                if(options.containsKey("--ontology-url")) { fetch(options.get("--ontology-url"),"ontology.obo"); manifest.setProperty("ontology_file","ontology.obo"); }
                else if(options.containsKey("--ontology")) {
                    Files.copy(Path.of(options.get("--ontology")),stage.resolve("ontology.obo"));
                    manifest.setProperty("sha256.ontology.obo",sha256(stage.resolve("ontology.obo")));
                    manifest.setProperty("ontology_file","ontology.obo");
                }
            } else switch(provider) {
                case "GO" -> {
                    String code=switch(species.toLowerCase(java.util.Locale.ROOT)) { case "human" -> "HUMAN"; case "mouse" -> "MOUSE"; case "rat" -> "RAT"; default -> throw new IOException("GO automatic download supports human, mouse, rat; use --source-url for another organism"); };
                    // Ontology data-version is NOT the GO annotation pipeline release.
                    if(release.equals("latest")) {
                        var match=Pattern.compile("\\d{4}-\\d{2}-\\d{2}").matcher(text("https://current.geneontology.org/metadata/release-date.json"));
                        if(!match.find()) throw new IOException("GO release date missing");
                        release=match.group(); manifest.setProperty("release",release);
                    }
                    String base="https://release.geneontology.org/"+release+"/";
                    fetch(base+"ontology/go-basic.obo","ontology.obo");
                    fetch(base+"annotations/gaf/"+code+"-uniprot.gaf.gz","annotations.gaf.gz");
                    configure("gaf","annotations.gaf.gz",options.getOrDefault("--gene-id-type","symbol"));
                    manifest.setProperty("ontology_file","ontology.obo");
                }
                case "Reactome" -> {
                    String ids=options.getOrDefault("--gene-id-type","ensembl").toLowerCase(java.util.Locale.ROOT);
                    String file=switch(ids) { case "ensembl" -> "Ensembl2Reactome_All_Levels.txt"; case "entrez" -> "NCBI2Reactome_All_Levels.txt"; case "uniprot" -> "UniProt2Reactome_All_Levels.txt"; default -> throw new IOException("Reactome gene IDs must be ensembl, entrez, or uniprot"); };
                    String base="https://reactome.org/download/"+(release.equals("latest")?"current":release)+"/";
                    fetch(base+file,"mapping.tsv"); fetch(base+"ReactomePathwaysRelation.txt","parents.tsv");
                    configure("reactome","mapping.tsv",ids); manifest.setProperty("parents_file","parents.tsv");
                }
                case "HPO" -> {
                    release=githubRelease("obophenotype/human-phenotype-ontology",release);
                    String base="https://github.com/obophenotype/human-phenotype-ontology/releases/download/"+release+"/";
                    fetch(base+"genes_to_phenotype.txt","hpo.tsv"); fetch(base+"hp.obo","ontology.obo");
                    configure("hpo","hpo.tsv","entrez"); manifest.setProperty("ontology_file","ontology.obo");
                }
                case "GWASCatalog" -> {
                    if(!release.equals("latest")) throw new IOException("historical GWAS Catalog snapshots require --source-url or --source-file with an extracted TSV");
                    fetch("https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip","catalog.zip");
                    extractSingleTable("catalog.zip","catalog.tsv"); configure("gwas","catalog.tsv","symbol");
                }
                case "Monarch:Mondo" -> {
                    String evidence=options.getOrDefault("--monarch-evidence","causal");
                    if(!Set.of("causal","correlated").contains(evidence)) throw new IOException("Monarch evidence must be causal or correlated");
                    manifest.setProperty("option.--monarch-evidence",evidence);
                    fetch("https://data.monarchinitiative.org/monarch-kg/"+release+"/tsv/all_associations/"+evidence+"_gene_to_disease_association.all.tsv.gz","monarch.tsv.gz");
                    String mondo=githubRelease("monarch-initiative/mondo","latest"); manifest.setProperty("mondo_release",mondo);
                    fetch("https://github.com/monarch-initiative/mondo/releases/download/"+mondo+"/mondo.obo","ontology.obo");
                    configure("monarch","monarch.tsv.gz","curie"); manifest.setProperty("ontology_file","ontology.obo");
                }
                case "OpenTargets" -> {
                    String root="https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/";
                    if(release.equals("latest")) {
                        var matcher=Pattern.compile("href=\"(\\d{2}\\.\\d{2})/\"").matcher(text(root));
                        List<String> versions=new ArrayList<>(); while(matcher.find()) versions.add(matcher.group(1));
                        release=versions.stream().max(String::compareTo).orElseThrow(()->new IOException("Open Targets release index unavailable"));
                    }
                    manifest.setProperty("release",release);
                    parquet(root+release+"/output/association_overall_direct/","associations");
                    parquet(root+release+"/output/disease/","disease");
                    configure("opentargets","associations","ensembl");
                }
                case "KEGG" -> {
                    if(!"true".equals(options.get("--kegg-academic"))) throw new IOException("KEGG public REST is for academic users: specify --kegg-academic true, or import an authorized local file");
                    if(!release.equals("latest")) throw new IOException("KEGG REST exposes current data; historical snapshots require local import");
                    String organism=switch(species.toLowerCase(java.util.Locale.ROOT)) { case "human" -> "hsa"; case "mouse" -> "mmu"; case "rat" -> "rno"; default -> species; };
                    if(!Set.of("hsa","mmu","rno").contains(organism)) throw new IOException("automatic KEGG Entrez mapping supports human, mouse, rat; import other organisms with their explicit namespace");
                    fetch("https://rest.kegg.jp/link/pathway/"+organism,"links.tsv");
                    Thread.sleep(400); fetch("https://rest.kegg.jp/list/pathway/"+organism,"names.tsv");
                    configure("kegg","links.tsv","entrez");
                }
                case "MSigDB" -> throw new IOException("MSigDB requires a registered download: use --source-file collection.gmt --gene-id-type symbol|entrez; https://www.gsea-msigdb.org/gsea/msigdb/index.jsp");
                default -> throw new IOException("Custom collections require --source-file or --source-url");
            }
            if(options.containsKey("--gene-id-type") && !options.get("--gene-id-type").equalsIgnoreCase(manifest.getProperty("gene_id_type")))
                throw new IOException("requested gene namespace is unavailable for this provider");
            try(var writer=Files.newBufferedWriter(stage.resolve("manifest.properties"))) { manifest.store(writer,"JLinAlg enrichment package; source hashes identify the exact snapshot"); }
            EnrichmentDatabase check=EnrichmentDatabase.load(stage,type,options);
            console.println("Validated "+check.genes.size()+" gene sets; namespace="+manifest.getProperty("gene_id_type"));
            Files.move(stage,destination); complete=true;
            console.println("Installed enrichment database: "+destination);
        } finally {
            if(!complete && Files.exists(stage)) {
                // Only remove our uniquely-created sibling staging directory.
                if(!stage.getParent().equals(destination.getParent()) || !stage.getFileName().toString().startsWith(".enrichment-download-")) throw new IOException("unsafe staging cleanup");
                try(var paths=Files.walk(stage)) { for(Path p:paths.sorted(Comparator.reverseOrder()).toList()) Files.deleteIfExists(p); }
            }
        }
    }
    private void configure(String format,String file,String ids) { manifest.setProperty("data_format",format); manifest.setProperty("data_file",file); manifest.setProperty("gene_id_type",ids); }
    private void extractSingleTable(String archive,String name) throws IOException {
        int tables=0;
        try(var zip=new java.util.zip.ZipInputStream(Files.newInputStream(stage.resolve(archive)))) {
            for(java.util.zip.ZipEntry entry; (entry=zip.getNextEntry())!=null;) {
                if(entry.isDirectory()) continue;
                if(++tables!=1 || !entry.getName().endsWith(".tsv")) throw new IOException("Catalog archive must contain exactly one TSV");
                // Never use an archive-supplied path as an extraction destination.
                try(var output=Files.newOutputStream(stage.resolve(name))) {
                    byte[] buffer=new byte[65536]; long count=0;
                    for(int n; (n=zip.read(buffer))!=-1;) {
                        count+=n; if(count>maximumBytes) throw new IOException("expanded Catalog table exceeds --max-download-bytes");
                        output.write(buffer,0,n);
                    }
                    if(count==0) throw new IOException("empty Catalog table");
                }
                manifest.setProperty("archive_entry."+name,entry.getName());
            }
        }
        if(tables!=1) throw new IOException("Catalog archive has no table");
        manifest.setProperty("sha256."+name,sha256(stage.resolve(name)));
    }
    private String githubRelease(String repository,String requested) throws IOException,InterruptedException {
        if(!requested.equals("latest")) return requested;
        Object json=SimpleJson.parse(text("https://api.github.com/repos/"+repository+"/releases/latest"));
        if(!(json instanceof Map<?,?> map) || !(map.get("tag_name") instanceof String tag) || !tag.matches("[A-Za-z0-9._-]+")) throw new IOException("GitHub release metadata missing");
        if(!repository.endsWith("/mondo")) manifest.setProperty("release",tag);
        return tag;
    }
    private void parquet(String url,String directory) throws IOException,InterruptedException {
        var matcher=Pattern.compile("href=\"([^\"/]+\\.parquet)\"").matcher(text(url));
        int count=0;
        while(matcher.find()) { String name=matcher.group(1); if(!name.matches("[A-Za-z0-9._-]+")) throw new IOException("unsafe Parquet filename"); fetch(url+name,directory+"/"+name); count++; }
        if(count==0) throw new IOException("no Parquet parts found at "+url);
    }
    private HttpRequest request(String url) throws IOException {
        URI uri=URI.create(url); if(!"https".equalsIgnoreCase(uri.getScheme()) || uri.getHost()==null || uri.getUserInfo()!=null) throw new IOException("database source must be HTTPS without embedded credentials");
        return HttpRequest.newBuilder(uri).timeout(Duration.ofMinutes(30)).header("User-Agent","JLinAlg-enrichment/1").GET().build();
    }
    private String text(String url) throws IOException,InterruptedException {
        var response=client.send(request(url),HttpResponse.BodyHandlers.ofInputStream());
        try(var input=response.body()) {
            if(response.statusCode()!=200) throw new IOException("HTTP "+response.statusCode()+" from "+url);
            byte[] bytes=input.readNBytes(4*1024*1024+1); if(bytes.length>4*1024*1024) throw new IOException("database index too large");
            return new String(bytes,java.nio.charset.StandardCharsets.UTF_8);
        }
    }
    private void fetch(String url,String name) throws IOException,InterruptedException {
        for (int attempt=1; ; attempt++) {
            try { fetchOnce(url,name); return; }
            catch (IOException e) {
                String message=String.valueOf(e.getMessage());
                if (attempt==3 || message.contains("exceeds") || message.startsWith("HTTP 4")) throw e;
                console.println("Retry "+attempt+" after source download failure: "+e.getMessage());
                Thread.sleep(attempt*1000L);
            }
        }
    }
    private void fetchOnce(String url,String name) throws IOException,InterruptedException {
        console.println("Downloading "+url);
        Path output=EnrichmentDatabase.safeChild(stage,name); Files.createDirectories(output.getParent());
        var response=client.send(request(url),HttpResponse.BodyHandlers.ofInputStream());
        try(var input=response.body()) {
            if(response.statusCode()!=200) throw new IOException("HTTP "+response.statusCode()+" downloading "+url);
            long declared=response.headers().firstValueAsLong("Content-Length").orElse(-1);
            if(declared>maximumBytes) throw new IOException("download exceeds --max-download-bytes");
            long total=0,next=128L*1024*1024;
            try(var file=Files.newOutputStream(output)) {
                byte[] buffer=new byte[65536];
                for(int n;(n=input.read(buffer))!=-1;) {
                    total+=n; if(total>maximumBytes) throw new IOException("download exceeds --max-download-bytes");
                    file.write(buffer,0,n);
                    if(total>=next) { console.println("  "+name+": "+total+" bytes"); next+=128L*1024*1024; }
                }
            }
            if(total==0 || (declared>=0 && total!=declared)) throw new IOException("empty or incomplete download: "+url);
        }
        manifest.setProperty("source."+name,response.uri().toString());
        manifest.setProperty("sha256."+name,sha256(output));
        response.headers().firstValue("Last-Modified").ifPresent(v->manifest.setProperty("last_modified."+name,v));
    }
    static String sha256(Path path) throws IOException {
        try {
            MessageDigest digest=MessageDigest.getInstance("SHA-256");
            try(var input=Files.newInputStream(path)) { byte[] b=new byte[65536]; for(int n;(n=input.read(b))!=-1;) digest.update(b,0,n); }
            return HexFormat.of().formatHex(digest.digest());
        } catch(NoSuchAlgorithmException e) { throw new IllegalStateException(e); }
    }
}
