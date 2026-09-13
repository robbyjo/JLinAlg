/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import org.jlinalg.enrichment.GeneSet;

/** Source-specific annotation parsing and ontology closure. No gene-list dependent filtering. */
final class EnrichmentDatabase {
    final Map<String, Set<String>> genes = new TreeMap<>();
    final Map<String, String> names = new HashMap<>(), aspects = new HashMap<>(), aliases = new HashMap<>();
    final Map<String, Set<String>> parents = new TreeMap<>();
    final Set<String> obsolete = new HashSet<>();
    final Properties metadata = new Properties();
    long rejectedRows;

    static String provider(String type) {
        String t = type.toUpperCase(Locale.ROOT);
        if (t.equals("GO") || Set.of("GO:BP","GO:MF","GO:CC").contains(t)) return "GO";
        if (t.startsWith("MSIGDB:") || t.equals("MSIGDB")) return "MSigDB";
        return switch(t) { case "REACTOME" -> "Reactome"; case "KEGG" -> "KEGG";
            case "MONARCH:MONDO" -> "Monarch:Mondo"; case "HPO" -> "HPO";
            case "GWASCATALOG", "GWAS" -> "GWASCatalog"; case "OPENTARGETS" -> "OpenTargets";
            case "CUSTOM" -> "Custom"; default -> throw new IllegalArgumentException("unknown enrichment collection: " + type); };
    }

    static EnrichmentDatabase load(Path path, String type, Map<String,String> options) throws IOException {
        EnrichmentDatabase db = new EnrichmentDatabase();
        String provider = provider(type), format = options.get("--db-format");
        Path data = path, ontology = options.containsKey("--ontology") ? Path.of(options.get("--ontology")) : null;
        if (Files.isDirectory(path)) {
            try (var r = Files.newBufferedReader(path.resolve("manifest.properties"))) { db.metadata.load(r); }
            if (!"1".equals(db.metadata.getProperty("format_version"))) throw new IOException("unsupported enrichment package version");
            if (!provider.equals(db.metadata.getProperty("provider"))) throw new IOException("database provider does not match --enrichment");
            String collection = db.metadata.getProperty("collection");
            if (provider.equals("MSigDB") && !type.equalsIgnoreCase(collection))
                throw new IOException("MSigDB package collection does not match --enrichment");
            for (String key : db.metadata.stringPropertyNames()) if (key.startsWith("sha256.")) {
                Path file = safeChild(path, key.substring(7));
                if (!EnrichmentDownloads.sha256(file).equals(db.metadata.getProperty(key)))
                    throw new IOException("database integrity failure: " + file.getFileName());
            }
            if (format != null && !format.equals(db.metadata.getProperty("data_format"))) throw new IOException("--db-format conflicts with package");
            format = db.metadata.getProperty("data_format");
            data = safeChild(path, db.metadata.getProperty("data_file"));
            if (ontology == null && db.metadata.containsKey("ontology_file")) ontology = safeChild(path, db.metadata.getProperty("ontology_file"));
            if (options.containsKey("--gene-id-type") && !options.get("--gene-id-type").equalsIgnoreCase(db.metadata.getProperty("gene_id_type")))
                throw new IOException("gene identifier namespace does not match database package");
            if (options.containsKey("--species") && !options.get("--species").equalsIgnoreCase(db.metadata.getProperty("species")))
                throw new IOException("species does not match database package; install a separate package for that species");
            if (options.containsKey("--monarch-evidence") && !options.get("--monarch-evidence").equals(db.metadata.getProperty("option.--monarch-evidence")))
                throw new IOException("Monarch evidence is fixed by the downloaded source; install a separate package");
        } else {
            db.metadata.setProperty("provider",provider);
            db.metadata.setProperty("sha256.input",EnrichmentDownloads.sha256(path));
            db.metadata.setProperty("gene_id_type",options.getOrDefault("--gene-id-type",switch(provider) {
                case "GO", "GWASCatalog" -> "symbol"; case "HPO", "KEGG" -> "entrez";
                case "Monarch:Mondo" -> "curie"; default -> "unspecified";
            }));
        }
        if (format == null) format = data.toString().matches("(?i).*\\.gmt(?:\\.gz)?$") ? "gmt" : switch(provider) {
            case "GWASCatalog" -> "gwas"; case "GO" -> "gaf"; case "HPO" -> "hpo";
            case "Reactome" -> "reactome"; case "Monarch:Mondo" -> "monarch";
            default -> "term2gene";
        };
        if (ontology != null) {
            db.obo(ontology);
            String prefix=switch(provider) { case "GO" -> "GO:"; case "HPO" -> "HP:"; case "Monarch:Mondo" -> "MONDO:"; default -> ""; };
            // Full OBO exports include imported ontologies and cross-domain relations.
            // Enrichment inheritance belongs only to the requested term hierarchy.
            if (!prefix.isEmpty()) {
                db.parents.keySet().removeIf(id -> !id.startsWith(prefix));
                db.parents.values().forEach(ps -> ps.removeIf(id -> !id.startsWith(prefix)));
            }
            db.metadata.setProperty("ontology_sha256",EnrichmentDownloads.sha256(ontology));
        }
        Map<String,String> settings = new HashMap<>();
        for (String key : db.metadata.stringPropertyNames()) if(key.startsWith("option.")) settings.put(key.substring(7),db.metadata.getProperty(key));
        settings.putAll(options);
        for (String key : List.of("--evidence", "--catalog-gene-source", "--catalog-p-threshold", "--trait-id-type", "--association-min-score"))
            if (settings.containsKey(key)) db.metadata.setProperty("effective."+key,settings.get(key));
        switch(format) {
            case "gmt" -> db.gmt(data);
            case "term2gene" -> db.term2gene(data);
            case "gaf" -> {
                if (ontology == null) throw new IOException("GO GAF requires --ontology or a downloaded GO package");
                db.gaf(data,settings);
            }
            case "gwas" -> db.gwas(data,settings);
            case "hpo" -> db.hpo(data);
            case "reactome" -> db.reactome(data,settings);
            case "monarch" -> db.monarch(data,settings);
            case "opentargets" -> db.openTargets(data,settings);
            case "kegg" -> db.kegg(data);
            default -> throw new IOException("unsupported database format: " + format);
        }
        if (Files.isDirectory(path) && db.metadata.containsKey("parents_file")) db.parentTable(safeChild(path,db.metadata.getProperty("parents_file")));
        boolean propagate = !settings.getOrDefault("--propagate","true").equals("false");
        if (propagate && !db.parents.isEmpty()) db.propagate();
        db.metadata.setProperty("propagation",propagate ? "is_a,part_of where applicable; source closure retained" : "none added; source closure retained");
        if (provider.equals("GO") && type.contains(":")) {
            String aspect = type.substring(type.indexOf(':')+1).toUpperCase(Locale.ROOT);
            db.genes.keySet().removeIf(id -> !aspect.equals(db.aspects.get(id)));
        }
        // Exclude uninformative ontology roots by ID, independently of observed hits.
        for(String root : List.of("GO:0008150","GO:0003674","GO:0005575","HP:0000001","MONDO:0000001")) db.genes.remove(root);
        db.metadata.setProperty("rejected_source_rows",Long.toString(db.rejectedRows));
        if (db.genes.isEmpty()) throw new IOException("database has no usable gene sets for the requested settings");
        return db;
    }
    static Path safeChild(Path directory,String name) throws IOException {
        if(name==null) throw new IOException("database manifest is incomplete");
        Path root=directory.toAbsolutePath().normalize(), child=root.resolve(name).normalize();
        if(!child.startsWith(root) || child.equals(root)) throw new IOException("unsafe database file path");
        return child;
    }
    List<GeneSet> sets() { return genes.entrySet().stream().map(e -> new GeneSet(e.getKey(),names.getOrDefault(e.getKey(),e.getKey()),e.getValue())).toList(); }
    void add(String id,String name,String gene) {
        if(missing(id) || missing(gene)) return;
        id=aliases.getOrDefault(id.trim(),id.trim());
        if(obsolete.contains(id)) return;
        genes.computeIfAbsent(id,ignored -> new TreeSet<>()).add(gene.trim());
        if(!missing(name)) names.putIfAbsent(id,name.trim());
    }
    static boolean missing(String s) { return EnrichmentExpression.missing(s) || s.equalsIgnoreCase("intergenic") || s.equalsIgnoreCase("NR"); }
    static BufferedReader reader(Path path) throws IOException {
        var input=Files.newInputStream(path);
        try { return new BufferedReader(new InputStreamReader(path.toString().endsWith(".gz") ? new GZIPInputStream(input) : input,StandardCharsets.UTF_8)); }
        catch(IOException e) { input.close(); throw e; }
    }
    @FunctionalInterface interface RowConsumer { void accept(Map<String,String> row) throws IOException; }
    static void table(Path path,RowConsumer consumer) throws IOException {
        table(path,consumer,false);
    }
    private static void table(Path path,RowConsumer consumer,boolean literalTsv) throws IOException {
        try(var r=reader(path)) {
            String first=r.readLine(); if(first==null) throw new IOException("empty table: "+path);
            first=first.replaceFirst("^\\uFEFF","");
            char delimiter=path.toString().matches("(?i).*\\.csv(?:\\.gz)?$") ? ',' : '\t';
            List<String> header=DelimitedData.parse(first,delimiter,1,path);
            Set<String> keys=new HashSet<>();
            for(int i=0;i<header.size();i++) { String key=header.get(i).trim(); header.set(i,key); if(key.isBlank() || !keys.add(key)) throw new IOException("duplicate or empty database column"); }
            long line=1;
            for(String text;(text=r.readLine())!=null;) {
                line++; if(text.isBlank()) continue;
                List<String> values=literalTsv && delimiter=='\t' ? List.of(text.split("\t",-1)) : DelimitedData.parse(text,delimiter,line,path);
                if(values.size()!=header.size()) throw new IOException("database row width mismatch at "+line);
                Map<String,String> row=new LinkedHashMap<>(); for(int i=0;i<header.size();i++) row.put(header.get(i),values.get(i).trim());
                consumer.accept(row);
            }
        }
    }
    static String field(Map<String,String> row,String... candidates) throws IOException {
        for(String candidate:candidates) if(row.containsKey(candidate)) return row.get(candidate);
        for(String candidate:candidates) for(var e:row.entrySet()) if(canonical(e.getKey()).equals(canonical(candidate))) return e.getValue();
        throw new IOException("database column missing: "+String.join(" / ",candidates));
    }
    static String optional(Map<String,String> row,String fallback,String... candidates) {
        try { return field(row,candidates); } catch(IOException e) { return fallback; }
    }
    private static String canonical(String x) { return x.toUpperCase(Locale.ROOT).replaceAll("[^A-Z0-9]",""); }
    private void term2gene(Path p) throws IOException {
        table(p,row -> add(field(row,"term_id"),optional(row,"","term_name"),field(row,"gene_id")));
    }
    private void gmt(Path p) throws IOException {
        try(var r=reader(p)) { for(String line;(line=r.readLine())!=null;) {
            if(line.isBlank() || line.startsWith("#")) continue;
            String[] x=line.split("\t",-1); if(x.length<3) throw new IOException("GMT needs term, description and members");
            for(int i=2;i<x.length;i++) add(x[0],x[1],x[i]);
        } }
    }
    private void gaf(Path p,Map<String,String> options) throws IOException {
        String idType=options.getOrDefault("--gene-id-type",metadata.getProperty("gene_id_type","symbol"));
        if(!Set.of("symbol","uniprot").contains(idType.toLowerCase(Locale.ROOT))) throw new IOException("GAF supports symbol or uniprot IDs");
        Set<String> evidence=options.containsKey("--evidence") ? Set.of(options.get("--evidence").split(",")) : Set.of();
        String species=options.getOrDefault("--species",metadata.getProperty("species","human"));
        String taxon=switch(species.toLowerCase(Locale.ROOT)) { case "human" -> "taxon:9606"; case "mouse" -> "taxon:10090"; case "rat" -> "taxon:10116"; default -> species; };
        try(var r=reader(p)) { for(String line;(line=r.readLine())!=null;) {
            if(line.startsWith("!")) { if(line.startsWith("!date-generated:") || line.startsWith("!go-version:")) metadata.setProperty(line.substring(1,line.indexOf(':')),line.substring(line.indexOf(':')+1).trim()); continue; }
            if(line.isBlank()) continue;
            String[] x=line.split("\t",-1); if(x.length<15) throw new IOException("GAF row has fewer than 15 fields");
            if (!x[12].split("\\|")[0].equals(taxon)) { rejectedRows++; continue; }
            if(List.of(x[3].split("\\|")).contains("NOT") || (!evidence.isEmpty() && !evidence.contains(x[6]))) { rejectedRows++; continue; }
            if(!names.containsKey(x[4]) && !aliases.containsKey(x[4])) { rejectedRows++; continue; }
            add(x[4],"",idType.equalsIgnoreCase("symbol") ? x[2] : x[1]);
        } }
    }
    private void gwas(Path p,Map<String,String> options) throws IOException {
        String source=options.getOrDefault("--catalog-gene-source","union").toLowerCase(Locale.ROOT);
        if(!Set.of("union","reported","mapped").contains(source)) throw new IOException("invalid catalog gene source");
        double threshold=Double.parseDouble(options.getOrDefault("--catalog-p-threshold","5e-8"));
        if(!(threshold>0 && threshold<=1)) throw new IOException("invalid catalog p threshold");
        String traits=options.getOrDefault("--trait-id-type","efo");
        if(!Set.of("efo","label").contains(traits)) throw new IOException("trait ID type must be efo or label");
        metadata.setProperty("catalog_gene_source",source); metadata.setProperty("catalog_p_threshold",Double.toString(threshold)); metadata.setProperty("trait_id_type",traits);
        table(p,row -> {
            String raw=field(row,"P-VALUE","P.VALUE"); double pv;
            try { pv=Double.parseDouble(raw); } catch(NumberFormatException e) { rejectedRows++; return; }
            if(!(pv>0 && pv<=1) || !Double.isFinite(pv)) { rejectedRows++; return; }
            if(pv>threshold) return;
            String label=field(row,"DISEASE/TRAIT","DISEASE.TRAIT"), term=traits.equals("efo") ? field(row,"MAPPED_TRAIT_URI") : label;
            if(missing(term)) { rejectedRows++; return; }
            Set<String> gs=new TreeSet<>();
            if(!source.equals("mapped")) splitGenes(gs,field(row,"REPORTED GENE(S)","REPORTED.GENE.S."));
            if(!source.equals("reported")) splitGenes(gs,field(row,"MAPPED_GENE","MAPPED GENE(S)"));
            // Use ontology IDs for each mapped trait; retain the study's original label as provenance in raw data.
            for(String t:term.split("[,;]\\s*")) {
                String id=t.trim().replaceFirst("^https?://[^ ]*/","").replaceFirst("^(EFO|MONDO|HP)_","$1:");
                String name=traits.equals("label") ? label : optional(row,label,"MAPPED_TRAIT");
                for(String gene:gs) add(id,name,gene);
            }
        },true); // Catalog TSV quotes are literal characters, including unmatched source quotes.
    }
    private static void splitGenes(Set<String> result,String raw) {
        for(String gene:raw.split("\\s+-\\s+|[,;]\\s*")) if(!missing(gene.trim())) result.add(gene.trim());
    }
    private void hpo(Path p) throws IOException {
        table(p,row -> add(field(row,"hpo_id"),field(row,"hpo_name"),field(row,"ncbi_gene_id")));
    }
    private void reactome(Path p,Map<String,String> options) throws IOException {
        String species=options.getOrDefault("--species",metadata.getProperty("species","human"));
        String organism=species.equalsIgnoreCase("human") ? "Homo sapiens" : species.equalsIgnoreCase("mouse") ? "Mus musculus" : species;
        try(var r=reader(p)) { for(String line;(line=r.readLine())!=null;) {
            String[] x=line.split("\t",-1); if(x.length<6) throw new IOException("Reactome mapping needs six fields");
            if(x[5].equalsIgnoreCase(organism)) add(x[1],x[3],x[0]);
        } }
    }
    private void monarch(Path p,Map<String,String> options) throws IOException {
        String species=options.getOrDefault("--species",metadata.getProperty("species","human"));
        String taxon=species.equalsIgnoreCase("human") ? "NCBITaxon:9606" : species.equalsIgnoreCase("mouse") ? "NCBITaxon:10090" : species;
        table(p,row -> {
            if(!field(row,"subject_taxon").equals(taxon)) return;
            String term=field(row,"object"); if(!term.startsWith("MONDO:")) { rejectedRows++; return; }
            add(term,optional(row,"","object_label"),field(row,"subject"));
        });
    }
    private void kegg(Path p) throws IOException {
        try(var r=reader(p)) { for(String line;(line=r.readLine())!=null;) {
            String[] x=line.split("\t"); if(x.length!=2) throw new IOException("KEGG link row needs two fields");
            String term=x[0].startsWith("path:") ? x[0] : x[1], gene=x[0].startsWith("path:") ? x[1] : x[0];
            add(term.replace("path:",""),"",gene.substring(gene.indexOf(':')+1));
        } }
        Path labels=p.resolveSibling("names.tsv");
        if(Files.exists(labels)) try(var r=reader(labels)) { for(String line;(line=r.readLine())!=null;) {
            String[] x=line.split("\t",2); if(x.length==2) names.put(x[0].replace("path:",""),x[1]);
        } }
    }
    private void openTargets(Path p,Map<String,String> options) throws IOException {
        double minimum=Double.parseDouble(options.getOrDefault("--association-min-score","0.5"));
        if(!(minimum>0 && minimum<=1)) throw new IOException("association minimum score must be in (0,1]");
        metadata.setProperty("association_min_score",Double.toString(minimum));
        String files=p.toAbsolutePath().toString().replace('\\','/')+"/*.parquet";
        try(var connection=java.sql.DriverManager.getConnection("jdbc:duckdb:")) {
            Set<String> columns=new HashSet<>();
            try(var schema=connection.prepareStatement("DESCRIBE SELECT * FROM read_parquet(?)")) {
                schema.setString(1,files);
                try(var rows=schema.executeQuery()) { while(rows.next()) columns.add(rows.getString(1)); }
            }
            String score=columns.contains("associationScore") ? "associationScore" : "score";
            if(!columns.containsAll(Set.of("targetId","diseaseId",score))) throw new IOException("Open Targets association schema lacks targetId, diseaseId or an overall score");
            metadata.setProperty("association_score_column",score);
            try(var query=connection.prepareStatement("SELECT targetId, diseaseId FROM read_parquet(?) WHERE "+score+" >= ?")) {
                query.setString(1,files); query.setDouble(2,minimum);
                try(var rows=query.executeQuery()) { while(rows.next()) add(rows.getString(2),"",rows.getString(1)); }
            }
            Path disease=p.resolveSibling("disease");
            if(Files.isDirectory(disease)) try(var q=connection.prepareStatement("SELECT id, name FROM read_parquet(?)")) {
                q.setString(1,disease.toAbsolutePath().toString().replace('\\','/')+"/*.parquet");
                try(var rows=q.executeQuery()) { while(rows.next()) names.put(rows.getString(1),rows.getString(2)); }
            }
        } catch(java.sql.SQLException e) { throw new IOException("cannot read Open Targets Parquet schema: "+e.getMessage(),e); }
    }
    private void obo(Path p) throws IOException {
        try(var r=reader(p)) {
            String id=null; boolean term=false;
            for(String line;(line=r.readLine())!=null;) {
                if(line.startsWith("data-version:")) metadata.setProperty("ontology_version",line.substring(13).trim());
                if(line.startsWith("[")) { term=line.equals("[Term]"); id=null; continue; }
                if(!term) continue;
                if(line.startsWith("id: ")) id=line.substring(4).trim();
                if(id==null) continue;
                if(line.startsWith("name: ")) names.put(id,line.substring(6).trim());
                if(line.startsWith("alt_id: ")) aliases.put(line.substring(8).trim(),id);
                if(line.equals("is_obsolete: true")) obsolete.add(id);
                if(line.startsWith("namespace: ")) aspects.put(id,switch(line.substring(11).trim()) {
                    case "biological_process" -> "BP"; case "molecular_function" -> "MF"; case "cellular_component" -> "CC"; default -> ""; });
                String parent=null;
                if(line.startsWith("is_a: ")) parent=line.substring(6).split(" ")[0];
                if(line.startsWith("relationship: part_of ")) parent=line.substring(22).split(" ")[0];
                if(parent!=null) parents.computeIfAbsent(id,ignored -> new TreeSet<>()).add(parent);
            }
        }
    }
    private void parentTable(Path p) throws IOException {
        try(var r=reader(p)) { for(String line;(line=r.readLine())!=null;) {
            String[] x=line.split("\t"); if(x.length!=2) throw new IOException("parent table needs parent and child");
            parents.computeIfAbsent(x[1],ignored -> new TreeSet<>()).add(x[0]);
        } }
    }
    private void propagate() throws IOException {
        Map<String,Set<String>> closure=new HashMap<>();
        for(String id:parents.keySet()) ancestors(id,closure,new HashSet<>());
        Map<String,Set<String>> original=new TreeMap<>(genes);
        for(var e:original.entrySet()) for(String parent:closure.getOrDefault(e.getKey(),Set.of()))
            if(!obsolete.contains(parent)) genes.computeIfAbsent(parent,ignored -> new TreeSet<>()).addAll(e.getValue());
    }
    private Set<String> ancestors(String id,Map<String,Set<String>> memo,Set<String> visiting) throws IOException {
        if(memo.containsKey(id)) return memo.get(id);
        if(!visiting.add(id)) throw new IOException("cycle in ontology: "+id);
        if(visiting.size()>1000) throw new IOException("ontology depth exceeds 1000");
        Set<String> result=new TreeSet<>();
        for(String parent:parents.getOrDefault(id,Set.of())) { result.add(parent); result.addAll(ancestors(parent,memo,visiting)); }
        visiting.remove(id); memo.put(id,result); return result;
    }
}
