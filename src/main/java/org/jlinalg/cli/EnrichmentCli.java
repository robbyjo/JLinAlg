/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;
import org.jlinalg.enrichment.EnrichmentAnalysis;
import org.jlinalg.enrichment.EnrichmentResult;

/** Input-derived universes, annotation joins, downloadable collections and enrichment reports. */
final class EnrichmentCli {
    private EnrichmentCli() { }
    private static final Set<String> OPTIONS=Set.of("--enrichment","--enrichment-db","--enrichment-method",
        "--input","--input-id","--input-cols","--gene-col","--gene-separator","--gene-id-type","--strip-gene-version",
        "--selection","--universe-selection","--missing-selection","--background","--background-id",
        "--annot","--annot-id","--annot-cols","--out","--output","--fdr-type","--min-set-size","--max-set-size",
        "--catalog-gene-source","--catalog-p-threshold","--trait-id-type","--species","--download","--db-release",
        "--source-url","--source-file","--db-format","--ontology","--ontology-url","--evidence","--propagate",
        "--monarch-evidence","--association-min-score","--kegg-academic","--max-download-bytes","--wallenius-max-updates");

    static int run(String[] args,PrintStream console,PrintStream error) {
        PrintStream log=null;
        try {
            if(Arrays.asList(args).contains("--help")) { console.println(help()); return 0; }
            Map<String,String> options=new LinkedHashMap<>();
            for(int i=0;i<args.length;i++) {
                String key=args[i];
                if(!OPTIONS.contains(key) || i+1==args.length) throw new IllegalArgumentException("unknown or incomplete enrichment option: "+key);
                key=key.equals("--output")?"--out":key.equals("--input-cols")?"--input-id":key;
                if(options.put(key,args[++i])!=null) throw new IllegalArgumentException("duplicate enrichment option: "+key);
            }
            String type=required(options,"--enrichment"); EnrichmentDatabase.provider(type);
            for(String key:List.of("--propagate","--strip-gene-version","--kegg-academic"))
                if(options.containsKey(key) && !Set.of("true","false").contains(options.get(key))) throw new IllegalArgumentException(key+" must be true or false");
            if(options.containsKey("--download")) {
                for(String key:List.of("--input","--selection","--universe-selection","--background","--out"))
                    if(options.containsKey(key)) throw new IllegalArgumentException("download and analysis modes cannot be combined: "+key);
                new EnrichmentDownloads(console,Long.parseLong(options.getOrDefault("--max-download-bytes","20000000000")))
                    .install(type,Path.of(options.get("--download")),options);
                return 0;
            }
            for(String key:List.of("--source-file","--source-url","--ontology-url","--db-release","--max-download-bytes","--kegg-academic"))
                if(options.containsKey(key)) throw new IllegalArgumentException(key+" requires --download");
            Path input=Path.of(required(options,"--input")).toAbsolutePath().normalize();
            Path database=Path.of(required(options,"--enrichment-db")).toAbsolutePath().normalize();
            Path output=Path.of(required(options,"--out")).toAbsolutePath().normalize();
            Path logPath=Path.of(output+".log"), auditPath=Path.of(output+".features.tsv"), metadataPath=Path.of(output+".manifest.properties");
            PipelinePaths.requireFreshOutputs(output,logPath,auditPath,metadataPath);
            Files.createDirectories(output.getParent());
            log=new PrintStream(Files.newOutputStream(logPath),true,java.nio.charset.StandardCharsets.UTF_8);
            String method=options.getOrDefault("--enrichment-method","ora").toLowerCase(Locale.ROOT);
            if(!Set.of("ora","gsameth").contains(method)) throw new IllegalArgumentException("enrichment method must be ora or gsameth");
            var correction=EnrichmentAnalysis.Fdr.valueOf(options.getOrDefault("--fdr-type","BH").toUpperCase(Locale.ROOT));
            info(console,log,"enrichment="+type+"; method="+method+"; term_fdr="+correction);
            EnrichmentDatabase db=EnrichmentDatabase.load(database,type,options);
            if(EnrichmentDatabase.provider(type).equals("GWASCatalog"))
                info(console,log,"catalog_gene_source="+db.metadata.getProperty("catalog_gene_source")+" (reported/mapped policy; deduplicated gene memberships)");
            for(String key:new TreeSet<>(db.metadata.stringPropertyNames())) info(console,log,"database."+key+"="+db.metadata.getProperty(key));
            Data selectedInput=read(input,options.get("--input-id"));
            if(options.containsKey("--selection") && options.containsKey("--background"))
                throw new IllegalArgumentException("use input-derived universe with --selection, or selected-only input with --background");
            if(!options.containsKey("--selection") && !options.containsKey("--background"))
                throw new IllegalArgumentException("provide --selection for a complete results table, or --background for a selected-only list");
            if(options.containsKey("--selection") && options.get("--input-id")==null)
                throw new IllegalArgumentException("table selections require --input-id (use no ID option only for headerless gene lists)");
            Data universe=options.containsKey("--background") ? read(Path.of(options.get("--background")),options.get("--background-id")) : selectedInput;
            EnrichmentExpression eligibility=new EnrichmentExpression(options.getOrDefault("--universe-selection","true"),universe.columns());
            EnrichmentExpression selection=new EnrichmentExpression(options.getOrDefault("--selection","true"),universe.columns());
            String missing=options.getOrDefault("--missing-selection","error");
            if(!Set.of("error","exclude").contains(missing)) throw new IllegalArgumentException("--missing-selection must be error or exclude");
            Map<String,Map<String,String>> eligible=new LinkedHashMap<>();
            int excluded=0, unevaluable=0;
            for(var row:universe.rows()) {
                Boolean value=eligibility.test(row,universe.rows().size());
                if(value==null) {
                    if(missing.equals("error")) throw new IllegalArgumentException("unevaluable universe expression for "+row.get(universe.idColumn()));
                    unevaluable++; continue;
                }
                if(!value) { excluded++; continue; }
                String id=row.get(universe.idColumn());
                if(EnrichmentExpression.missing(id)) throw new IllegalArgumentException("blank feature ID in universe");
                if(eligible.put(id,row)!=null) throw new IllegalArgumentException("duplicate eligible feature ID: "+id+"; select one contrast first");
            }
            if(eligible.isEmpty()) throw new IllegalArgumentException("eligible universe is empty");
            int provisionalCount=eligible.size();
            for(String id:new ArrayList<>(eligible.keySet())) {
                if(selection.test(eligible.get(id),provisionalCount)==null) {
                    if(missing.equals("error")) throw new IllegalArgumentException("unevaluable selection for "+id+"; use an eligibility filter or --missing-selection exclude");
                    eligible.remove(id); unevaluable++;
                }
            }
            int testCount=eligible.size(); if(testCount==0) throw new IllegalArgumentException("no evaluable eligible features");
            Set<String> selected=new LinkedHashSet<>();
            if(options.containsKey("--background")) {
                for(var row:selectedInput.rows()) {
                    String id=row.get(selectedInput.idColumn());
                    if(!eligible.containsKey(id)) throw new IllegalArgumentException("selected feature outside eligible background: "+id);
                    if(!selected.add(id)) throw new IllegalArgumentException("duplicate selected feature: "+id);
                }
            } else for(var e:eligible.entrySet()) if(Boolean.TRUE.equals(selection.test(e.getValue(),testCount))) selected.add(e.getKey());
            Map<String,Set<String>> mapping=new LinkedHashMap<>();
            Map<String,Map<String,Set<String>>> extra=new HashMap<>();
            String separator=options.getOrDefault("--gene-separator",";");
            if(separator.isEmpty()) throw new IllegalArgumentException("gene separator must be nonempty");
            List<String> annotationColumns=List.of();
            if(options.containsKey("--annot")) {
                Data annotation=read(Path.of(options.get("--annot")),options.getOrDefault("--annot-id",universe.idColumn()));
                String gene=required(options,"--gene-col");
                if(!annotation.columns().contains(gene)) throw new IllegalArgumentException("gene annotation column absent: "+gene);
                annotationColumns=options.containsKey("--annot-cols") ? List.of(options.get("--annot-cols").split(",")) : List.of();
                if(annotationColumns.equals(List.of("all"))) annotationColumns=annotation.columns().stream().sorted().toList();
                if(!annotation.columns().containsAll(annotationColumns)) throw new IllegalArgumentException("requested annotation columns absent");
                for(var row:annotation.rows()) {
                    String id=row.get(annotation.idColumn()); if(!eligible.containsKey(id)) continue;
                    mapping.computeIfAbsent(id,ignored->new TreeSet<>()).addAll(geneIds(row.get(gene),separator,options));
                    for(String column:annotationColumns) extra.computeIfAbsent(id,ignored->new HashMap<>()).computeIfAbsent(column,ignored->new TreeSet<>()).add(row.get(column));
                }
            } else {
                if(method.equals("gsameth")) throw new IllegalArgumentException("gsameth requires --annot and --gene-col with eligible CpG-to-gene mappings");
                if(options.containsKey("--annot-cols") || options.containsKey("--annot-id")) throw new IllegalArgumentException("annotation options require --annot");
                String gene=options.getOrDefault("--gene-col",universe.idColumn());
                if(!universe.columns().contains(gene)) throw new IllegalArgumentException("input gene column absent: "+gene);
                for(var e:eligible.entrySet()) mapping.put(e.getKey(),geneIds(e.getValue().get(gene),separator,options));
            }
            for(String id:eligible.keySet()) mapping.putIfAbsent(id,Set.of());
            Set<String> allGenes=new TreeSet<>(), selectedGenes=new TreeSet<>();
            int unmapped=0,selectedUnmapped=0;
            for(var e:mapping.entrySet()) {
                if(e.getValue().isEmpty()) { unmapped++; if(selected.contains(e.getKey())) selectedUnmapped++; }
                allGenes.addAll(e.getValue()); if(selected.contains(e.getKey())) selectedGenes.addAll(e.getValue());
            }
            if(allGenes.isEmpty()) throw new IllegalArgumentException("no eligible features map to genes");
            Set<String> collectionGenes=new TreeSet<>(); db.genes.values().forEach(collectionGenes::addAll);
            long represented=allGenes.stream().filter(collectionGenes::contains).count();
            if(represented==0) throw new IllegalArgumentException("no universe genes match database; check --gene-id-type and mapping");
            info(console,log,"universe_selection="+options.getOrDefault("--universe-selection","true"));
            info(console,log,"selection="+options.getOrDefault("--selection","explicit selected-only input"));
            info(console,log,"eligible_features="+testCount+"; bonferroni_default_m="+testCount+"; selected_features="+selected.size()
                +"; excluded_features="+excluded+"; unevaluable_excluded="+unevaluable);
            info(console,log,"universe_genes="+allGenes.size()+"; selected_genes="+selectedGenes.size()+"; unmapped_features="+unmapped
                +"; selected_unmapped="+selectedUnmapped+"; genes_without_collection_membership="+(allGenes.size()-represented));
            info(console,log,"background_policy=retain all mapped eligible genes, including genes without collection membership");
            if(method.equals("gsameth")) info(console,log,"gsameth=equivalent CpG coverage; fractional weights capped at one; rank tricube span=0.5; lexicographic ID ties; inclusive floor(weighted_overlap); no central fallback");
            int min=Integer.parseInt(options.getOrDefault("--min-set-size","10")),max=Integer.parseInt(options.getOrDefault("--max-set-size","500"));
            info(console,log,"term_size_range="+min+".."+max+" (after universe intersection)");
            List<EnrichmentResult> results=method.equals("ora") ? EnrichmentAnalysis.ora(allGenes,selectedGenes,db.sets(),min,max,correction)
                : EnrichmentAnalysis.gsameth(mapping,selected,db.sets(),min,max,correction,Long.parseLong(options.getOrDefault("--wallenius-max-updates","100000000")));
            if(results.isEmpty()) throw new IllegalArgumentException("no terms pass background-size filters");
            Properties manifest=new Properties(); manifest.putAll(db.metadata);
            options.forEach((k,v)->manifest.setProperty("command."+k,v));
            manifest.setProperty("enrichment_method",method); manifest.setProperty("term_fdr",correction.name());
            manifest.setProperty("min_set_size",Integer.toString(min)); manifest.setProperty("max_set_size",Integer.toString(max));
            manifest.setProperty("background_policy","all mapped eligible genes including genes without collection membership");
            manifest.setProperty("bonferroni_default_m",Integer.toString(testCount));
            manifest.setProperty("eligible_features",Integer.toString(testCount)); manifest.setProperty("selected_features",Integer.toString(selected.size()));
            manifest.setProperty("universe_genes",Integer.toString(allGenes.size())); manifest.setProperty("tested_terms",Integer.toString(results.size()));
            manifest.setProperty("input_sha256",EnrichmentDownloads.sha256(input));
            for(String key:List.of("--annot","--background")) if(options.containsKey(key)) manifest.setProperty(key.substring(2)+"_sha256",EnrichmentDownloads.sha256(Path.of(options.get(key))));
            writeResults(output,results,method);
            try(var writer=Files.newBufferedWriter(auditPath)) {
                List<String> header=new ArrayList<>(List.of("feature_id","selected","gene_ids"));
                for(String column:annotationColumns) header.add("annotation."+column);
                writeRow(writer,header,'\t');
                for(String id:eligible.keySet()) {
                    List<String> row=new ArrayList<>(List.of(id,Boolean.toString(selected.contains(id)),String.join(";",mapping.get(id))));
                    for(String column:annotationColumns) row.add(String.join(";",extra.getOrDefault(id,Map.of()).getOrDefault(column,Set.of())));
                    writeRow(writer,row,'\t');
                }
            }
            try(var writer=Files.newBufferedWriter(metadataPath)) { manifest.store(writer,"JLinAlg enrichment run"); }
            info(console,log,"complete: tested_terms="+results.size()+"; zero-hit terms retained; output="+output);
            return 0;
        } catch(InterruptedException e) {
            Thread.currentThread().interrupt(); error.println("jlinalg: enrichment download interrupted"); if(log!=null) log.println("FAILED: interrupted"); return 1;
        } catch(IOException | RuntimeException e) {
            error.println("jlinalg: "+e.getMessage()); if(log!=null) log.println("FAILED: "+e.getMessage()); return 2;
        } finally { if(log!=null) log.close(); }
    }
    private static Set<String> geneIds(String value,String separator,Map<String,String> options) {
        Set<String> result=new TreeSet<>();
        if(value==null) return result;
        for(String raw:value.split(Pattern.quote(separator))) {
            String gene=raw.trim(); if(EnrichmentDatabase.missing(gene)) continue;
            if("true".equals(options.get("--strip-gene-version"))) gene=gene.replaceFirst("\\.\\d+$","");
            result.add(gene);
        }
        return result;
    }
    private record Data(List<Map<String,String>> rows,Set<String> columns,String idColumn) { }
    private static Data read(Path path,String idColumn) throws IOException {
        List<Map<String,String>> rows=new ArrayList<>();
        if(idColumn==null) {
            try(var r=EnrichmentDatabase.reader(path)) { for(String line;(line=r.readLine())!=null;) {
                if(line.isBlank()) continue;
                if(line.contains("\t") || line.contains(",")) throw new IOException("table input requires an ID column option");
                rows.add(Map.of("id",line.trim()));
            } }
            return new Data(rows,Set.of("id"),"id");
        }
        EnrichmentDatabase.table(path,rows::add);
        if(rows.isEmpty()) throw new IOException("empty input table: "+path);
        if(!rows.get(0).containsKey(idColumn)) throw new IOException("ID column is absent: "+idColumn);
        return new Data(rows,Set.copyOf(rows.get(0).keySet()),idColumn);
    }
    private static void writeResults(Path output,List<EnrichmentResult> results,String method) throws IOException {
        Path temporary=Files.createTempFile(output.getParent(),".enrichment-", ".tmp");
        char delimiter=output.toString().toLowerCase(Locale.ROOT).endsWith(".csv")?',':'\t';
        try {
            try(var writer=Files.newBufferedWriter(temporary)) {
                writeRow(writer,List.of("term_id","term_name","method","universe_genes","selected_genes","term_genes","overlap_genes",
                    "weighted_overlap","expected_overlap_unweighted","fold_enrichment_unweighted","odds_ratio_unweighted","selection_odds",
                    "p_value","log_p_value","FDR","log_FDR","hit_genes"),delimiter);
                for(var r:results) writeRow(writer,List.of(r.id(),r.name(),method,Integer.toString(r.universeSize()),Integer.toString(r.selectedSize()),
                    Integer.toString(r.setSize()),Integer.toString(r.overlap()),number(r.weightedOverlap()),number(r.expectedOverlap()),number(r.foldEnrichment()),
                    number(r.oddsRatio()),number(r.selectionOdds()),number(r.pValue()),number(r.logPValue()),number(r.adjustedPValue()),number(r.logAdjustedPValue()),String.join(";",r.hitGenes())),delimiter);
            }
            Files.move(temporary,output);
        } finally { Files.deleteIfExists(temporary); }
    }
    static void writeRow(java.io.BufferedWriter writer,List<String> cells,char delimiter) throws IOException {
        for(int i=0;i<cells.size();i++) { if(i>0) writer.write(delimiter); writer.write('"'); writer.write(cells.get(i).replace("\"","\"\"")); writer.write('"'); } writer.newLine();
    }
    private static String number(double x) { return Double.isNaN(x)?"NA":Double.toString(x); }
    private static String required(Map<String,String> options,String key) { String value=options.get(key); if(value==null || value.isBlank()) throw new IllegalArgumentException("required option: "+key); return value; }
    private static void info(PrintStream console,PrintStream log,String message) { console.println(message); log.println(message); }
    static String help() { return """
        Enrichment (all expressions must be shell-quoted):
          --enrichment TYPE --enrichment-db PATH --input results.csv --input-id probe_id
          --selection "FDR < 0.05" [--universe-selection "status == 'OK' && n_studies >= 3"]
          [--annot annotation.csv --annot-id probe_id --gene-col gene_name --gene-separator ";"]
          --out enrichment.tsv [--enrichment-method ora|gsameth] [--fdr-type BH|BY|NONE]
        Selected-only lists: --input genes.txt --background tested-genes.txt (no --selection).
        Headered background tables require --background-id. --input-cols aliases --input-id.
        Defaults: ora; BH; term size 10..500; GWAS reported/mapped union; Catalog p<=5e-8;
          trait IDs=efo; missing-selection=error; propagation=true; gene-separator=";".
        --selection "p_value < bonferroni(0.05)" uses eligible feature count BEFORE gene mapping.
          bonferroni(0.05,850000) supplies the original correction-family size explicitly.
        Expression operators: && || ! < <= > >= == != + - * /; parentheses;
          abs(), is_missing(), is_finite(); backtick-quoted column names; quoted strings.
        --missing-selection exclude removes unevaluable rows before resolving default Bonferroni m.
        Other analysis options: --min-set-size --max-set-size --gene-id-type --strip-gene-version true
          --catalog-gene-source union|reported|mapped --catalog-p-threshold --trait-id-type efo|label
          --db-format gmt|term2gene|gaf|gwas|reactome|hpo|monarch|kegg --ontology FILE
          --evidence EXP,IDA,... --propagate true|false --association-min-score 0.5
          --wallenius-max-updates 100000000 --annot-cols c1,c2,... (or all).
        Downloads:
          --enrichment TYPE --download NEW_DIRECTORY [--species human] [--db-release latest]
          [--source-file FILE | --source-url HTTPS_URL] [--db-format gmt|term2gene|...]
          [--ontology FILE | --ontology-url HTTPS_URL] [--gene-id-type NAMESPACE]
          [--monarch-evidence causal|correlated] [--kegg-academic true]
          [--max-download-bytes 20000000000] (per source file).
        Reports: output, output.log, output.features.tsv and output.manifest.properties.
        Existing outputs/packages are never overwritten. See the enrichment vignette.
        """+EnrichmentDownloads.choices(); }
}
