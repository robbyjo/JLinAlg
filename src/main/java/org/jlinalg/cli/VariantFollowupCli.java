/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.zip.GZIPInputStream;
import static org.jlinalg.cli.FollowupSupport.*;

/** Local/remote annotation, external consequences and explicit evidence scoring. */
final class VariantFollowupCli {
    private VariantFollowupCli() { }
    private static final Set<String> COMMON=Set.of("out","input","genome-build","timeout","max-download-mb","token-env");
    static int run(String command,String[] args,PrintStream out,PrintStream err) {
        if(Arrays.asList(args).contains("--help")) { out.println(help(command)); return 0; }
        try {
            Set<String> allowed=new HashSet<>(COMMON);
            allowed.addAll(switch(command) {
                case "variant-db" -> Set.of("source-file","source-url","source-sha256","release","license","compression","format","chrom-column","pos-column","ref-column","alt-column");
                case "variant-annotate" -> Set.of("database","backend","endpoint");
                case "variant-consequence" -> Set.of("engine","annotations","format","executable","cache","cache-version","fasta","database","protocol","operation","endpoint","perl");
                case "variant-score" -> Set.of("components");
                default -> throw new IllegalArgumentException("Unknown variant operation");
            });
            Map<String,String> options=options(args,allowed);
            transaction(options,command,(dir,manifest)-> {
                switch(command) {
                    case "variant-db" -> install(options,dir,manifest);
                    case "variant-annotate" -> annotate(options,dir,manifest);
                    case "variant-consequence" -> consequence(options,dir,manifest);
                    case "variant-score" -> score(options,dir,manifest);
                    default -> throw new IllegalArgumentException(command);
                }
            });
            out.println(command+" completed: "+options.get("out")); return 0;
        } catch(Exception e) { err.println("jlinalg: "+e.getMessage()); return 2; }
    }
    static String help(String command) { return """
        Variant follow-up (outputs are NEW directories; no implicit downloads):
          variant-db --source-file annotations.tsv [or --source-url URL --source-sha256 SHA256]
            --genome-build GRCh38 --release VERSION --license LICENSE --out DATABASE
            [--compression none|gzip --format canonical|table|vcf]
            [--chrom-column CHROM --pos-column POS --ref-column REF --alt-column ALT]
          variant-annotate --input variants.tsv --genome-build GRCh38 --database DATABASE
            --out DIRECTORY [--backend local|post --endpoint URL]
          variant-consequence --input variants.tsv --genome-build GRCh38 --engine import
            --annotations FILE --format canonical|vep|annovar-vcf --out DIRECTORY
          variant-consequence --input variants.tsv --genome-build GRCh38 --engine vep
            --executable VEP --cache DIRECTORY --cache-version VERSION --fasta FASTA --out DIRECTORY
          variant-consequence --input variants.tsv --genome-build GRCh38 --engine annovar
            --executable table_annovar.pl --perl PERL --database DIRECTORY
            --protocol refGene --operation g --out DIRECTORY
          variant-consequence --input variants.tsv --genome-build GRCh38 --engine vep-rest
            --out DIRECTORY [--endpoint URL]
          variant-score --input evidence.tsv --components components.tsv --out DIRECTORY
        Canonical variants: genome_build,chrom,pos,ref,alt (biallelic VCF coordinates).
        Optional variant_id must equal genome_build:chrom:pos:ref:alt. Chromosome prefix chr is removed.
        Reference-normalize indels upstream; no implicit liftover, strand flip or reference repair.
        Local snapshots/POST responses: same five key columns plus annotation columns; repeated
        variants preserve transcript records. POST sends the canonical input TSV and expects TSV.
        Components: column,weight,minimum,maximum,direction (higher|lower). Missing component -> NA score.
        HTTP controls: --timeout SECONDS --max-download-mb N --token-env VARIABLE_NAME.
        All result directories contain a manifest with input/output SHA256 and settings.
        """; }
    static String build(Map<String,String> o) {
        String value=required(o,"genome-build");
        if(!Set.of("GRCh37","GRCh38").contains(value)) throw new IllegalArgumentException("genome-build must be GRCh37 or GRCh38"); return value;
    }
    static String key(DelimitedData data,String[] row,String expected) {
        String build=row[data.column("genome_build")];
        if(!build.equals(expected)) throw new IllegalArgumentException("Genome build mismatch: "+build+" vs "+expected);
        String chr=row[data.column("chrom")].replaceFirst("^chr","");
        if(!chr.matches("[A-Za-z0-9_.-]+"))throw new IllegalArgumentException("Invalid chromosome");
        long pos=Long.parseLong(row[data.column("pos")]);
        String ref=row[data.column("ref")].toUpperCase(Locale.ROOT),alt=row[data.column("alt")].toUpperCase(Locale.ROOT);
        if(pos<1 || !ref.matches("[ACGT]+") || !alt.matches("[ACGT]+") || ref.equals(alt))
            throw new IllegalArgumentException("Require positive position and distinct biallelic ACGT ref/alt (no symbolic variants)");
        String key=build+":"+chr+":"+pos+":"+ref+":"+alt;
        if(data.header().contains("variant_id") && !row[data.column("variant_id")].equals(key))
            throw new IllegalArgumentException("variant_id must match canonical alleles: "+key);
        return key;
    }
    private static Map<String,String[]> queries(DelimitedData data,String expected) {
        Map<String,String[]> queries=new LinkedHashMap<>();
        for(String[] row:data.rows()) { String id=key(data,row,expected); if(queries.put(id,row)!=null)throw new IllegalArgumentException("Duplicate input variant: "+id); }
        return queries;
    }
    private static void install(Map<String,String> o,Path dir,Map<String,Object> manifest) throws Exception {
        String build=build(o); required(o,"release"); required(o,"license");
        if(o.containsKey("source-file")==o.containsKey("source-url"))throw new IllegalArgumentException("Specify exactly one of source-file/source-url");
        Path source=dir.resolve("download");
        if(o.containsKey("source-file"))Files.copy(Path.of(o.get("source-file")),source);
        else { required(o,"source-sha256"); request(o.get("source-url"),null,"application/octet-stream",source,o); }
        String actual=hash(source); manifest.put("source_sha256",actual);
        if(o.containsKey("source-sha256") && !actual.equalsIgnoreCase(o.get("source-sha256")))throw new IOException("Source checksum mismatch");
        String compression=o.getOrDefault("compression","none"); Path table=dir.resolve("annotations.tsv");
        if(compression.equals("gzip")) {
            try(InputStream in=new GZIPInputStream(Files.newInputStream(source))) { Files.copy(in,table); }
            Files.delete(source);
        } else if(compression.equals("none"))Files.move(source,table); else throw new IllegalArgumentException("compression must be none or gzip");
        String format=o.getOrDefault("format","canonical");
        if(!format.equals("canonical")) {
            Path original=dir.resolve("source-table.txt");Files.move(table,original);
            try(BufferedReader in=Files.newBufferedReader(original);BufferedWriter out=Files.newBufferedWriter(table)) {
                if(format.equals("table")) {
                    String line=in.readLine();if(line==null)throw new IOException("Empty source table");
                    List<String> header=DelimitedData.parse(line,'\t',1,original);
                    int[] columns=new int[4];String[] keys={"chrom-column","pos-column","ref-column","alt-column"};
                    for(int i=0;i<4;i++){columns[i]=header.indexOf(required(o,keys[i]));if(columns[i]<0)throw new IOException("Source key column absent: "+o.get(keys[i]));}
                    List<String> names=new ArrayList<>(List.of("genome_build","chrom","pos","ref","alt"));names.addAll(header.stream().map(h->"source."+h).toList());
                    FollowupSupport.row(out,names.toArray(String[]::new));
                    while((line=in.readLine())!=null){if(line.isBlank())continue;List<String> values=DelimitedData.parse(line,'\t',0,original);if(values.size()!=header.size())throw new IOException("Source row width mismatch");List<String> result=new ArrayList<>(List.of(build));for(int c:columns)result.add(values.get(c));result.addAll(values);FollowupSupport.row(out,result.toArray(String[]::new));}
                } else if(format.equals("vcf")) {
                    FollowupSupport.row(out,new String[]{"genome_build","chrom","pos","ref","alt","source_id","info"});boolean header=false;
                    for(String line;(line=in.readLine())!=null;){if(line.startsWith("#CHROM"))header=true;if(line.startsWith("#")||line.isBlank())continue;String[] fields=line.split("\t",-1);if(!header||fields.length<8)throw new IOException("Invalid source VCF");FollowupSupport.row(out,new String[]{build,fields[0],fields[1],fields[3],fields[4],fields[2],fields[7]});}
                    if(!header)throw new IOException("VCF header absent");
                } else throw new IllegalArgumentException("format must be canonical, table or vcf");
            }
            Files.delete(original);
        }
        long count=stream(table,(data,row)->key(data,row,build));
        if(count==0)throw new IOException("Annotation database has no rows");
        manifest.put("rows",count); manifest.put("genome_build",build); manifest.put("release",o.get("release")); manifest.put("license",o.get("license"));
    }
    @FunctionalInterface interface RowAction { void accept(DelimitedData header,String[] row) throws Exception; }
    // Header access through a tiny parsed table; annotation database itself is streamed.
    private static long stream(Path path,RowAction action) throws Exception {
        Path stub=Files.createTempFile("jlinalg-header-",".tsv");
        try(BufferedReader in=Files.newBufferedReader(path)) {
            String first=in.readLine(); if(first==null)throw new IOException("Empty annotation table");
            List<String> names=DelimitedData.parse(first,'\t',1,path);
            Files.writeString(stub,first+"\n"+String.join("\t",Collections.nCopies(names.size(),"x"))+"\n");
            DelimitedData header=DelimitedData.read(stub); long count=0;
            for(String line;(line=in.readLine())!=null;) {
                if(line.isBlank())continue; List<String> row=DelimitedData.parse(line,'\t',count+2,path);
                if(row.size()!=names.size())throw new IOException("Annotation row width mismatch");
                action.accept(header,row.toArray(String[]::new)); count++;
            }
            return count;
        } finally {Files.deleteIfExists(stub);}
    }
    private static void annotate(Map<String,String> o,Path dir,Map<String,Object> manifest) throws Exception {
        String build=build(o); DelimitedData input=DelimitedData.read(Path.of(required(o,"input")));
        Map<String,String[]> queries=queries(input,build); Path annotations;
        switch(o.getOrDefault("backend","local")) {
            case "local" -> {
                Path db=Path.of(required(o,"database")); annotations=db.resolve("annotations.tsv");
                Map<String,Object> metadata;
                try(var reader=Files.newBufferedReader(db.resolve("manifest.yaml"))) {metadata=ProjectConfiguration.mapping(ProjectConfiguration.yaml().load(reader),"manifest");}
                if(!build.equals(metadata.get("genome_build")))throw new IllegalArgumentException("Database genome build mismatch");
                String expected=String.valueOf(ProjectConfiguration.mapping(metadata.get("output_sha256"),"hashes").get("annotations.tsv"));
                if(!hash(annotations).equals(expected))throw new IOException("Database checksum mismatch");
                manifest.put("database_manifest_sha256",hash(db.resolve("manifest.yaml"))); manifest.put("database_release",metadata.get("release"));
            }
            case "post" -> {
                annotations=dir.resolve("response.tsv");
                List<String[]> rows=queries.keySet().stream().map(id->id.split(":")) .toList();
                Path request=dir.resolve("request.tsv");table(request,List.of("genome_build","chrom","pos","ref","alt"),rows);
                request(required(o,"endpoint"),Files.readAllBytes(request),"text/tab-separated-values",annotations,o);
            }
            default -> throw new IllegalArgumentException("backend must be local or post");
        }
        join(input,queries,annotations,build,dir,manifest);
    }
    private static void join(DelimitedData input,Map<String,String[]> queries,Path annotations,String build,Path dir,Map<String,Object> manifest) throws Exception {
        List<String[]> rows=new ArrayList<>(); Set<String> found=new HashSet<>(); List<String> names=new ArrayList<>(List.of("variant_id","status"));
        try(var reader=Files.newBufferedReader(annotations)) {
            String line=reader.readLine();if(line==null)throw new IOException("Annotation header is absent");
            List<String> columns=DelimitedData.parse(line,'\t',1,annotations);
            if(!columns.containsAll(List.of("genome_build","chrom","pos","ref","alt")))throw new IOException("Annotation key columns are absent");
            names.addAll(columns.stream().map(s->"annotation."+s).toList());
        }
        stream(annotations,(header,row)-> {
            String key=key(header,row,build);
            if(queries.containsKey(key)) { String[] result=new String[row.length+2]; result[0]=key; result[1]="annotated";System.arraycopy(row,0,result,2,row.length);rows.add(result);found.add(key); }
        });
        for(String key:queries.keySet())if(!found.contains(key)) {String[] row=new String[names.size()];Arrays.fill(row,"NA");row[0]=key;row[1]="not_found";rows.add(row);}
        table(dir.resolve("annotations.tsv"),names,rows);
        manifest.put("matched_variants",found.size());manifest.put("unmatched_variants",queries.size()-found.size());
    }
    private static void consequence(Map<String,String> o,Path dir,Map<String,Object> manifest) throws Exception {
        String build=build(o); DelimitedData input=DelimitedData.read(Path.of(required(o,"input")));
        Map<String,String[]> queries=queries(input,build); String engine=required(o,"engine"); Path raw;
        if(engine.equals("import")) {
            raw=Path.of(required(o,"annotations"));
            if(o.getOrDefault("format","canonical").equals("canonical")) {join(input,queries,raw,build,dir,manifest);return;}
            if(o.get("format").equals("vep"))vep(raw,queries,dir);
            else if(o.get("format").equals("annovar-vcf"))annovar(raw,queries,dir);
            else throw new IllegalArgumentException("format must be canonical, vep or annovar-vcf");
            Files.copy(raw,dir.resolve("imported-annotations.txt")); return;
        }
        if(engine.equals("vep-rest")) {vepRest(o,dir,queries,manifest);return;}
        Path vcf=dir.resolve("input.vcf");
        List<String> records=new ArrayList<>(List.of("##fileformat=VCFv4.2","##reference="+build,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"));
        for(String id:queries.keySet()) {String[] parts=id.split(":");records.add(String.join("\t",parts[1],parts[2],id,parts[3],parts[4],".",".","."));}
        Files.write(vcf,records);
        List<String> command=new ArrayList<>();
        if(engine.equals("vep")) {
            command.add(executable(required(o,"executable")));
            command.addAll(List.of("--input_file",vcf.toString(),"--output_file",dir.resolve("vep.tsv").toString(),"--format","vcf","--tab","--offline","--cache","--dir_cache",Path.of(required(o,"cache")).toAbsolutePath().toString(),"--cache_version",required(o,"cache-version"),"--assembly",build,"--fasta",Path.of(required(o,"fasta")).toAbsolutePath().toString(),"--everything","--no_stats","--force_overwrite"));
        } else if(engine.equals("annovar")) {
            command.add(executable(o.getOrDefault("perl","perl")));command.add(Path.of(required(o,"executable")).toAbsolutePath().toString());
            command.addAll(List.of(vcf.toString(),Path.of(required(o,"database")).toAbsolutePath().toString(),"-buildver",build.equals("GRCh38")?"hg38":"hg19","-out",dir.resolve("annovar").toString(),"-protocol",required(o,"protocol"),"-operation",required(o,"operation"),"-nastring",".","-vcfinput"));
        } else throw new IllegalArgumentException("engine must be import, vep, annovar or vep-rest");
        manifest.put("external_command",command);process(command,dir,integer(o,"timeout",3600,1,86400));
        if(engine.equals("vep"))vep(dir.resolve("vep.tsv"),queries,dir);
        else annovar(dir.resolve("annovar."+(build.equals("GRCh38")?"hg38":"hg19")+"_multianno.vcf"),queries,dir);
    }
    private static final List<String> CONSEQUENCES=List.of("variant_id","gene","transcript","consequence","detail","status");
    private static void finish(List<String[]> rows,Map<String,String[]> queries,Path dir) throws IOException {
        Set<String> found=new HashSet<>();for(String[] row:rows) {if(!queries.containsKey(row[0]))throw new IOException("Unknown variant in annotation response: "+row[0]);found.add(row[0]);}
        for(String id:queries.keySet())if(!found.contains(id))rows.add(new String[]{id,"NA","NA","NA","NA","not_found"});
        table(dir.resolve("consequences.tsv"),CONSEQUENCES,rows);
    }
    private static void vep(Path file,Map<String,String[]> queries,Path dir) throws IOException {
        List<String> header=null;List<String[]> rows=new ArrayList<>();
        for(String line:Files.readAllLines(file)) {
            if(line.startsWith("##")||line.isBlank())continue;
            if(line.startsWith("#Uploaded_variation")) {header=Arrays.asList(line.substring(1).split("\t",-1));continue;}
            if(header==null)throw new IOException("VEP output requires #Uploaded_variation TSV header");
            String[] r=line.split("\t",-1);if(r.length!=header.size())throw new IOException("VEP row width mismatch");
            rows.add(new String[]{field(header,r,"Uploaded_variation"),field(header,r,"Gene"),field(header,r,"Feature"),field(header,r,"Consequence"),header.contains("Extra")?field(header,r,"Extra"):"NA","annotated"});
        }
        if(header==null)throw new IOException("Missing VEP header");finish(rows,queries,dir);
    }
    private static String field(List<String> header,String[] row,String name) throws IOException {int i=header.indexOf(name);if(i<0)throw new IOException("Missing annotation column: "+name);return row[i];}
    private static void annovar(Path file,Map<String,String[]> queries,Path dir) throws IOException {
        List<String[]> rows=new ArrayList<>();boolean header=false;
        for(String line:Files.readAllLines(file)) {
            if(line.startsWith("#CHROM"))header=true;
            if(line.startsWith("#")||line.isBlank())continue;
            String[] r=line.split("\t",-1);if(r.length<8)throw new IOException("Invalid ANNOVAR VCF");
            String key=queries.keySet().stream().filter(k->k.equals(r[2])).findFirst().orElseThrow(()->new IOException("ANNOVAR must preserve canonical VCF IDs"));
            String[] parts=key.split(":");
            if(!r[0].replaceFirst("^chr","").equals(parts[1])||!r[1].equals(parts[2])||!r[3].equals(parts[3])||!r[4].equals(parts[4]))throw new IOException("ANNOVAR returned changed alleles");
            Map<String,String> info=new LinkedHashMap<>();for(String item:r[7].split(";")) {String[] kv=item.split("=",2);if(kv.length==2)info.put(kv[0],kv[1]);}
            boolean found=false;
            for(String name:info.keySet())if(name.startsWith("Gene.")) {String protocol=name.substring(5);rows.add(new String[]{key,info.get(name),"NA",info.getOrDefault("ExonicFunc."+protocol,info.getOrDefault("Func."+protocol,"NA")),r[7],"annotated"});found=true;}
            if(!found)rows.add(new String[]{key,"NA","NA","NA",r[7],"annotated"});
        }
        if(!header)throw new IOException("Missing VCF header");finish(rows,queries,dir);
    }
    private static void vepRest(Map<String,String> o,Path dir,Map<String,String[]> queries,Map<String,Object> manifest) throws Exception {
        String endpoint=o.getOrDefault("endpoint",build(o).equals("GRCh38")?"https://rest.ensembl.org/vep/human/region":"https://grch37.rest.ensembl.org/vep/human/region");
        List<String> ids=new ArrayList<>(queries.keySet());List<String[]> rows=new ArrayList<>();
        for(int start=0;start<ids.size();start+=200) {
            List<String> lines=new ArrayList<>();Map<String,String> sent=new HashMap<>();
            for(String id:ids.subList(start,Math.min(start+200,ids.size()))) {String[] p=id.split(":");String line=String.join(" ",p[1],p[2],id,p[3],p[4],".",".",".");lines.add("\""+line+"\"");sent.put(line,id);}
            Path response=dir.resolve("response-"+start+".json");
            request(endpoint,("{\"variants\":["+String.join(",",lines)+"]}").getBytes(StandardCharsets.UTF_8),"application/json",response,o);
            Object parsed;try(var reader=Files.newBufferedReader(response)){parsed=ProjectConfiguration.yaml().load(reader);}
            if(!(parsed instanceof List<?> records))throw new IOException("Expected VEP JSON array");
            Set<String> returned=new HashSet<>();
            for(Object value:records) {
                Map<String,Object> item=ProjectConfiguration.mapping(value,"VEP result");String id=sent.get(String.valueOf(item.get("input")));
                if(id==null || !returned.add(id))throw new IOException("Unmatched/duplicate VEP response input");
                if(!build(o).equals(String.valueOf(item.get("assembly_name"))))throw new IOException("VEP response genome build mismatch");
                Object consequences=item.get("transcript_consequences");
                if(consequences instanceof List<?> transcripts && !transcripts.isEmpty())for(Object t:transcripts) {
                    Map<String,Object> tx=ProjectConfiguration.mapping(t,"transcript");
                    rows.add(new String[]{id,String.valueOf(tx.getOrDefault("gene_id","NA")),String.valueOf(tx.getOrDefault("transcript_id","NA")),String.valueOf(tx.getOrDefault("consequence_terms","NA")),"response-"+start+".json","annotated"});
                } else rows.add(new String[]{id,"NA","NA",String.valueOf(item.getOrDefault("most_severe_consequence","NA")),"response-"+start+".json","annotated"});
            }
            if(returned.size()!=sent.size())throw new IOException("Incomplete VEP batch response");
        }
        manifest.put("endpoint",endpoint);finish(rows,queries,dir);
    }
    private static void score(Map<String,String> o,Path dir,Map<String,Object> manifest) throws Exception {
        DelimitedData data=DelimitedData.read(Path.of(required(o,"input"))),components=DelimitedData.read(Path.of(required(o,"components")));
        List<String> header=new ArrayList<>(data.header());Set<String> used=new HashSet<>();
        int k=components.rows().size();int[] index=new int[k];double[] weights=new double[k],low=new double[k],high=new double[k];boolean[] reverse=new boolean[k];double sum=0;
        for(int c=0;c<k;c++) {
            String[] row=components.rows().get(c);String name=row[components.column("column")];
            if(!used.add(name))throw new IllegalArgumentException("Duplicate score component");index[c]=data.column(name);
            weights[c]=finite(row[components.column("weight")],"weight");low[c]=finite(row[components.column("minimum")],"minimum");high[c]=finite(row[components.column("maximum")],"maximum");
            String direction=row[components.column("direction")];if(!Set.of("higher","lower").contains(direction))throw new IllegalArgumentException("direction must be higher or lower");reverse[c]=direction.equals("lower");
            if(weights[c]<=0 || high[c]<=low[c] || !Double.isFinite(high[c]-low[c]))throw new IllegalArgumentException("Require positive weights and a finite positive maximum - minimum");sum+=weights[c];header.add("contribution."+name);
        }
        if(!Double.isFinite(sum))throw new IllegalArgumentException("Weight sum overflow");
        header.addAll(List.of("priority_score","weight_coverage","score_status"));
        if(new HashSet<>(header).size()!=header.size())throw new IllegalArgumentException("Input collides with score output columns");
        List<String[]> rows=new ArrayList<>();
        for(String[] row:data.rows()) {
            String[] result=Arrays.copyOf(row,header.size());double score=0,covered=0;boolean complete=true;
            for(int c=0;c<k;c++) {
                String raw=row[index[c]];
                if(Set.of("NA",".","").contains(raw)) {result[row.length+c]="NA";complete=false;continue;}
                double value=finite(raw,"component");if(value<low[c]||value>high[c])throw new IllegalArgumentException("Component outside declared bounds: "+raw);
                double scaled=(value-low[c])/(high[c]-low[c]);if(reverse[c])scaled=1-scaled;
                double contribution=weights[c]/sum*scaled;score+=contribution;covered+=weights[c];result[row.length+c]=Double.toString(contribution);
            }
            result[header.size()-3]=complete?Double.toString(score):"NA";result[header.size()-2]=Double.toString(covered/sum);result[header.size()-1]=complete?"complete":"missing_evidence";rows.add(result);
        }
        table(dir.resolve("scores.tsv"),header,rows);manifest.put("interpretation","Weighted prioritization index, not a causal or pathogenicity probability; rows are not collapsed across transcripts/genes.");
    }
}
