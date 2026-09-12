/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.IOException;
import java.nio.file.*;
import java.util.*;

/** Strict labeled input/output helpers for the summary xWAS commands. */
final class XwasFiles {
    private XwasFiles() { }
    static Map<String,String> options(String[] args,String... allowed) {
        Set<String> keys=Set.of(allowed);Map<String,String> result=new LinkedHashMap<>();
        for(int i=0;i<args.length;i++) {
            String key=args[i];
            if(!keys.contains(key)||i+1==args.length||args[i+1].startsWith("--")||result.put(key,args[++i])!=null)
                throw new IllegalArgumentException("unknown, duplicate, or incomplete option: "+key);
        }
        return result;
    }
    static String required(Map<String,String> o,String k) {
        String v=o.get(k);if(v==null||v.isBlank())throw new IllegalArgumentException("required option: "+k);return v;
    }
    static Path path(Map<String,String> o,String k){return Path.of(required(o,k)).toAbsolutePath().normalize();}
    static double number(String x) {
        double v=Double.parseDouble(x);if(!Double.isFinite(v))throw new IllegalArgumentException("finite number required: "+x);return v;
    }
    static double number(DelimitedData t,String[] row,String column){return number(row[t.column(column)]);}
    static String id(String value) {
        if(value.isBlank()||!value.equals(value.trim())||value.contains("\t")||value.contains("\n")||value.contains("\r")||value.contains("\""))
            throw new IllegalArgumentException("identifiers must be nonblank, single-line, unpadded and without quotes");return value;
    }
    static Map<String,String[]> indexed(DelimitedData t,String key) {
        int c=t.column(key);Map<String,String[]> map=new LinkedHashMap<>();
        for(String[] row:t.rows())if(map.put(id(row[c]),row)!=null)throw new IllegalArgumentException("duplicate "+key+": "+row[c]);
        return map;
    }
    static List<String> names(String x) {
        List<String> names=List.of(x.split(",",-1));Set<String> seen=new HashSet<>();
        for(String name:names)if(!seen.add(id(name)))throw new IllegalArgumentException("duplicate name: "+name);
        return names;
    }
    static double[] matrix(Path path,List<String> labels) throws IOException {
        DelimitedData t=DelimitedData.read(path);
        if(t.header().size()!=labels.size()+1)throw new IllegalArgumentException("matrix column count differs from labels: "+path);
        Map<String,String[]> rows=indexed(t,t.header().get(0));
        if(!rows.keySet().equals(new HashSet<>(labels)))throw new IllegalArgumentException("matrix row labels differ: "+path);
        int n=labels.size();double[] a=new double[n*n];
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)a[i*n+j]=number(rows.get(labels.get(i))[t.column(labels.get(j))]);
        return a;
    }
    static String matrix(List<String> labels,double[] a) {
        String key="row";while(labels.contains(key))key="_"+key;
        StringBuilder text=new StringBuilder(key+"\t"+String.join("\t",labels)+"\n");
        int n=labels.size();
        for(int i=0;i<n;i++){text.append(labels.get(i));for(int j=0;j<n;j++)text.append('\t').append(a[i*n+j]);text.append('\n');}
        return text.toString();
    }
    /** Preflight every destination, stage all contents, then publish without replacement. */
    static void publish(Map<Path,String> files) throws IOException {
        for(Path path:files.keySet())if(!path.toString().endsWith(".tsv"))
            throw new IllegalArgumentException("xWAS output files must end in .tsv");
        PipelinePaths.requireFreshOutputs(files.keySet().toArray(Path[]::new));
        Map<Path,Path> staged=new LinkedHashMap<>();List<Path> published=new ArrayList<>();
        try {
            for(var entry:files.entrySet()) {
                Path target=entry.getKey().toAbsolutePath();
                Files.createDirectories(target.getParent());
                Path tmp=Files.createTempFile(target.getParent(),".xwas-",".tmp");staged.put(target,tmp);
                Files.writeString(tmp,entry.getValue());
            }
            // Publish the primary result last; roll back only newly created outputs on I/O failure.
            List<Path> targets=new ArrayList<>(staged.keySet());
            Collections.rotate(targets,-1);
            for(Path target:targets){Files.move(staged.get(target),target);published.add(target);}
        } catch(IOException failure) {
            for(Path path:published)try{Files.deleteIfExists(path);}catch(IOException cleanup){failure.addSuppressed(cleanup);}
            throw failure;
        } finally {for(Path p:staged.values())Files.deleteIfExists(p);}
    }
    static List<String> pairs(List<String> traits) {
        List<String> result=new ArrayList<>();
        for(int i=0;i<traits.size();i++)for(int j=i;j<traits.size();j++)result.add(traits.get(i)+":"+traits.get(j));
        if(new HashSet<>(result).size()!=result.size())throw new IllegalArgumentException("ambiguous trait pair labels");
        return result;
    }
}
