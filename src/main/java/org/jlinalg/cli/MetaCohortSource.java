/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

/** External sort and streaming reader for one cohort; bounded merge fan-in. */
final class MetaCohortSource implements AutoCloseable {
    record Row(String id, double effect, double se) { }
    private final BufferedReader reader;
    private Row current;

    MetaCohortSource(Path input, String idColumn, String effectColumn, String seColumn,
            Path scratch, int chunkRows) throws IOException {
        List<Path> runs = new ArrayList<>();
        try (BufferedReader in = open(input)) {
            char separator = delimiter(input);
            String first = in.readLine();
            if (first == null) throw new IOException("empty cohort file: " + input);
            List<String> header = DelimitedData.parse(first.replaceFirst("^\uFEFF", ""), separator, 1, input);
            header.replaceAll(String::trim);
            if (new HashSet<>(header).size() != header.size() || header.contains(""))
                throw new IOException("cohort columns must be unique and nonblank: " + input);
            int id = column(header,idColumn), effect = column(header,effectColumn), se = column(header,seColumn);
            List<Row> rows = new ArrayList<>(chunkRows);
            long lineNumber = 1;
            for (String line; (line=in.readLine())!=null;) {
                lineNumber++;
                if (line.isBlank()) continue;
                List<String> fields = DelimitedData.parse(line,separator,lineNumber,input);
                if (fields.size()!=header.size()) throw new IOException("wrong field count at " + input + ":" + lineNumber);
                String key=fields.get(id).trim();
                if(key.isEmpty() || key.chars().anyMatch(Character::isISOControl))
                    throw new IOException("invalid feature ID at " + input + ":" + lineNumber);
                double beta=number(fields.get(effect),input,lineNumber), error=number(fields.get(se),input,lineNumber);
                if (!Double.isNaN(error) && (!(error>0) || !Double.isFinite(error*error) || error*error==0))
                    throw new IOException("SE must have a finite positive variance at " + input + ":" + lineNumber);
                if (Double.isNaN(beta) || Double.isNaN(error)) { beta=Double.NaN; error=Double.NaN; }
                rows.add(new Row(key,beta,error));
                if(rows.size()==chunkRows) { runs.add(spill(rows,scratch)); rows.clear(); }
            }
            if(!rows.isEmpty() || runs.isEmpty()) runs.add(spill(rows,scratch));
        }
        while(runs.size()>1) {
            List<Path> next=new ArrayList<>();
            for(int i=0;i<runs.size();i+=32) {
                List<Path> group=runs.subList(i,Math.min(i+32,runs.size()));
                Path merged=Files.createTempFile(scratch,"merge-",".tsv");
                merge(group,merged);
                next.add(merged);
                for(Path path:group) Files.delete(path);
            }
            runs=next;
        }
        reader=Files.newBufferedReader(runs.get(0));
        advance();
    }

    Row current() { return current; }
    void advance() throws IOException { current=read(reader); }
    @Override public void close() throws IOException { reader.close(); }

    static BufferedReader open(Path input) throws IOException {
        InputStream stream=Files.newInputStream(input);
        try {
            if(input.toString().toLowerCase(Locale.ROOT).endsWith(".gz")) stream=new GZIPInputStream(stream);
            return new BufferedReader(new InputStreamReader(stream,java.nio.charset.StandardCharsets.UTF_8));
        } catch(IOException failure) { stream.close(); throw failure; }
    }
    static char delimiter(Path input) {
        String name=input.toString().toLowerCase(Locale.ROOT);
        return name.endsWith(".csv") || name.endsWith(".csv.gz") ? ',' : '\t';
    }
    static int column(List<String> header,String name) {
        int index=header.indexOf(name);
        if(index<0) throw new IllegalArgumentException("column is absent: " + name);
        return index;
    }
    static double number(String text,Path path,long line) throws IOException {
        String value=text.trim();
        if(value.isEmpty() || value.equals(".") || value.equalsIgnoreCase("NA") || value.equalsIgnoreCase("NaN")) return Double.NaN;
        try {
            double number=Double.parseDouble(value);
            if(!Double.isFinite(number)) throw new NumberFormatException();
            return number;
        } catch(NumberFormatException failure) { throw new IOException("invalid number '"+value+"' at "+path+":"+line); }
    }
    private static Path spill(List<Row> rows,Path scratch) throws IOException {
        rows.sort(Comparator.comparing(Row::id));
        Path path=Files.createTempFile(scratch,"cohort-",".tsv");
        try(BufferedWriter out=Files.newBufferedWriter(path)) {
            String previous=null;
            for(Row row:rows) {
                duplicate(previous,row.id()); previous=row.id(); write(out,row);
            }
        }
        return path;
    }
    private static void merge(List<Path> paths,Path output) throws IOException {
        List<BufferedReader> readers=new ArrayList<>();
        record Head(int source,Row row) { }
        PriorityQueue<Head> queue=new PriorityQueue<>(Comparator.comparing(h->h.row().id()));
        try(BufferedWriter out=Files.newBufferedWriter(output)) {
            for(Path path:paths) {
                BufferedReader in=Files.newBufferedReader(path); readers.add(in);
                Row row=read(in); if(row!=null) queue.add(new Head(readers.size()-1,row));
            }
            String previous=null;
            while(!queue.isEmpty()) {
                Head head=queue.remove(); duplicate(previous,head.row().id()); previous=head.row().id();
                write(out,head.row()); Row next=read(readers.get(head.source()));
                if(next!=null) queue.add(new Head(head.source(),next));
            }
        } finally { for(BufferedReader in:readers) in.close(); }
    }
    private static void duplicate(String previous,String id) throws IOException {
        if(id.equals(previous)) throw new IOException("duplicate feature ID within cohort: " + id);
    }
    private static Row read(BufferedReader in) throws IOException {
        String line=in.readLine(); if(line==null) return null;
        String[] fields=line.split("\t",-1);
        return new Row(fields[0],Double.parseDouble(fields[1]),Double.parseDouble(fields[2]));
    }
    private static void write(BufferedWriter out,Row row) throws IOException {
        out.write(row.id()+"\t"+row.effect()+"\t"+row.se()+"\n");
    }
}
