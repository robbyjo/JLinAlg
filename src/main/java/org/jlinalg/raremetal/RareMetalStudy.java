/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.raremetal;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.zip.GZIPInputStream;
import htsjdk.tribble.readers.TabixReader;

/** Biallelic RAREMETALWORKER/rvtests summary reader. Indexed files support
 * bounded regional access; unindexed files are rescanned without retaining a
 * genome-wide matrix. Covariances are restored to the score-information scale.
 * Regional reads are synchronized; independent numerical work can run in parallel.
 */
public final class RareMetalStudy implements AutoCloseable {
    private final Path scores, covariance;
    private final Map<String,Integer> header=new HashMap<>();
    private final TabixReader scoreIndex,covIndex;
    private final Map<String,String> scoreChromosomes,covChromosomes;
    private int samples;
    private String genomeBuild;
    private boolean scoreScale;
    private double minimumCallRate,minimumHwe;
    private long cacheLimit,cacheSize;
    private final LinkedHashMap<String,CacheEntry> cache=new LinkedHashMap<>(16,.75f,true);
    private record CacheEntry(Object value,long bytes) { }
    /** Bound estimated retained payload bytes; zero disables caching. */
    public synchronized RareMetalStudy cacheBytes(long bytes) {
        if(bytes<0)throw new IllegalArgumentException("negative cache budget");
        cacheLimit=bytes;cache.clear();cacheSize=0;return this;
    }
    private void remember(String key,Object value,long bytes) {
        if(bytes>cacheLimit||cacheLimit==0)return;
        while(!cache.isEmpty()&&cacheSize+bytes>cacheLimit) {
            var first=cache.entrySet().iterator();var entry=first.next();cacheSize-=entry.getValue().bytes;first.remove();
        }
        cache.put(key,new CacheEntry(value,bytes));cacheSize+=bytes;
    }

    public RareMetalStudy(Path scores,Path covariance,double minimumCallRate,double minimumHwe) throws IOException {
        this.scores=scores; this.covariance=covariance;
        this.minimumCallRate=minimumCallRate; this.minimumHwe=minimumHwe;
        if(minimumCallRate<0 || minimumCallRate>1 || minimumHwe<0 || minimumHwe>1
                || !Double.isFinite(minimumCallRate) || !Double.isFinite(minimumHwe))
            throw new IllegalArgumentException("QC cutoffs must be in [0,1]");
        try(BufferedReader r=open(scores)) {
            String line;
            while((line=r.readLine())!=null) {
                if(line.startsWith("##AnalyzedSamples=")) samples=Integer.parseInt(line.substring(line.indexOf('=')+1).trim());
                if(line.startsWith("##GenomeBuild=")) genomeBuild=line.substring(line.indexOf('=')+1).trim();
                if(line.equals("##CovarianceScale=score")) scoreScale=true;
                String clean=line.replaceFirst("^#+","");
                if(clean.startsWith("CHROM")) {
                    String[] names=fields(clean);
                    for(int i=0;i<names.length;i++) header.put(names[i],i);
                    break;
                }
            }
        }
        for(String required:List.of("CHROM","POS","REF","ALT","N_INFORMATIVE","U_STAT","SQRT_V_STAT"))
            if(!header.containsKey(required)) throw new IOException("missing score column "+required+" in "+scores);
        if(samples<=0) throw new IOException("positive ##AnalyzedSamples is required: "+scores);
        TabixReader s=null,c=null;
        try {
            if(Files.exists(Path.of(scores+".tbi"))) s=new TabixReader(scores.toString());
            if(covariance!=null && Files.exists(Path.of(covariance+".tbi"))) c=new TabixReader(covariance.toString());
            scoreIndex=s; covIndex=c;
            scoreChromosomes=chromosomes(s); covChromosomes=chromosomes(c);
        } catch(IOException|RuntimeException e) { if(s!=null)s.close(); if(c!=null)c.close(); throw e; }
    }
    public int samples(){return samples;}
    /** Declared build, or null for historical files without a build header. */
    public String genomeBuild(){return genomeBuild;}
    public boolean indexed(){return scoreIndex!=null && (covariance==null || covIndex!=null);}

    /** Complete sequential scan; caller closes the cursor. Input must be position sorted. */
    public Cursor cursor() throws IOException {return new Cursor(open(scores));}

    /** At most one biallelic record per position; ambiguous multiallelic rows fail. */
    public synchronized Map<Long,Score> region(String chromosome,long start,long end) throws IOException {
        String key="scores:"+chromosome+":"+start+":"+end;
        CacheEntry cached=cache.get(key);
        if(cached!=null) {
            @SuppressWarnings("unchecked") Map<Long,Score> value=(Map<Long,Score>)cached.value;
            return value;
        }
        Map<Long,Score> result=new HashMap<>();
        lines(scores,scoreIndex,scoreChromosomes,chromosome,start,end,line->{
            Score value=parse(line);
            if(result.put(value.position(),value)!=null) throw new IOException("duplicate/multiallelic score position: "+value.id());
        });
        Map<Long,Score> immutable=Map.copyOf(result);
        long bytes=128+2L*key.length();
        for(Score value:result.values())bytes+=256+2L*(value.chromosome().length()+value.reference().length()+value.alternate().length());
        remember(key,immutable,bytes);
        return immutable;
    }

    /** Covariance block in the order of supplied, informative score records. */
    public synchronized double[] covariance(List<Score> records) throws IOException {
        String key="cov:"+records.stream().map(s->s.id()+":"+s.variance()).collect(java.util.stream.Collectors.joining(";"));
        CacheEntry cached=cache.get(key);
        if(cached!=null)return ((double[])cached.value).clone();
        int n=records.size(); double[] matrix=new double[Math.multiplyExact(n,n)];
        if(n==0)return matrix;
        for(int i=0;i<n;i++) matrix[i*n+i]=records.get(i).variance();
        if(n==1)return matrix;
        if(covariance==null) throw new IOException("group tests require a covariance file: "+scores);
        Map<Long,Integer> wanted=new HashMap<>();
        long low=Long.MAX_VALUE,high=0; String chromosome=records.get(0).chromosome();
        for(int i=0;i<n;i++) {
            Score s=records.get(i);
            if(!s.chromosome().equals(chromosome) || wanted.put(s.position(),i)!=null)
                throw new IOException("covariance block must contain unique positions on one chromosome");
            low=Math.min(low,s.position()); high=Math.max(high,s.position());
        }
        boolean[] found=new boolean[matrix.length]; boolean[] rows=new boolean[n];
        lines(covariance,covIndex,covChromosomes,chromosome,low,high,line->{
            String[] f=fields(line); long position=Long.parseLong(f[1]); Integer i=wanted.get(position);
            if(i==null)return;
            if(rows[i])throw new IOException("duplicate covariance row at "+position);
            rows[i]=true;
            // RMW: CHROM CURRENT_POS MARKERS_IN_WINDOW COV_MATRICES.
            // rvtests: CHROM START_POS END_POS NUM_MARKER MARKER_POS COV.
            int markers=f.length==4?2:f.length==6?4:-1;
            if(markers<0)throw new IOException("unsupported covariance layout at "+position);
            String[] positions=f[markers].split(","), values=f[markers+1].split(",");
            if(positions.length!=values.length)throw new IOException("covariance position/value count mismatch");
            for(int z=0;z<positions.length;z++) {
                Integer j=wanted.get(Long.parseLong(positions[z])); if(j==null)continue;
                double value=Double.parseDouble(values[z])*(scoreScale?1:samples);
                if(!Double.isFinite(value))throw new IOException("nonfinite covariance");
                if(i.equals(j)) {
                    double expected=records.get(i).variance();
                    if(Math.abs(value-expected)>2e-4*Math.max(expected,Math.abs(value)))
                        throw new IOException("covariance diagonal/score variance disagree at "+position);
                } else {
                    if(found[i*n+j] && Math.abs(matrix[i*n+j]-value)>2e-4*Math.max(Math.abs(value),Math.abs(matrix[i*n+j])))
                        throw new IOException("inconsistent duplicate covariance pair");
                    matrix[i*n+j]=matrix[j*n+i]=value;
                }
                found[i*n+j]=found[j*n+i]=true;
            }
        });
        for(int i=0;i<n;i++) for(int j=i;j<n;j++) if(!found[i*n+j])
            throw new IOException("missing covariance for "+records.get(i).id()+" / "+records.get(j).id()+"; covariance window may be too short");
        long bytes=8L*matrix.length+128+2L*key.length();
        if(cacheLimit>0&&bytes<=cacheLimit)remember(key,matrix.clone(),bytes);
        return matrix;
    }

    private Score parse(String line) throws IOException {
        try {
            String[] f=fields(line);
            String ref=text(f,"REF"),alt=text(f,"ALT");
            if(alt.equals(".")||alt.equals("0"))alt="";
            if(!allele(ref) || !alt.isEmpty() && (!allele(alt) || ref.equals(alt)))
                throw new IllegalArgumentException("only normalized biallelic A/C/G/T alleles are supported");
            double n=number(f,"N_INFORMATIVE"),af=header.containsKey("ALL_AF")?number(f,"ALL_AF"):number(f,"AF");
            // Historical RMW writes the sole observed allele in REF and leaves
            // ALT blank at monomorphic sites, sometimes with ALL_AF=1.
            double homRef=number(f,"N_REF"),het=number(f,"N_HET"),homAlt=number(f,"N_ALT");
            boolean soleHomozygote=het==0&&(homRef>0&&homAlt==0||homAlt>0&&homRef==0)
                && !Double.isFinite(number(f,"U_STAT"))&&!Double.isFinite(number(f,"SQRT_V_STAT"));
            boolean invalidAlleles=alt.isEmpty()&&af!=0&&af!=1&&!soleHomozygote;
            if(alt.isEmpty()&&!invalidAlleles)af=0;
            if(!Double.isFinite(n)||n<0||n>samples)throw new IllegalArgumentException("invalid informative sample count");
            if(!Double.isNaN(af)&&(!Double.isFinite(af)||af<0||af>1))throw new IllegalArgumentException("invalid allele frequency");
            double u=number(f,"U_STAT"),root=number(f,"SQRT_V_STAT");
            if(Double.isInfinite(u)||Double.isInfinite(root)||root<0)throw new IllegalArgumentException("invalid score or variance");
            boolean qc=(minimumCallRate==0 || number(f,"CALL_RATE")>=minimumCallRate)
                && (minimumHwe==0 || number(f,"HWE_PVALUE")>=minimumHwe);
            if(invalidAlleles)af=Double.NaN;
            String status=invalidAlleles?"invalid_alleles":!qc?"qc_excluded":n==0||Double.isNaN(af)?"missing":
                af==0||af==1?"monomorphic":!Double.isFinite(u)||!Double.isFinite(root)||root==0?"missing":"ok";
            long position=Long.parseLong(text(f,"POS"));
            if(position<1||position>Integer.MAX_VALUE)throw new IllegalArgumentException("position outside supported genomic range");
            if(Double.isFinite(root)&&!Double.isFinite(root*root))throw new IllegalArgumentException("score variance overflow");
            return new Score(chromosome(text(f,"CHROM")),position,ref,alt,n,af,
                status.equals("ok")?u:Double.NaN,status.equals("ok")?root*root:Double.NaN,status,number(f,"CALL_RATE"),number(f,"HWE_PVALUE"));
        } catch(RuntimeException e){throw new IOException("invalid score row in "+scores+": "+e.getMessage(),e);}
    }
    private String text(String[] f,String name){return f[header.get(name)];}
    private static boolean allele(String s){if(s.isEmpty())return false;for(int i=0;i<s.length();i++)if("ACGT".indexOf(s.charAt(i))<0)return false;return true;}
    private double number(String[] f,String name){Integer i=header.get(name);if(i==null)return Double.NaN;String s=f[i];return s.equals("NA")||s.equals(".")||s.equals("NaN")?Double.NaN:Double.parseDouble(s);}
    private static String[] fields(String line){return line.indexOf('\t')>=0?line.split("\t",-1):line.trim().split("\\s+");}
    public static String chromosome(String value){return value.startsWith("chr")?value.substring(3):value;}
    private static Map<String,String> chromosomes(TabixReader reader) throws IOException {
        Map<String,String> result=new HashMap<>();
        if(reader!=null)for(String c:reader.getChromosomes()) if(result.put(chromosome(c),c)!=null)
            throw new IOException("ambiguous chromosome aliases in index");
        return result;
    }
    @FunctionalInterface private interface LineConsumer{void accept(String line)throws IOException;}
    private static void lines(Path path,TabixReader index,Map<String,String> chromosomes,String chr,long start,long end,LineConsumer consumer)throws IOException {
        if(index!=null) {
            String actual=chromosomes.get(chromosome(chr));if(actual==null)return;
            if(start<1||end>Integer.MAX_VALUE)throw new IOException("tabix coordinates outside supported range");
            TabixReader.Iterator it=index.query(actual,(int)start-1,(int)end);
            if(it!=null)for(String line;(line=it.next())!=null;)consumer.accept(line);
        } else try(BufferedReader r=open(path)) {
            for(String line;(line=r.readLine())!=null;) {
                if(line.isBlank()||line.startsWith("#")||line.startsWith("CHROM"))continue;
                String[] f=fields(line);long pos=Long.parseLong(f[1]);
                if(chromosome(f[0]).equals(chromosome(chr))&&pos>=start&&pos<=end)consumer.accept(line);
            }
        }
    }
    public static BufferedReader open(Path path)throws IOException {
        InputStream stream=Files.newInputStream(path);
        try {return new BufferedReader(new InputStreamReader(path.toString().endsWith(".gz")?new GZIPInputStream(stream):stream,java.nio.charset.StandardCharsets.UTF_8));}
        catch(IOException|RuntimeException e){stream.close();throw e;}
    }
    @Override public void close(){if(scoreIndex!=null)scoreIndex.close();if(covIndex!=null)covIndex.close();}
    public final class Cursor implements AutoCloseable {
        private final BufferedReader reader; private Score previous;
        private Cursor(BufferedReader reader){this.reader=reader;}
        public Score next()throws IOException {
            for(String line;(line=reader.readLine())!=null;) {
                if(line.isBlank()||line.startsWith("#")||line.startsWith("CHROM"))continue;
                Score score=parse(line);
                if(previous!=null && comparePosition(previous,score)>=0)
                    throw new IOException("scores must have unique, sorted biallelic positions: "+score.id());
                previous=score;return score;
            }
            return null;
        }
        @Override public void close()throws IOException{reader.close();}
    }
    /** Natural chromosome ordering: numeric autosomes, X, Y, MT, then other contigs. */
    public static int comparePosition(Score a,Score b){int c=compareChromosome(a.chromosome,b.chromosome);return c!=0?c:Long.compare(a.position,b.position);}
    public static int compareChromosome(String a,String b){int x=rank(a),y=rank(b);return x!=y?Integer.compare(x,y):a.equals(b)?0:a.compareTo(b);}
    private static int rank(String c){try{return Integer.parseInt(c);}catch(NumberFormatException e){return c.equals("X")?1000:c.equals("Y")?1001:c.equals("MT")||c.equals("M")?1002:1003;}}
    public record Score(String chromosome,long position,String reference,String alternate,double samples,double frequency,double score,double variance,String status,double callRate,double hwePValue) {
        /** Source-compatible construction for programmatically supplied summaries. */
        public Score(String chromosome,long position,String reference,String alternate,double samples,double frequency,double score,double variance,String status) {
            this(chromosome,position,reference,alternate,samples,frequency,score,variance,status,Double.NaN,Double.NaN);
        }
        public String id(){return chromosome+":"+position+":"+reference+":"+alternate;}
        public boolean informative(){return status.equals("ok");}
        public int orientation(String ref,String alt){
            if(alternate.isEmpty()&&status.equals("monomorphic"))return reference.equals(ref)?1:reference.equals(alt)?-1:0;
            return reference.equals(ref)&&alternate.equals(alt)?1:reference.equals(alt)&&alternate.equals(ref)?-1:0;
        }
    }
}
