/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import org.jlinalg.network.NetworkAnalysis;
import static org.jlinalg.cli.FollowupSupport.*;

/** Network follow-up CLI: native association methods and audited R adapters. */
final class NetworkCli {
    private NetworkCli() { }
    static int run(String[] args,PrintStream output,PrintStream error) {
        if(Arrays.asList(args).contains("--help")){output.println(help());return 0;}
        try {
            Map<String,String> o=options(args,Set.of("method","out","matrix","matrix-b","edges","hits","directed","restart","radius","permutations","seed","lambda","rule","bootstraps","iterations","tolerance","rscript","r-library","timeout","power","min-module-size","merge-cut-height","edge-threshold","pheno","trait","reference-matrix","regulators","trees","activities","solver","solver-path","max-features","max-samples"));
            String method=required(o,"method");
            Set<String> specific=switch(method) {
                case "reference" -> Set.of("edges","hits","directed","restart","radius","permutations","seed","iterations","tolerance");
                case "sparse" -> Set.of("matrix","lambda","rule","bootstraps","seed","iterations","tolerance","max-features","max-samples");
                case "differential" -> Set.of("matrix","matrix-b","max-features","max-samples");
                case "wgcna" -> Set.of("matrix","rscript","r-library","timeout","seed","power","min-module-size","merge-cut-height","edge-threshold","pheno","trait","reference-matrix","permutations","max-features","max-samples");
                case "regulatory" -> Set.of("matrix","rscript","r-library","timeout","seed","regulators","trees","edge-threshold","max-features","max-samples");
                case "signaling" -> Set.of("edges","activities","rscript","r-library","timeout","solver","solver-path","seed");
                default -> throw new IllegalArgumentException("Unknown network method: "+method);
            };
            for(String key:o.keySet())if(!Set.of("method","out").contains(key)&&!specific.contains(key))throw new IllegalArgumentException("--"+key+" is not applicable to "+method);
            transaction(o,"network-"+method,(dir,manifest)-> {
                switch(method) {
                    case "reference" -> reference(o,dir,manifest);
                    case "sparse" -> sparse(o,dir,manifest);
                    case "differential" -> differential(o,dir,manifest);
                    default -> adapter(method,o,dir,manifest);
                }
                graphml(o,dir,method);
            });
            output.println("network "+method+" completed: "+o.get("out"));return 0;
        } catch(Exception failure){error.println("jlinalg: "+failure.getMessage());return 2;}
    }
    static String help(){return """
        Network follow-up; --out is a NEW directory; all methods write manifests.
          network --method reference --edges edges.tsv --hits hits.tsv --out DIR
            [--directed true|false --restart 0.5 --radius 1 --permutations 999 --seed 1]
          network --method sparse --matrix matrix.tsv --lambda 0.1 --out DIR
            [--rule and|or --bootstraps 0 --seed 1 --iterations 10000 --tolerance 1e-8]
          network --method differential --matrix group-a.tsv --matrix-b group-b.tsv --out DIR
          network --method wgcna --matrix matrix.tsv --power 6 --out DIR
            [--min-module-size 30 --merge-cut-height 0.25 --edge-threshold 0.1]
            [--pheno phenotype.tsv --trait COLUMN --reference-matrix replication.tsv --permutations 100]
          network --method regulatory --matrix expression.tsv --regulators regulators.tsv --out DIR
            [--trees 1000 --edge-threshold 0 --seed 1]
          network --method signaling --edges signed-directed.tsv --activities activities.tsv --out DIR
            [--solver lpSolve|cbc|cplex --solver-path EXECUTABLE]
        R adapters: --rscript EXECUTABLE [--r-library DIRECTORY --timeout 3600].
        wgcna uses WGCNA; regulatory uses GENIE3; signaling uses inverse CARNIVAL.
        Matrix: sample_id followed by numerical feature columns, independent rows, no missing values.
        Networks from matrices require --max-features (default 2000, hard limit 2000) and
        --max-samples (default 10000) bounds. Normalize/adjust covariates upstream.
        Reference edges: source,target,weight,type; no duplicate pairs/self edges; positive weights.
        Hits/regulators: gene column. Signaling edges: source,target,sign (-1 or 1).
        Activities: gene,activity (-1 or 1), representing regulator activities, not raw expression.
        Association edges, hubs and inferred regulators are not experimentally validated drivers.
        Each result includes network.graphml for Cytoscape or another GraphML viewer.
        """;}
    private static void graphml(Map<String,String> o,Path dir,String method) throws Exception {
        Path edgeFile=dir.resolve("edges.tsv");List<String[]> rows=new ArrayList<>();List<String> header;
        try(var in=Files.newBufferedReader(edgeFile)) {
            header=DelimitedData.parse(in.readLine(),'\t',1,edgeFile);
            for(String line;(line=in.readLine())!=null;)if(!line.isBlank())rows.add(DelimitedData.parse(line,'\t',0,edgeFile).toArray(String[]::new));
        }
        int s=header.indexOf("source"),t=header.indexOf("target");if(s<0||t<0)throw new IOException("Edges require canonical source/target columns");
        Map<String,Integer> nodes=new LinkedHashMap<>();
        if(o.containsKey("matrix")) {
            Path matrix=Path.of(o.get("matrix"));
            try(var reader=Files.newBufferedReader(matrix)) {
                List<String> names=DelimitedData.parse(reader.readLine(),matrix.toString().toLowerCase(Locale.ROOT).endsWith(".csv")?',':'\t',1,matrix);
                for(String gene:names.subList(1,names.size()))nodes.put(gene,nodes.size());
            }
        }
        for(String[] row:rows){nodes.computeIfAbsent(row[s],key->nodes.size());nodes.computeIfAbsent(row[t],key->nodes.size());}
        boolean directed=Set.of("regulatory","signaling").contains(method)||method.equals("reference")&&bool(o,"directed",false);
        try(var stream=Files.newOutputStream(dir.resolve("network.graphml"))) {
            var xml=javax.xml.stream.XMLOutputFactory.newFactory().createXMLStreamWriter(stream,"UTF-8");
            xml.writeStartDocument("UTF-8","1.0");xml.writeStartElement("graphml");xml.writeDefaultNamespace("http://graphml.graphdrawing.org/xmlns");
            xml.writeEmptyElement("key");xml.writeAttribute("id","label");xml.writeAttribute("for","node");xml.writeAttribute("attr.name","gene");xml.writeAttribute("attr.type","string");
            for(int i=0;i<header.size();i++)if(i!=s&&i!=t){xml.writeEmptyElement("key");xml.writeAttribute("id","k"+i);xml.writeAttribute("for","edge");xml.writeAttribute("attr.name",header.get(i));xml.writeAttribute("attr.type","string");}
            xml.writeStartElement("graph");xml.writeAttribute("id","network");xml.writeAttribute("edgedefault",directed?"directed":"undirected");
            for(var entry:nodes.entrySet()){xml.writeStartElement("node");xml.writeAttribute("id","n"+entry.getValue());xml.writeStartElement("data");xml.writeAttribute("key","label");xml.writeCharacters(entry.getKey());xml.writeEndElement();xml.writeEndElement();}
            int id=0;for(String[] row:rows){xml.writeStartElement("edge");xml.writeAttribute("id","e"+id++);xml.writeAttribute("source","n"+nodes.get(row[s]));xml.writeAttribute("target","n"+nodes.get(row[t]));for(int i=0;i<header.size();i++)if(i!=s&&i!=t){xml.writeStartElement("data");xml.writeAttribute("key","k"+i);xml.writeCharacters(row[i]);xml.writeEndElement();}xml.writeEndElement();}
            xml.writeEndElement();xml.writeEndElement();xml.writeEndDocument();xml.close();
        }
    }
    record Matrix(List<String> ids,List<String> genes,double[][] values) { }
    private static Matrix matrix(Path file,Map<String,String> o) throws IOException {
        DelimitedData data=DelimitedData.read(file);if(!data.header().get(0).equals("sample_id"))throw new IllegalArgumentException("Matrix first column must be sample_id");
        int p=data.header().size()-1,n=data.rows().size();
        if(p<2||p>integer(o,"max-features",2000,2,2000)||n<4||n>integer(o,"max-samples",10000,4,1000000)||(long)n*p>20_000_000L)throw new IllegalArgumentException("Matrix exceeds sample/feature/work bounds or has <4 samples / <2 features");
        List<String> ids=new ArrayList<>();Set<String> seen=new HashSet<>();double[][] values=new double[n][p];
        for(int i=0;i<n;i++){String[] row=data.rows().get(i);if(row[0].isBlank()||!seen.add(row[0]))throw new IllegalArgumentException("Duplicate/blank sample_id");ids.add(row[0]);for(int j=0;j<p;j++)values[i][j]=finite(row[j+1],"matrix value");}
        NetworkAnalysis.standardize(values);return new Matrix(ids,data.header().subList(1,p+1),values);
    }
    private static void writeMatrix(Matrix matrix,Path file) throws IOException {
        List<String> header=new ArrayList<>(List.of("sample_id"));header.addAll(matrix.genes);List<String[]> rows=new ArrayList<>();
        for(int i=0;i<matrix.ids.size();i++){String[] row=new String[header.size()];row[0]=matrix.ids.get(i);for(int j=0;j<matrix.genes.size();j++)row[j+1]=Double.toString(matrix.values[i][j]);rows.add(row);}table(file,header,rows);
    }
    private static Set<String> genes(Path path) throws IOException {
        DelimitedData data=DelimitedData.read(path);int col=data.column("gene");Set<String> values=new LinkedHashSet<>();
        for(String[] row:data.rows())if(row[col].isBlank()||!values.add(row[col]))throw new IllegalArgumentException("Duplicate/blank gene");return values;
    }
    private static void reference(Map<String,String> o,Path dir,Map<String,Object> manifest) throws Exception {
        DelimitedData edges=DelimitedData.read(Path.of(required(o,"edges")));int s=edges.column("source"),t=edges.column("target"),w=edges.column("weight"),type=edges.column("type");
        boolean directed=bool(o,"directed",false);Map<String,Integer> index=new LinkedHashMap<>();Set<String> pairs=new HashSet<>();
        for(String[] row:edges.rows()) {
            if(row[s].isBlank()||row[t].isBlank()||row[s].equals(row[t])||row[type].isBlank()||finite(row[w],"edge weight")<=0)throw new IllegalArgumentException("Require distinct nonblank nodes, edge type and positive weights");
            String pair=directed||row[s].compareTo(row[t])<0?row[s]+"\0"+row[t]:row[t]+"\0"+row[s];
            if(!pairs.add(pair))throw new IllegalArgumentException("Duplicate edge pair; select/aggregate an evidence layer explicitly");
            index.computeIfAbsent(row[s],k->index.size());index.computeIfAbsent(row[t],k->index.size());
        }
        if(index.size()>100000||edges.rows().size()>1000000)throw new IllegalArgumentException("Reference network exceeds 100000 nodes / 1000000 edges");
        Set<String> hits=genes(Path.of(required(o,"hits")));List<String[]> unmapped=new ArrayList<>();for(String hit:hits)if(!index.containsKey(hit))unmapped.add(new String[]{hit,"absent_from_network"});
        table(dir.resolve("unmapped.tsv"),List.of("gene","status"),unmapped);
        int n=index.size();List<Map<Integer,Double>> graph=new ArrayList<>();for(int i=0;i<n;i++)graph.add(new LinkedHashMap<>());
        for(String[] row:edges.rows()){int a=index.get(row[s]),b=index.get(row[t]);double weight=Double.parseDouble(row[w]);graph.get(a).put(b,weight);if(!directed)graph.get(b).put(a,weight);}
        double[] seeds=new double[n];for(String hit:hits)if(index.containsKey(hit))seeds[index.get(hit)]=1;
        double[] score=NetworkAnalysis.propagate(graph,seeds,number(o,"restart",.5),number(o,"tolerance",1e-10),integer(o,"iterations",10000,1,1000000));
        int radius=integer(o,"radius",1,1,3),permutations=integer(o,"permutations",999,0,100000),seed=integer(o,"seed",1,0,Integer.MAX_VALUE);
        List<Set<Integer>> neighborhoods=new ArrayList<>();long entries=0;
        for(int i=0;i<n;i++) {
            Set<Integer> seen=new LinkedHashSet<>(Set.of(i)),front=new LinkedHashSet<>(seen);
            for(int r=0;r<radius;r++){Set<Integer> next=new LinkedHashSet<>();for(int node:front)next.addAll(graph.get(node).keySet());next.removeAll(seen);seen.addAll(next);front=next;}
            seen.remove(i);entries+=seen.size();if(entries*(permutations+1L)>100_000_000L)throw new IllegalArgumentException("Neighborhood permutation work exceeds 100 million visits; reduce radius/permutations/network");neighborhoods.add(seen);
        }
        int[] observed=new int[n],exceed=new int[n];for(int i=0;i<n;i++)for(int j:neighborhoods.get(i))if(seeds[j]>0)observed[i]++;
        Map<Integer,List<Integer>> bins=new TreeMap<>();for(int i=0;i<n;i++){int degree=graph.get(i).size();int bin=degree==0?0:32-Integer.numberOfLeadingZeros(degree);bins.computeIfAbsent(bin,k->new ArrayList<>()).add(i);}
        Random random=new Random(seed);
        for(int b=0;b<permutations;b++) {
            boolean[] shuffled=new boolean[n];
            for(List<Integer> bin:bins.values()){int count=0;for(int i:bin)if(seeds[i]>0)count++;List<Integer> copy=new ArrayList<>(bin);Collections.shuffle(copy,random);for(int j=0;j<count;j++)shuffled[copy.get(j)]=true;}
            for(int i=0;i<n;i++){int count=0;for(int j:neighborhoods.get(i))if(shuffled[j])count++;if(count>=observed[i])exceed[i]++;}
        }
        double[] p=new double[n];for(int i=0;i<n;i++)p[i]=(exceed[i]+1.0)/(permutations+1.0);double[] q=NetworkAnalysis.bh(p);
        List<String[]> rows=new ArrayList<>();for(var e:index.entrySet()){int i=e.getValue();rows.add(new String[]{e.getKey(),Boolean.toString(seeds[i]>0),Integer.toString(graph.get(i).size()),Double.toString(score[i]),Integer.toString(neighborhoods.get(i).size()),Integer.toString(observed[i]),permutations>0?Double.toString(p[i]):"NA",permutations>0?Double.toString(q[i]):"NA"});}
        table(dir.resolve("nodes.tsv"),List.of("gene","seed","degree","propagation_score","neighborhood_size","neighborhood_hits","driver_p","driver_bh"),rows);
        table(dir.resolve("edges.tsv"),edges.header(),edges.rows());
        manifest.put("interpretation","Propagation is a ranking. Candidate-driver p values permute seed labels within floor(log2(out-degree)) bins, excluding each candidate from its neighborhood; BH covers all network nodes. No causal identification.");
        manifest.put("unmapped_hits",unmapped.size());manifest.put("resolved",Map.of("directed",directed,"radius",radius,"permutations",permutations,"seed",seed,"restart",number(o,"restart",.5)));
    }
    private static boolean selected(double a,double b,String rule){return rule.equals("and")?Math.abs(a)>1e-10&&Math.abs(b)>1e-10:Math.abs(a)>1e-10||Math.abs(b)>1e-10;}
    private static void sparse(Map<String,String> o,Path dir,Map<String,Object> manifest) throws Exception {
        Matrix m=matrix(Path.of(required(o,"matrix")),o);double lambda=Double.parseDouble(required(o,"lambda"));String rule=o.getOrDefault("rule","and");
        if(!Set.of("and","or").contains(rule))throw new IllegalArgumentException("rule must be and or or");
        int iterations=integer(o,"iterations",10000,1,1000000),boots=integer(o,"bootstraps",0,0,1000),seed=integer(o,"seed",1,0,Integer.MAX_VALUE);double tolerance=number(o,"tolerance",1e-8);
        int n=m.ids.size(),p=m.genes.size();if((long)n*p*p*(boots+1)>500_000_000L)throw new IllegalArgumentException("Sparse network work exceeds 500 million sample-feature-pair visits");
        double[][] beta=NetworkAnalysis.neighborhoodLasso(m.values,lambda,iterations,tolerance);int[][] counts=new int[p][p];Random random=new Random(seed);
        for(int b=0;b<boots;b++){double[][] sample=new double[n][];for(int i=0;i<n;i++)sample[i]=m.values[random.nextInt(n)].clone();double[][] fit=NetworkAnalysis.neighborhoodLasso(sample,lambda,iterations,tolerance);for(int i=0;i<p;i++)for(int j=i+1;j<p;j++)if(selected(fit[i][j],fit[j][i],rule))counts[i][j]++;}
        List<String[]> rows=new ArrayList<>();for(int i=0;i<p;i++)for(int j=i+1;j<p;j++)if(selected(beta[i][j],beta[j][i],rule))rows.add(new String[]{m.genes.get(i),m.genes.get(j),Double.toString(beta[j][i]),Double.toString(beta[i][j]),Double.toString((Math.abs(beta[j][i])+Math.abs(beta[i][j]))/2),"conditional_association",boots>0?Double.toString((double)counts[i][j]/boots):"NA"});
        table(dir.resolve("edges.tsv"),List.of("source","target","target_on_source","source_on_target","weight","type","bootstrap_frequency"),rows);
        manifest.put("resolved",Map.of("lambda",lambda,"rule",rule,"bootstraps",boots,"seed",seed,"iterations",iterations,"tolerance",tolerance));manifest.put("interpretation","Standardized Gaussian neighborhood LASSO; edge magnitudes are regression weights, not partial correlations or p values. Bootstrap frequency is stability, not posterior probability.");
    }
    private static Matrix reorder(Matrix second,List<String> genes) {
        if(!new HashSet<>(second.genes).equals(new HashSet<>(genes)))throw new IllegalArgumentException("Matrices must have identical feature sets");
        double[][] values=new double[second.ids.size()][genes.size()];for(int j=0;j<genes.size();j++){int source=second.genes.indexOf(genes.get(j));for(int i=0;i<values.length;i++)values[i][j]=second.values[i][source];}return new Matrix(second.ids,genes,values);
    }
    private static void independent(Matrix a,Matrix b){Set<String> overlap=new HashSet<>(a.ids);overlap.retainAll(b.ids);if(!overlap.isEmpty())throw new IllegalArgumentException("Independent cohorts required: overlapping sample IDs (paired/repeated data unsupported)");}
    private static void differential(Map<String,String> o,Path dir,Map<String,Object> manifest) throws Exception {
        Matrix a=matrix(Path.of(required(o,"matrix")),o),b=reorder(matrix(Path.of(required(o,"matrix-b")),o),a.genes);independent(a,b);
        double[][] ca=NetworkAnalysis.correlation(a.values),cb=NetworkAnalysis.correlation(b.values);List<String[]> rows=new ArrayList<>();List<Double> pv=new ArrayList<>();
        double se=Math.sqrt(1.0/(a.ids.size()-3)+1.0/(b.ids.size()-3));
        for(int i=0;i<a.genes.size();i++)for(int j=i+1;j<a.genes.size();j++) {
            double ra=ca[i][j],rb=cb[i][j];boolean valid=Math.abs(ra)<1-1e-12&&Math.abs(rb)<1-1e-12;
            double z=valid?(.5*Math.log((1+ra)/(1-ra))-.5*Math.log((1+rb)/(1-rb)))/se:Double.NaN;
            double p=valid?Math.min(1,2*jdistlib.Normal.cumulative(Math.abs(z),0,1,false,false)):1;
            pv.add(p);rows.add(new String[]{a.genes.get(i),a.genes.get(j),Double.toString(ra),Double.toString(rb),Double.toString(ra-rb),valid?Double.toString(z):"NA",valid?Double.toString(p):"NA","NA",valid?"ok":"degenerate_correlation"});
        }
        double[] q=NetworkAnalysis.bh(pv.stream().mapToDouble(Double::doubleValue).toArray());for(int i=0;i<rows.size();i++)if(rows.get(i)[8].equals("ok"))rows.get(i)[7]=Double.toString(q[i]);
        table(dir.resolve("edges.tsv"),List.of("source","target","r_a","r_b","difference","z","p","bh","status"),rows);
        manifest.put("interpretation","Independent-group Fisher-z tests assume approximately bivariate-normal observations. BH family includes every feature pair (degenerate pairs counted at p=1, reported NA). No paired/clustered inference.");manifest.put("sample_counts",List.of(a.ids.size(),b.ids.size()));
    }
    private static void adapter(String method,Map<String,String> o,Path dir,Map<String,Object> manifest) throws Exception {
        Map<String,String> settings=new LinkedHashMap<>(o);settings.putIfAbsent("seed","1");integer(settings,"seed",1,0,Integer.MAX_VALUE);
        if(method.equals("wgcna")||method.equals("regulatory")) {
            Matrix m=matrix(Path.of(required(o,"matrix")),o);writeMatrix(m,dir.resolve("matrix.tsv"));
            if(method.equals("wgcna")) {
                if(m.ids.size()<15)throw new IllegalArgumentException("WGCNA requires at least 15 independent samples (20+ recommended)");
                double power=Double.parseDouble(required(o,"power"));if(!Double.isFinite(power)||power<=0)throw new IllegalArgumentException("power must be positive");
                settings.putIfAbsent("min-module-size","30");integer(settings,"min-module-size",30,2,m.genes.size());
                settings.putIfAbsent("merge-cut-height","0.25");double cut=number(settings,"merge-cut-height",.25);if(cut<0||cut>1)throw new IllegalArgumentException("merge-cut-height must be in [0,1]");
                settings.putIfAbsent("permutations","100");integer(settings,"permutations",100,10,10000);
                if(o.containsKey("pheno")!=o.containsKey("trait"))throw new IllegalArgumentException("pheno and trait must be supplied together");
                if(o.containsKey("pheno")) {
                    DelimitedData pheno=DelimitedData.read(Path.of(o.get("pheno")));int id=pheno.column("sample_id"),trait=pheno.column(o.get("trait"));Map<String,String> values=new HashMap<>();
                    for(String[] row:pheno.rows()){finite(row[trait],"trait");if(values.put(row[id],row[trait])!=null)throw new IllegalArgumentException("Duplicate phenotype sample ID");}
                    List<String[]> rows=new ArrayList<>();for(String sample:m.ids){if(!values.containsKey(sample))throw new IllegalArgumentException("Missing phenotype sample: "+sample);rows.add(new String[]{sample,values.get(sample)});}table(dir.resolve("pheno.tsv"),List.of("sample_id","trait"),rows);
                }
                if(o.containsKey("reference-matrix")){Matrix ref=reorder(matrix(Path.of(o.get("reference-matrix")),o),m.genes);independent(m,ref);if(ref.ids.size()<15)throw new IllegalArgumentException("Preservation requires >=15 independent replication samples");writeMatrix(ref,dir.resolve("reference.tsv"));}
            } else {
                Set<String> regulators=genes(Path.of(required(o,"regulators")));if(!m.genes.containsAll(regulators))throw new IllegalArgumentException("All regulators must be measured features");
                List<String[]> rows=regulators.stream().map(s->new String[]{s}).toList();table(dir.resolve("regulators.tsv"),List.of("gene"),rows);
                settings.putIfAbsent("trees","1000");integer(settings,"trees",1000,1,100000);
            }
            settings.putIfAbsent("edge-threshold",method.equals("wgcna")?"0.1":"0");double threshold=number(settings,"edge-threshold",0);if(threshold<0||threshold>1)throw new IllegalArgumentException("edge-threshold must be in [0,1]");
        } else {
            DelimitedData edges=DelimitedData.read(Path.of(required(o,"edges")));int s=edges.column("source"),t=edges.column("target"),sign=edges.column("sign");Set<String> nodes=new HashSet<>(),pairs=new HashSet<>();List<String[]> rows=new ArrayList<>();
            for(String[] row:edges.rows()){if(!Set.of("-1","1").contains(row[sign])||row[s].isBlank()||row[t].isBlank()||row[s].equals(row[t])||!pairs.add(row[s]+"\0"+row[t]))throw new IllegalArgumentException("Signaling requires unique signed directed edges without self loops");nodes.add(row[s]);nodes.add(row[t]);rows.add(new String[]{row[s],row[sign],row[t]});}table(dir.resolve("prior.tsv"),List.of("source","interaction","target"),rows);
            DelimitedData acts=DelimitedData.read(Path.of(required(o,"activities")));int gene=acts.column("gene"),act=acts.column("activity");Set<String> seen=new HashSet<>();rows=new ArrayList<>();
            for(String[] row:acts.rows()){if(!nodes.contains(row[gene])||!seen.add(row[gene])||!Set.of("-1","1").contains(row[act]))throw new IllegalArgumentException("Activities require unique prior-network genes and activity +/-1");rows.add(new String[]{row[gene],row[act]});}table(dir.resolve("activities.tsv"),List.of("gene","activity"),rows);
            settings.putIfAbsent("solver","lpSolve");if(!Set.of("lpSolve","cbc","cplex").contains(settings.get("solver")))throw new IllegalArgumentException("solver must be lpSolve, cbc or cplex");
            if(!settings.get("solver").equals("lpSolve"))settings.put("solver-path",Path.of(required(o,"solver-path")).toAbsolutePath().toString());
        }
        if(o.containsKey("r-library"))settings.put("r-library",Path.of(o.get("r-library")).toAbsolutePath().toString());
        List<String[]> values=settings.entrySet().stream().map(e->new String[]{e.getKey(),e.getValue()}).toList();table(dir.resolve("settings.tsv"),List.of("key","value"),values);
        Path script=dir.resolve("adapter.R");try(InputStream in=NetworkCli.class.getResourceAsStream("/org/jlinalg/network/"+method+".R")){if(in==null)throw new IOException("Missing packaged adapter");Files.copy(in,script);}
        List<String> command=List.of(executable(required(o,"rscript")),"--vanilla",script.toString());manifest.put("external_command",command);manifest.put("resolved",settings);
        process(command,dir,integer(o,"timeout",3600,1,86400));
        for(String file:List.of("edges.tsv","sessionInfo.txt"))if(!Files.isRegularFile(dir.resolve(file)))throw new IOException("Adapter did not produce "+file);
        manifest.put("interpretation","External adapter output; inspect engine diagnostics and sessionInfo.txt. Inferred modules, regulatory edges and contextual signaling are hypotheses, not experimentally validated drivers.");
    }
}
