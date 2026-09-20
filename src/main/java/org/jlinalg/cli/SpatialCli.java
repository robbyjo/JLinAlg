/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import org.jlinalg.spatial.SpatialStatistics;
import org.jlinalg.spatial.SpatialStatistics.Edge;
import org.jlinalg.singlecell.SampleInference;
import org.jlinalg.network.NetworkAnalysis;
import static org.jlinalg.cli.FollowupSupport.*;

/** Physical tissue graphs and independent-sample spatial comparisons. */
final class SpatialCli {
    private SpatialCli() { }
    static int run(String[] args,PrintStream out,PrintStream err) {
        if(Arrays.asList(args).contains("--help")){out.println(help());return 0;}
        try {
            Set<String> allowed=new HashSet<>(SingleCellCli.INPUT);allowed.addAll(SingleCellCli.DESIGN);
            allowed.addAll(Set.of("units","resolution","coordinate-system","graph","radius","k","edges","permutations","seed","feature-list","distance-column"));
            Map<String,String> o=options(args,allowed);String method=required(o,"method");
            if(!Set.of("graph","autocorrelation","neighborhoods","compare","gradient").contains(method))throw new IllegalArgumentException("Unknown spatial method: "+method);
            Set<String> specific=new HashSet<>(Set.of("units","resolution","coordinate-system"));
            if(!method.equals("gradient"))specific.addAll(Set.of("graph","radius","k","edges"));
            if(Set.of("autocorrelation","neighborhoods").contains(method))specific.addAll(Set.of("permutations","seed"));
            if(Set.of("compare","gradient").contains(method))specific.addAll(SingleCellCli.DESIGN);
            if(Set.of("autocorrelation","gradient").contains(method))specific.add("feature-list");
            if(method.equals("gradient"))specific.add("distance-column");
            for(String key:o.keySet())if(!SingleCellCli.INPUT.contains(key) && !specific.contains(key))throw new IllegalArgumentException("--"+key+" is not applicable to "+method);
            transaction(o,"spatial-"+method,(dir,manifest)->execute(method,o,dir,manifest));
            out.println("spatial "+method+" completed: "+o.get("out"));return 0;
        }catch(Exception failure){err.println("jlinalg: "+failure.getMessage());return 2;}
    }
    static String help(){return """
        spatial --method graph|autocorrelation|neighborhoods|compare|gradient
          --counts sparse.tsv --cells cells.tsv --samples samples.tsv --features features.tsv
          --units micrometer --coordinate-system registered-v1 --resolution cell|spot --out NEW_DIR
        Shared single-cell count/QC contract, plus cell columns:
          specimen_id,section_id,compartment,stratum,x,y [z]. Supplied coordinates are already transformed.
        Graph: --graph radius|knn|adjacency --radius POSITIVE_DISTANCE [--k 6]
          adjacency instead requires --edges TSV (source,target); same section/compartment only.
        Radius is also a maximum gap for kNN; undirected union graph. No cross-section edges.
        Autocorrelation: per-section Moran I / Geary C on log1p(10000 counts/library).
          --permutations 999 --seed 1 [--feature-list TSV with feature_id].
        Neighborhoods: random-label enrichment/depletion within section/compartment/stratum.
        Compare: sample-level edge-pair fractions, --reference A --tested B [--paired true --covariates age].
          Cell-resolved inputs only for neighborhood/compare; inferred spot mixtures are unsupported.
        Gradient: --distance-column COLUMN [--feature-list TSV] and the same inference options;
          compare descriptive within-sample linear expression/distance slopes across independent donors.
        Maps, edges, QC, sample identities, seeds, input hashes and null definitions are exported.
        """;}
    record Geometry(List<Integer> original,double[][] xyz,String[] section,String[] compartment,String[] strata) { }
    private static Geometry geometry(CellData data,boolean[] keep) {
        List<Integer> original=new ArrayList<>();for(int i=0;i<keep.length;i++)if(keep[i])original.add(i);
        int n=original.size(),dim=data.cells.header().contains("z")?3:2;
        double[][] xyz=new double[n][dim];String[] section=new String[n],comp=new String[n],strata=new String[n];
        Map<String,String> specimenSamples=new HashMap<>(),sectionSpecimens=new HashMap<>();
        for(int i=0;i<n;i++) {
            int source=original.get(i);String sample=data.cell(source,"sample_id"),specimen=data.cell(source,"specimen_id"),slice=data.cell(source,"section_id");
            String region=data.cell(source,"compartment"),stratum=data.cell(source,"stratum");
            if(List.of(specimen,slice,region,stratum).stream().anyMatch(String::isBlank))throw new IllegalArgumentException("Blank spatial hierarchy/stratum");
            String prior=specimenSamples.putIfAbsent(specimen,sample);if(prior!=null&&!prior.equals(sample))throw new IllegalArgumentException("specimen_id must belong to exactly one sample");
            prior=sectionSpecimens.putIfAbsent(slice,specimen);if(prior!=null&&!prior.equals(specimen))throw new IllegalArgumentException("section_id must be globally unique and belong to one specimen");
            section[i]=slice;comp[i]=slice+"\0"+region;strata[i]=comp[i]+"\0"+stratum;
            xyz[i][0]=finite(data.cell(source,"x"),"x");xyz[i][1]=finite(data.cell(source,"y"),"y");if(dim==3)xyz[i][2]=finite(data.cell(source,"z"),"z");
        }
        return new Geometry(original,xyz,section,comp,strata);
    }
    private static List<Edge> graph(CellData data,Geometry g,Map<String,String> o) throws IOException {
        String method=o.getOrDefault("graph","radius");
        if(method.equals("adjacency")) {
            if(o.containsKey("radius") || o.containsKey("k"))throw new IllegalArgumentException("Adjacency uses supplied edges; radius/k not applicable");
            DelimitedData input=CellData.bounded(Path.of(required(o,"edges")));int sc=input.column("source"),tc=input.column("target");
            Map<String,Integer> index=new HashMap<>();for(int i=0;i<g.original.size();i++)index.put(data.ids.get(g.original.get(i)),i);
            Set<String> all=new HashSet<>(data.ids);Set<List<String>> seen=new HashSet<>();List<Edge> result=new ArrayList<>();
            for(String[] row:input.rows()) {
                String a=row[sc],b=row[tc];if(a.compareTo(b)>0){String swap=a;a=b;b=swap;}
                if(!all.contains(a)||!all.contains(b)||a.equals(b)||!seen.add(List.of(a,b)))throw new IllegalArgumentException("Unknown/self/duplicate adjacency edge");
                Integer s=index.get(a),t=index.get(b);if(s==null||t==null)continue; // audit QC output explains removed nodes
                if(!g.compartment[s].equals(g.compartment[t]))throw new IllegalArgumentException("Adjacency edge crosses section/compartment");
                double d=0;for(int j=0;j<g.xyz[s].length;j++)d=Math.hypot(d,g.xyz[s][j]-g.xyz[t][j]);
                if(!Double.isFinite(d))throw new IllegalArgumentException("Adjacency distance overflows; rescale coordinates");
                result.add(new Edge(Math.min(s,t),Math.max(s,t),d));
                if(result.size()>1_000_000)throw new IllegalArgumentException("Graph exceeds one million edges");
            }
            return result;
        }
        if(o.containsKey("edges") || method.equals("radius")&&o.containsKey("k"))throw new IllegalArgumentException("Edges/k not applicable to graph method");
        return SpatialStatistics.graph(g.xyz,g.compartment,method,finite(required(o,"radius"),"radius"),integer(o,"k",6,1,1000));
    }
    private static void execute(String method,Map<String,String> o,Path dir,Map<String,Object> manifest) throws Exception {
        String units=required(o,"units"),resolution=required(o,"resolution");required(o,"coordinate-system");
        if(!Set.of("micrometer","millimeter").contains(units)||!Set.of("cell","spot").contains(resolution))throw new IllegalArgumentException("Require physical units micrometer/millimeter and cell/spot resolution");
        if(Set.of("neighborhoods","compare").contains(method)&&!resolution.equals("cell"))throw new IllegalArgumentException("Neighborhood cell-label inference requires cell-resolved observations");
        CellData data=CellData.read(o);
        if(Set.of("compare","gradient").contains(method))CellDesign.validate(data,o);
        boolean[] keep=data.qc(o,dir,manifest);Geometry g=geometry(data,keep);
        if(g.original.size()<3)throw new IllegalArgumentException("Fewer than three retained spatial observations");
        manifest.put("spatial_contract",Map.of("units",units,"resolution",resolution,"coordinate_system",o.get("coordinate-system"),"registration","supplied, no transformation estimated","boundary_policy","no cross-section or cross-compartment edges; radius caps gaps; interior masks require supplied adjacency"));
        if(method.equals("gradient")){gradient(data,g,o,dir,manifest);return;}
        List<Edge> edges=graph(data,g,o);Map<String,List<Edge>> sectionEdges=partitionEdges(g,edges);
        manifest.put("graph",Map.of("method",o.getOrDefault("graph","radius"),"edges",edges.size(),"observations",g.original.size(),"weight","binary symmetric, each unordered pair exported once"));
        switch(method) {
            case "autocorrelation" -> autocorrelation(data,g,edges,sectionEdges,o,dir,manifest);
            case "neighborhoods" -> neighborhoods(data,g,edges,sectionEdges,o,dir,manifest);
            case "compare" -> compare(data,g,edges,o,dir,manifest);
            default -> { }
        }
        writeGraph(data,g,edges,sectionEdges,dir,units);
    }
    private static List<Integer> features(CellData data,Map<String,String> o) throws IOException {
        if(!o.containsKey("feature-list")){List<Integer> result=new ArrayList<>();for(int j=0;j<data.genes.size();j++)result.add(j);return result;}
        DelimitedData table=CellData.bounded(Path.of(o.get("feature-list")));Map<String,Integer> unique=CellData.unique(table,"feature_id");List<Integer> out=new ArrayList<>();
        for(String name:unique.keySet()){int j=data.featureIndex.getOrDefault(name,-1);if(j<0)throw new IllegalArgumentException("Feature not measured: "+name);out.add(j);}return out;
    }
    private static Map<String,List<Integer>> sections(Geometry g) {
        Map<String,List<Integer>> groups=new LinkedHashMap<>();for(int i=0;i<g.section.length;i++)groups.computeIfAbsent(g.section[i],s->new ArrayList<>()).add(i);return groups;
    }
    private static Map<String,List<Edge>> partitionEdges(Geometry g,List<Edge> edges) {
        Map<String,List<Edge>> result=new HashMap<>();
        for(String section:g.section)result.computeIfAbsent(section,key->new ArrayList<>());
        for(Edge edge:edges)result.get(g.section[edge.source()]).add(edge);
        return result;
    }
    private static void autocorrelation(CellData data,Geometry g,List<Edge> edges,Map<String,List<Edge>> sectionEdges,Map<String,String> o,Path dir,Map<String,Object> manifest) throws IOException {
        List<Integer> features=features(data,o);int b=integer(o,"permutations",999,1,100000),seed=integer(o,"seed",1,0,Integer.MAX_VALUE);
        CellData.checkOutputSize((long)sectionEdges.size()*features.size());
        if((long)(edges.size()+g.original.size())*features.size()*b>100_000_000)throw new IllegalArgumentException("Spatial feature permutations exceed 100 million visits; select features or reduce permutations");
        List<String[]> rows=new ArrayList<>();List<Double> pv=new ArrayList<>();
        for(var entry:sections(g).entrySet()) {
            List<Integer> nodes=entry.getValue();Map<Integer,Integer> index=new HashMap<>();for(int i=0;i<nodes.size();i++)index.put(nodes.get(i),i);
            List<Edge> local=sectionEdges.get(entry.getKey()).stream().map(e->new Edge(index.get(e.source()),index.get(e.target()),e.distance())).toList();
            String[] strata=nodes.stream().map(i->g.strata[i]).toArray(String[]::new);
            for(int j:features) {
                double[] values=nodes.stream().mapToDouble(i->data.normalized(g.original.get(i),j)).toArray();
                try {
                    var fit=SpatialStatistics.test(values,local,strata,b,seed);
                    rows.add(new String[]{entry.getKey(),data.genes.get(j),Integer.toString(nodes.size()),Integer.toString(local.size()),Double.toString(fit.moran()),Double.toString(fit.geary()),Double.toString(fit.moranP()),Double.toString(fit.gearyP()),"NA","NA","tested"});pv.add(fit.moranP());pv.add(fit.gearyP());
                }catch(IllegalArgumentException insufficient){rows.add(new String[]{entry.getKey(),data.genes.get(j),Integer.toString(nodes.size()),Integer.toString(local.size()),"NA","NA","NA","NA","NA","NA",insufficient.getMessage()});pv.add(1.0);pv.add(1.0);}
            }
        }
        double[] q=NetworkAnalysis.bh(pv.stream().mapToDouble(Double::doubleValue).toArray());
        for(int i=0;i<rows.size();i++)if(rows.get(i)[10].equals("tested")){rows.get(i)[8]=Double.toString(q[2*i]);rows.get(i)[9]=Double.toString(q[2*i+1]);}
        table(dir.resolve("autocorrelation.tsv"),List.of("section_id","feature_id","observations","edges","moran_i","geary_c","moran_p","geary_p","moran_bh","geary_bh","status"),rows);
        manifest.put("null","Random exchangeable normalized feature labels within section/compartment/stratum. Two-sided |I+1/(n-1)| and |C-1|; p=(exceedances+1)/(B+1). BH jointly covers both statistics, all sections and all requested features. Pattern within a section is not a replicated condition effect.");
        manifest.put("permutations",b);manifest.put("seed",seed);
    }
    private static List<String> types(CellData data,Geometry g) {
        Set<String> result=new TreeSet<>();for(int i:g.original)result.add(data.cell(i,"cell_type"));
        if(result.size()>50)throw new IllegalArgumentException("Neighborhood inference limited to 50 labels");return new ArrayList<>(result);
    }
    private static int[] pairCounts(double[] labels,List<Edge> edges,int k) {
        int[] counts=new int[k*k];for(Edge e:edges){int a=(int)labels[e.source()],b=(int)labels[e.target()];counts[Math.min(a,b)*k+Math.max(a,b)]++;}return counts;
    }
    private static void neighborhoods(CellData data,Geometry g,List<Edge> edges,Map<String,List<Edge>> sectionEdges,Map<String,String> o,Path dir,Map<String,Object> manifest) throws IOException {
        List<String> types=types(data,g);int k=types.size(),b=integer(o,"permutations",999,1,100000),seed=integer(o,"seed",1,0,Integer.MAX_VALUE);
        CellData.checkOutputSize((long)sectionEdges.size()*k*(k+1)/2);
        if((long)(g.original.size()+edges.size()+sections(g).size()*k*k)*b>100_000_000)throw new IllegalArgumentException("Neighborhood permutations exceed 100 million visits");
        List<String[]> rows=new ArrayList<>();List<Double> pv=new ArrayList<>();Random random=new Random(seed);
        for(var section:sections(g).entrySet()) {
            List<Integer> nodes=section.getValue();Map<Integer,Integer> index=new HashMap<>();for(int i=0;i<nodes.size();i++)index.put(nodes.get(i),i);
            List<Edge> local=sectionEdges.get(section.getKey()).stream().map(e->new Edge(index.get(e.source()),index.get(e.target()),e.distance())).toList();
            double[] labels=nodes.stream().mapToDouble(i->types.indexOf(data.cell(g.original.get(i),"cell_type"))).toArray();
            var strata=SpatialStatistics.groups(nodes.stream().map(i->g.strata[i]).toArray(String[]::new));
            int[] observed=pairCounts(labels,local,k),up=new int[k*k],down=new int[k*k];double[] mean=new double[k*k];
            double[] shuffled=new double[labels.length];
            for(int iteration=0;iteration<b;iteration++) {
                System.arraycopy(labels,0,shuffled,0,labels.length);SpatialStatistics.shuffle(shuffled,strata,random);int[] counts=pairCounts(shuffled,local,k);
                for(int i=0;i<counts.length;i++){mean[i]+=(double)counts[i]/b;if(counts[i]>=observed[i])up[i]++;if(counts[i]<=observed[i])down[i]++;}
            }
            for(int a=0;a<k;a++)for(int c=a;c<k;c++) {
                int j=a*k+c;double p=(up[j]+1.0)/(b+1),d=(down[j]+1.0)/(b+1);
                rows.add(new String[]{section.getKey(),types.get(a),types.get(c),Integer.toString(observed[j]),Double.toString(mean[j]),local.isEmpty()?"NA":Double.toString(p),local.isEmpty()?"NA":Double.toString(d),"NA","NA",local.isEmpty()?"no_edges":"tested"});pv.add(local.isEmpty()?1:p);pv.add(local.isEmpty()?1:d);
            }
        }
        double[] q=NetworkAnalysis.bh(pv.stream().mapToDouble(Double::doubleValue).toArray());for(int i=0;i<rows.size();i++)if(rows.get(i)[9].equals("tested")){rows.get(i)[7]=Double.toString(q[2*i]);rows.get(i)[8]=Double.toString(q[2*i+1]);}
        table(dir.resolve("neighborhoods.tsv"),List.of("section_id","type_a","type_b","edges","null_mean","enrichment_p","depletion_p","enrichment_bh","depletion_bh","status"),rows);
        manifest.put("null","Cell-type labels permuted within section/compartment/stratum, preserving label counts and graph geometry. BH covers all sections, unordered label pairs and both tails; expression/proximity alone do not identify signaling.");manifest.put("seed",seed);manifest.put("permutations",b);
    }
    private static void compare(CellData data,Geometry g,List<Edge> edges,Map<String,String> o,Path dir,Map<String,Object> manifest) throws IOException {
        List<String> types=types(data,g);int k=types.size();Map<String,int[]> counts=new LinkedHashMap<>();Map<String,Integer> totals=new LinkedHashMap<>();
        for(int node:g.original){String sample=data.cell(node,"sample_id");counts.putIfAbsent(sample,new int[k*k]);totals.putIfAbsent(sample,0);}
        for(Edge edge:edges) {
            int s=g.original.get(edge.source()),t=g.original.get(edge.target());String sample=data.cell(s,"sample_id");int a=types.indexOf(data.cell(s,"cell_type")),b=types.indexOf(data.cell(t,"cell_type"));
            counts.get(sample)[Math.min(a,b)*k+Math.max(a,b)]++;totals.merge(sample,1,Integer::sum);
        }
        if(totals.values().stream().anyMatch(n->n==0))throw new IllegalArgumentException("Every analyzed sample must have edges; zero-edge samples cannot supply a spatial fraction");
        List<String> ids=new ArrayList<>(counts.keySet());SampleInference.Design design=CellDesign.create(data,ids,o);CellDesign.write(design,ids,dir.resolve("design.tsv"));
        CellData.checkOutputSize((long)ids.size()*k*(k+1)/2);
        CellDesign.fitWork(design,(long)k*(k+1)/2);
        List<String[]> rows=new ArrayList<>(),results=new ArrayList<>();
        for(int a=0;a<k;a++)for(int b=a;b<k;b++) {
            double[] y=new double[ids.size()];
            for(int i=0;i<ids.size();i++) {
                int n=totals.get(ids.get(i)),v=counts.get(ids.get(i))[a*k+b];double f=(double)v/n;y[i]=Math.asin(Math.sqrt(f));
                rows.add(new String[]{ids.get(i),types.get(a),types.get(b),Integer.toString(v),Integer.toString(n),Double.toString(f)});
            }
            results.add(SingleCellCli.fitRow(types.get(a),types.get(b),y,design));
        }
        table(dir.resolve("sample-relationships.tsv"),List.of("sample_id","type_a","type_b","pair_edges","all_edges","fraction"),rows);SingleCellCli.adjustWrite(dir.resolve("results.tsv"),results);
        manifest.put("estimand","Equal-sample Gaussian comparison of arcsine-square-root edge-pair fractions. All sections within a sample pooled before inference. Relative edge composition depends on abundance and geometry; not abundance-adjusted interaction or causality.");
    }
    private static void gradient(CellData data,Geometry g,Map<String,String> o,Path dir,Map<String,Object> manifest) throws IOException {
        String column=required(o,"distance-column");List<Integer> features=features(data,o);Map<String,List<Integer>> groups=new LinkedHashMap<>();
        double[] distances=new double[data.ids.size()];
        for(int i:g.original){double d=finite(data.cell(i,column),"distance");if(d<0)throw new IllegalArgumentException("Distance must be nonnegative physical units");distances[i]=d;groups.computeIfAbsent(data.cell(i,"sample_id"),s->new ArrayList<>()).add(i);}
        if((long)g.original.size()*features.size()>20_000_000)throw new IllegalArgumentException("Gradient exceeds 20 million observation-feature visits");
        List<String> ids=new ArrayList<>(groups.keySet());SampleInference.Design design=CellDesign.create(data,ids,o);CellDesign.write(design,ids,dir.resolve("design.tsv"));
        CellData.checkOutputSize((long)ids.size()*features.size());
        CellDesign.fitWork(design,features.size());
        List<String[]> rows=new ArrayList<>(),results=new ArrayList<>();
        for(int j:features) {
            double[] slopes=new double[ids.size()];
            for(int s=0;s<ids.size();s++) {
                var nodes=groups.get(ids.get(s));if(nodes.size()<3)throw new IllegalArgumentException("Gradient needs >=3 observations/sample");
                double[] x=new double[nodes.size()],y=new double[nodes.size()];
                for(int i=0;i<nodes.size();i++){x[i]=distances[nodes.get(i)];y[i]=data.normalized(nodes.get(i),j);}
                slopes[s]=SpatialStatistics.linearSlope(x,y);
                rows.add(new String[]{ids.get(s),data.genes.get(j),Integer.toString(nodes.size()),Double.toString(slopes[s])});
            }
            results.add(SingleCellCli.fitRow("distance_slope",data.genes.get(j),slopes,design));
        }
        table(dir.resolve("sample-slopes.tsv"),List.of("sample_id","feature_id","observations","slope_per_physical_unit"),rows);SingleCellCli.adjustWrite(dir.resolve("results.tsv"),results);
        manifest.put("estimand","Tested minus reference mean sample-specific linear log1p-normalized expression/distance slope. All section observations pooled within sample; spatially correlated observations only define descriptive slopes, with no cell-level p values. Gaussian independent-sample slope errors assumed. Supplied boundaries treated as fixed; no nonlinear/interface uncertainty model.");
    }
    private static void writeGraph(CellData data,Geometry g,List<Edge> edges,Map<String,List<Edge>> sectionEdges,Path dir,String units) throws IOException {
        List<String[]> rows=new ArrayList<>();int[] degrees=new int[g.original.size()];
        for(Edge e:edges){degrees[e.source()]++;degrees[e.target()]++;rows.add(new String[]{data.ids.get(g.original.get(e.source())),data.ids.get(g.original.get(e.target())),Double.toString(e.distance()),g.section[e.source()]});}
        table(dir.resolve("edges.tsv"),List.of("source","target","distance_"+units,"section_id"),rows);rows=new ArrayList<>();
        for(int i=0;i<degrees.length;i++)rows.add(new String[]{data.ids.get(g.original.get(i)),data.cell(g.original.get(i),"sample_id"),g.section[i],Integer.toString(degrees[i]),degrees[i]==0?"isolated":"connected"});
        table(dir.resolve("nodes.tsv"),List.of("obs_id","sample_id","section_id","degree","status"),rows);
        // Separate SVG per section: no artificial overlay of unrelated coordinate frames.
        List<String[]> maps=new ArrayList<>();int count=0;
        for(var section:sections(g).entrySet()) {
            String name="section-"+(++count)+".svg";maps.add(new String[]{section.getKey(),name});var nodes=section.getValue();
            double xmin=Double.POSITIVE_INFINITY,xmax=Double.NEGATIVE_INFINITY,ymin=xmin,ymax=xmax;
            for(int i:nodes){xmin=Math.min(xmin,g.xyz[i][0]);xmax=Math.max(xmax,g.xyz[i][0]);ymin=Math.min(ymin,g.xyz[i][1]);ymax=Math.max(ymax,g.xyz[i][1]);}
            if(!Double.isFinite(xmax-xmin)||!Double.isFinite(ymax-ymin))throw new IllegalArgumentException("Coordinate extent overflows; rescale physical units");
            double scale=500/Math.max(1,Math.max(xmax-xmin,ymax-ymin));
            try(var out=Files.newBufferedWriter(dir.resolve(name))) {
                out.write("<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 560 580\"><rect width=\"560\" height=\"580\" fill=\"white\"/><text x=\"20\" y=\"20\">"+xml(section.getKey())+" (XY projection)</text>");
                for(Edge e:sectionEdges.get(section.getKey()))out.write("<line stroke=\"#ccc\" x1=\""+(30+(g.xyz[e.source()][0]-xmin)*scale)+"\" y1=\""+(40+(g.xyz[e.source()][1]-ymin)*scale)+"\" x2=\""+(30+(g.xyz[e.target()][0]-xmin)*scale)+"\" y2=\""+(40+(g.xyz[e.target()][1]-ymin)*scale)+"\"/>");
                for(int i:nodes)out.write("<circle fill=\""+(degrees[i]==0?"#d94a38":"#236b8e")+"\" r=\"3\" cx=\""+(30+(g.xyz[i][0]-xmin)*scale)+"\" cy=\""+(40+(g.xyz[i][1]-ymin)*scale)+"\"><title>"+xml(data.ids.get(g.original.get(i)))+"; "+xml(data.cell(g.original.get(i),"cell_type"))+"</title></circle>");
                out.write("<text x=\"20\" y=\"565\">x: "+xmin+" to "+xmax+"; y: "+ymin+" to "+ymax+" "+units+"</text></svg>");
            }
        }
        table(dir.resolve("maps.tsv"),List.of("section_id","file"),maps);
    }
    private static String xml(String v){return v.replace("&","&amp;").replace("<","&lt;").replace(">","&gt;").replace("\"","&quot;");}
}
