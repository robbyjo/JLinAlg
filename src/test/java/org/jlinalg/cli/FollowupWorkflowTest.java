/* Copyright (C) 2026 JLinAlg contributors; SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import com.sun.net.httpserver.HttpServer;
import java.net.InetSocketAddress;
import org.jlinalg.network.NetworkAnalysis;

class FollowupWorkflowTest {
    @TempDir Path dir;
    Path file(String name,String text)throws Exception {Path p=dir.resolve(name);Files.createDirectories(p.getParent());Files.writeString(p,text);return p;}
    String output;
    int run(String... args){var b=new ByteArrayOutputStream();var s=new PrintStream(b);List<String> all=new ArrayList<>(List.of(args));all.add("--no-config");int status=JLinAlgCli.run(all.toArray(String[]::new),s,s);output=b.toString();return status;}
    static final String VARIANTS="genome_build\tchrom\tpos\tref\talt\nGRCh38\t1\t101\tA\tG\nGRCh38\t1\t202\tC\tT\n";
    @Test void projectOverridesLocalListsReplaceAndCliWins()throws Exception {
        file("home/.jlinalg/config.yaml","schema_version: 1\ncommands:\n  network:\n    method: sparse\n    lambda: 0.2\n    matrix: {path: local.tsv}\n    regulators: [A, B]\n");
        file("project/jlinalg.yaml","schema_version: 1\ncommands:\n  network:\n    lambda: 0.1\n    regulators: [C]\n    matrix: {path: cohort.tsv}\n");
        var resolved=ProjectConfiguration.resolve(new String[]{"network","--lambda","0.3"},dir.resolve("project"),dir.resolve("home"));
        List<String> a=Arrays.asList(resolved.arguments());assertEquals(1,Collections.frequency(a,"--lambda"));assertEquals("0.3",a.get(a.indexOf("--lambda")+1));
        assertEquals("C",a.get(a.indexOf("--regulators")+1));assertFalse(a.contains("A"));assertEquals(dir.resolve("project/cohort.tsv").toString(),a.get(a.indexOf("--matrix")+1));
        assertEquals(2,((List<?>)resolved.provenance().get("sources")).size());
        assertArrayEquals(new String[]{"network","--help"},ProjectConfiguration.resolve(new String[]{"network","--no-config","--help"},dir.resolve("project"),dir.resolve("home")).arguments());
    }
    @Test void unsafeDuplicateAndUnknownYamlFail()throws Exception {
        Path config=file("bad.yaml","schema_version: 1\ncommands: !!java.net.URL [https://example.org]\n");
        assertThrows(RuntimeException.class,()->ProjectConfiguration.resolve(new String[]{"network","--config",config.toString()},dir,dir));
        Files.writeString(config,"schema_version: 1\nschema_version: 1\n");assertThrows(RuntimeException.class,()->ProjectConfiguration.resolve(new String[]{"network","--config",config.toString()},dir,dir));
        Files.writeString(config,"schema_version: 1\ncommnads: {}\n");assertThrows(RuntimeException.class,()->ProjectConfiguration.resolve(new String[]{"network","--config",config.toString()},dir,dir));
    }
    @Test void localAnnotationPreservesTranscriptsMissingAndRejectsCorruption()throws Exception {
        Path input=file("variants.tsv",VARIANTS),source=file("source.tsv","genome_build\tchrom\tpos\tref\talt\tgene\ttranscript\nGRCh38\t1\t101\tA\tG\tGENE1\tTX1\nGRCh38\t1\t101\tA\tG\tGENE1\tTX2\n");
        Path db=dir.resolve("db"),out=dir.resolve("annotated");
        assertEquals(0,run("variant-db","--source-file",source.toString(),"--genome-build","GRCh38","--release","fixture-1","--license","CC0","--out",db.toString()),output);
        assertEquals(0,run("variant-annotate","--input",input.toString(),"--genome-build","GRCh38","--database",db.toString(),"--out",out.toString()),output);
        var table=DelimitedData.read(out.resolve("annotations.tsv"));assertEquals(3,table.rows().size());assertEquals("not_found",table.rows().get(2)[1]);assertTrue(Files.exists(Path.of(out+".config.yaml")));
        Files.writeString(db.resolve("annotations.tsv"),"corrupt");assertEquals(2,run("variant-annotate","--input",input.toString(),"--genome-build","GRCh38","--database",db.toString(),"--out",dir.resolve("broken").toString()));assertFalse(Files.exists(dir.resolve("broken")));
    }
    @Test void postAndVepRestUseValidatedResponseAndAtomicFailure()throws Exception {
        HttpServer server=HttpServer.create(new InetSocketAddress("127.0.0.1",0),0);
        server.createContext("/lookup",exchange->{byte[] body=exchange.getRequestBody().readAllBytes();assertEquals("POST",exchange.getRequestMethod());exchange.sendResponseHeaders(200,body.length);exchange.getResponseBody().write(body);exchange.close();});
        server.createContext("/vep",exchange->{exchange.getRequestBody().readAllBytes();byte[] body=("[{\"input\":\"1 101 GRCh38:1:101:A:G A G . . .\",\"assembly_name\":\"GRCh38\",\"transcript_consequences\":[{\"gene_id\":\"G\",\"transcript_id\":\"T\",\"consequence_terms\":[\"missense_variant\"]}]},{\"input\":\"1 202 GRCh38:1:202:C:T C T . . .\",\"assembly_name\":\"GRCh38\",\"most_severe_consequence\":\"intergenic_variant\"}]").getBytes();exchange.sendResponseHeaders(200,body.length);exchange.getResponseBody().write(body);exchange.close();});server.start();
        try {
            Path input=file("variants.tsv",VARIANTS);String base="http://127.0.0.1:"+server.getAddress().getPort();
            assertEquals(0,run("variant-annotate","--input",input.toString(),"--genome-build","GRCh38","--backend","post","--endpoint",base+"/lookup","--out",dir.resolve("post").toString()),output);
            assertEquals(0,run("variant-consequence","--input",input.toString(),"--genome-build","GRCh38","--engine","vep-rest","--endpoint",base+"/vep","--out",dir.resolve("rest").toString()),output);
            assertTrue(Files.readString(dir.resolve("rest/consequences.tsv")).contains("missense_variant"));
        } finally {server.stop(0);}
    }
    @Test void vepAndAnnovarImportsPreserveAllelesAndGeneContext()throws Exception {
        Path input=file("variants.tsv",VARIANTS),vep=file("vep.tsv","##VEP=fixture\n#Uploaded_variation\tGene\tFeature\tConsequence\tExtra\nGRCh38:1:101:A:G\tG\tTX1\tmissense_variant\tIMPACT=MODERATE\nGRCh38:1:101:A:G\tG\tTX2\tintron_variant\tIMPACT=MODIFIER\n");
        assertEquals(0,run("variant-consequence","--input",input.toString(),"--genome-build","GRCh38","--engine","import","--annotations",vep.toString(),"--format","vep","--out",dir.resolve("vep").toString()),output);
        assertEquals(3,DelimitedData.read(dir.resolve("vep/consequences.tsv")).rows().size());
        Path anno=file("anno.vcf","##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t101\tGRCh38:1:101:A:G\tA\tG\t.\t.\tGene.refGene=G;Func.refGene=exonic;ExonicFunc.refGene=nonsynonymous_SNV\n");
        assertEquals(0,run("variant-consequence","--input",input.toString(),"--genome-build","GRCh38","--engine","import","--annotations",anno.toString(),"--format","annovar-vcf","--out",dir.resolve("annovar").toString()),output);
        assertTrue(Files.readString(dir.resolve("annovar/consequences.tsv")).contains("nonsynonymous_SNV"));
        Files.writeString(anno,Files.readString(anno).replace("\tA\tG\t","\tA\tT\t"));
        assertEquals(2,run("variant-consequence","--input",input.toString(),"--genome-build","GRCh38","--engine","import","--annotations",anno.toString(),"--format","annovar-vcf","--out",dir.resolve("mismatch").toString()));
    }
    @Test void scoreRetainsComponentsAndDoesNotReweightMissingEvidence()throws Exception {
        Path input=file("evidence.tsv","variant_id\tpip\tfrequency\na\t0.8\t0.1\nb\t0.9\tNA\n"),weights=file("weights.tsv","column\tweight\tminimum\tmaximum\tdirection\npip\t3\t0\t1\thigher\nfrequency\t1\t0\t1\tlower\n");
        assertEquals(0,run("variant-score","--input",input.toString(),"--components",weights.toString(),"--out",dir.resolve("score").toString()),output);
        var t=DelimitedData.read(dir.resolve("score/scores.tsv"));assertEquals(.825,Double.parseDouble(t.rows().get(0)[t.column("priority_score")]),1e-14);assertEquals("NA",t.rows().get(1)[t.column("priority_score")]);assertEquals("0.75",t.rows().get(1)[t.column("weight_coverage")]);
        Files.writeString(weights,"column\tweight\tminimum\tmaximum\tdirection\npip\t1\t0\t1\thigher\nfrequency\t1e-30\t0\t1\tlower\n");
        assertEquals(0,run("variant-score","--input",input.toString(),"--components",weights.toString(),"--out",dir.resolve("tiny-weight").toString()),output);
        var tiny=DelimitedData.read(dir.resolve("tiny-weight/scores.tsv"));assertEquals("missing_evidence",tiny.rows().get(1)[tiny.column("score_status")]);
    }
    @Test void referencePropagationAndPermutationAreReproducible()throws Exception {
        Path edges=file("edges.tsv","source\ttarget\tweight\ttype\nA\tB\t1\tphysical\nB\tC\t1\tphysical\nC\tD\t1\tphysical\n"),hits=file("hits.tsv","gene\nA\nUNKNOWN\n");
        for(String name:List.of("one","two"))assertEquals(0,run("network","--method","reference","--edges",edges.toString(),"--hits",hits.toString(),"--permutations","99","--out",dir.resolve(name).toString()),output);
        assertEquals(Files.readString(dir.resolve("one/nodes.tsv")),Files.readString(dir.resolve("two/nodes.tsv")));assertTrue(Files.readString(dir.resolve("one/unmapped.tsv")).contains("UNKNOWN"));
        var t=DelimitedData.read(dir.resolve("one/nodes.tsv"));double mass=t.rows().stream().mapToDouble(r->Double.parseDouble(r[t.column("propagation_score")])).sum();assertEquals(1,mass,1e-9);
    }
    @Test void numericNetworkBoundariesAndAnalyticLasso() {
        double[][] x={{-1,-1},{-1,-1},{1,1},{1,1}};
        double[][] coefficients=NetworkAnalysis.neighborhoodLasso(x,.2,10000,1e-10);assertEquals(.8,coefficients[0][1],1e-8);assertEquals(.8,coefficients[1][0],1e-8);
        assertThrows(IllegalArgumentException.class,()->NetworkAnalysis.correlation(new double[][]{{1,1},{1,2},{1,3},{1,4}}));
        assertArrayEquals(new double[]{.03,.04,.04},NetworkAnalysis.bh(new double[]{.01,.04,.03}),1e-14);
        var graph=List.of(Map.of(1,1.0),Map.<Integer,Double>of());double[] score=NetworkAnalysis.propagate(graph,new double[]{1,0},.5,1e-12,10000);assertArrayEquals(new double[]{2.0/3,1.0/3},score,1e-10);
    }
    @Test void differentialAlignsFeaturesRejectsOverlapAndMarksDegeneracy()throws Exception {
        Path a=file("a.tsv","sample_id\tA\tB\na1\t1\t1\na2\t2\t2\na3\t3\t3\na4\t4\t4\na5\t5\t5\n"),b=file("b.tsv","sample_id\tB\tA\nb1\t1\t2\nb2\t3\t3\nb3\t2\t1\nb4\t5\t4\nb5\t4\t5\n");
        assertEquals(0,run("network","--method","differential","--matrix",a.toString(),"--matrix-b",b.toString(),"--out",dir.resolve("diff").toString()),output);assertTrue(Files.readString(dir.resolve("diff/edges.tsv")).contains("degenerate_correlation"));
        assertEquals(2,run("network","--method","differential","--matrix",a.toString(),"--matrix-b",a.toString(),"--out",dir.resolve("paired").toString()));assertTrue(output.contains("overlapping"));
    }
    @Test void missingAdapterExecutableFailsWithoutPublishingPartialResults()throws Exception {
        Path edges=file("signs.tsv","source\ttarget\tsign\nA\tB\t1\n"),acts=file("acts.tsv","gene\tactivity\nB\t1\n");
        assertEquals(2,run("network","--method","signaling","--edges",edges.toString(),"--activities",acts.toString(),"--rscript",dir.resolve("absent.exe").toString(),"--out",dir.resolve("failed").toString()));assertFalse(Files.exists(dir.resolve("failed")));
    }
    @Test void nativeNetworksMatchIndependentRAndGlmnetFixtures()throws Exception {
        Path root=Path.of("src/test/resources/network-reference");
        var a=DelimitedData.read(root.resolve("a.tsv"));int n=a.rows().size(),p=a.header().size()-1;double[][] x=new double[n][p];
        for(int i=0;i<n;i++)for(int j=0;j<p;j++)x[i][j]=Double.parseDouble(a.rows().get(i)[j+1]);
        double[][] actual=NetworkAnalysis.neighborhoodLasso(x,.15,100000,1e-12);var expected=DelimitedData.read(root.resolve("lasso.tsv"));
        for(int i=0;i<p;i++)for(int j=0;j<p;j++)assertEquals(Double.parseDouble(expected.rows().get(i)[j+1]),actual[i][j],1e-8,"coefficient "+i+","+j);
        assertEquals(0,run("network","--method","differential","--matrix",root.resolve("a.tsv").toString(),"--matrix-b",root.resolve("b.tsv").toString(),"--out",dir.resolve("reference-diff").toString()),output);
        var results=DelimitedData.read(dir.resolve("reference-diff/edges.tsv"));var reference=DelimitedData.read(root.resolve("differential.tsv"));Map<String,String[]> map=new HashMap<>();for(String[] row:reference.rows())map.put(row[0]+":"+row[1],row);
        for(String[] row:results.rows()) {String[] ref=map.get(row[0]+":"+row[1]);for(String column:List.of("r_a","r_b","z","p","bh"))assertEquals(Double.parseDouble(ref[reference.column(column)]),Double.parseDouble(row[results.column(column)]),1e-11,column);}
    }
    @Test void sourceSchemaMappingsAndVcfSnapshotsAreSearchable()throws Exception {
        Path source=file("scores.tsv","#chr\tposition\treference\talternate\tCADD\n1\t101\tA\tG\t25\n"),input=file("queries.tsv",VARIANTS);
        assertEquals(0,run("variant-db","--source-file",source.toString(),"--format","table","--chrom-column","#chr","--pos-column","position","--ref-column","reference","--alt-column","alternate","--genome-build","GRCh38","--release","test","--license","test","--out",dir.resolve("mapped").toString()),output);
        assertEquals(0,run("variant-annotate","--input",input.toString(),"--genome-build","GRCh38","--database",dir.resolve("mapped").toString(),"--out",dir.resolve("query").toString()),output);
        assertTrue(Files.readString(dir.resolve("query/annotations.tsv")).contains("annotation.source.CADD"));
        Path vcf=file("source.vcf","##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t101\trs1\tA\tG\t.\t.\tCLNSIG=Pathogenic\n");
        assertEquals(0,run("variant-db","--source-file",vcf.toString(),"--format","vcf","--genome-build","GRCh38","--release","test","--license","test","--out",dir.resolve("vcf-db").toString()),output);
        assertTrue(Files.readString(dir.resolve("vcf-db/annotations.tsv")).contains("CLNSIG=Pathogenic"));
    }
    @Test void graphmlIncludesIsolatedFeaturesAndEscapesIdentifiers()throws Exception {
        Path matrix=file("matrix.tsv","sample_id\tA&B\tC\na\t-1\t-1\nb\t-1\t1\nc\t1\t-1\nd\t1\t1\n");
        assertEquals(0,run("network","--method","sparse","--matrix",matrix.toString(),"--lambda","0.2","--out",dir.resolve("graph").toString()),output);
        var factory=javax.xml.parsers.DocumentBuilderFactory.newInstance();factory.setFeature("http://apache.org/xml/features/disallow-doctype-decl",true);
        var doc=factory.newDocumentBuilder().parse(dir.resolve("graph/network.graphml").toFile());
        assertEquals(2,doc.getElementsByTagName("node").getLength());assertEquals(0,doc.getElementsByTagName("edge").getLength());assertTrue(Files.readString(dir.resolve("graph/network.graphml")).contains("A&amp;B"));
    }
}
