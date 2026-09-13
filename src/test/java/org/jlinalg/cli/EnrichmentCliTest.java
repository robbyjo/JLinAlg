/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class EnrichmentCliTest {
    @TempDir Path dir;
    Path file(String name,String text) throws Exception { Path p=dir.resolve(name); Files.writeString(p,text); return p; }
    @Test void selectionGrammarIsSafeAndExplicit() {
        var expression=new EnrichmentExpression("status == 'OK' && p < bonferroni(0.05) && abs(beta) >= 2",Set.of("status","p","beta"));
        assertTrue(expression.test(Map.of("status","OK","p","0.001","beta","-2"),20));
        assertFalse(expression.test(Map.of("status","OK","p","0.01","beta","-2"),20));
        assertNull(expression.test(Map.of("status","OK","p","NA","beta","-2"),20));
        assertThrows(IllegalArgumentException.class,()->new EnrichmentExpression("unknown < 1",Set.of()));
        assertThrows(IllegalArgumentException.class,()->new EnrichmentExpression("System.exit(0)",Set.of()));
        assertThrows(IllegalArgumentException.class,()->new EnrichmentExpression("p < bonferroni(.05,0)",Set.of("p")).test(Map.of("p",".1"),10));
        assertTrue(new EnrichmentExpression("is_missing(p) || `effect size` > 0",Set.of("p","effect size")).test(Map.of("p","NA","effect size","2"),10));
        assertFalse(new EnrichmentExpression("is_finite(p) && p < .05",Set.of("p")).test(Map.of("p","Infinity"),10));
        assertNull(new EnrichmentExpression("!(p < .05)",Set.of("p")).test(Map.of("p","NA"),10));
    }
    @Test void analysisDerivedUniverseAndPreMappingBonferroni() throws Exception {
        Path input=file("input.csv","probe,p,status\np1,0.009,OK\np2,0.02,OK\np3,0.5,OK\np4,0.5,OK\np5,0.5,OK\nbad,NA,FAIL\n");
        Path annotation=file("annot.csv","IlmnID,gene,chr\np1,A;B,1\np2,A,1\np3,C,2\np4,D,3\n");
        Path sets=file("sets.gmt","T\tTerm\tA\tB\nZ\tZero\tD\n");
        Path out=dir.resolve("out.tsv");
        ByteArrayOutputStream bytes=new ByteArrayOutputStream(); PrintStream console=new PrintStream(bytes);
        String[] args={"--enrichment","Custom","--enrichment-db",sets.toString(),"--input",input.toString(),"--input-id","probe",
            "--universe-selection","status == 'OK'","--selection","p < bonferroni(.05)","--annot",annotation.toString(),
            "--annot-id","IlmnID","--gene-col","gene","--annot-cols","chr","--min-set-size","1","--out",out.toString()};
        assertEquals(0,JLinAlgCli.run(args,console,console),bytes.toString());
        assertTrue(bytes.toString().contains("bonferroni_default_m=5"));
        assertTrue(bytes.toString().contains("selected_features=1"));
        assertTrue(bytes.toString().contains("unmapped_features=1"));
        var results=DelimitedData.read(out);
        assertEquals(2,results.rows().size());
        assertEquals("2",results.rows().get(0)[results.column("overlap_genes")]);
        assertTrue(Files.readString(Path.of(out+".features.tsv")).contains("annotation.chr"));
        String original=Files.readString(out);
        assertEquals(2,JLinAlgCli.run(args,console,console)); assertEquals(original,Files.readString(out));
    }
    @Test void missingValuesRequireExplicitExclusion() throws Exception {
        Path input=file("input.tsv","id\tp\nA\t0.001\nB\tNA\nC\t0.5\n");
        Path sets=file("sets.gmt","T\tTerm\tA\n");
        List<String> args=new ArrayList<>(List.of("--enrichment","Custom","--enrichment-db",sets.toString(),"--input",input.toString(),"--input-id","id","--selection","p < .05","--min-set-size","1","--out",dir.resolve("out.tsv").toString()));
        ByteArrayOutputStream b=new ByteArrayOutputStream(); PrintStream s=new PrintStream(b);
        assertEquals(2,JLinAlgCli.run(args.toArray(String[]::new),s,s)); assertTrue(b.toString().contains("unevaluable selection"));
        args.set(args.size()-1,dir.resolve("good.tsv").toString()); args.addAll(List.of("--missing-selection","exclude"));
        assertEquals(0,JLinAlgCli.run(args.toArray(String[]::new),s,s),b.toString());
        assertTrue(b.toString().contains("unevaluable_excluded=1"));
    }
    @Test void catalogUnionRetainsAllEligibleTraitsAndCleansMappings() throws Exception {
        Path catalog=file("catalog.tsv","P-VALUE\tSNPS\tDISEASE/TRAIT\tMAPPED_TRAIT_URI\tREPORTED GENE(S)\tMAPPED_GENE\n"
            +"1e-9\trs1\tTrait1\thttp://www.ebi.ac.uk/efo/EFO_0001\tA\tB, C - D\n"
            +"1e-9\trs2\tTrait2\thttp://www.ebi.ac.uk/efo/EFO_0002\tE\tNR\n"
            +"1000000000\trs3\tBad\thttp://www.ebi.ac.uk/efo/EFO_0003\tA\tA\n"
            +"1e-6\tchr11:102751102\"-?\tLiteral quote\thttp://www.ebi.ac.uk/efo/EFO_0003\tA\tA\n");
        var db=EnrichmentDatabase.load(catalog,"GWASCatalog",Map.of());
        assertEquals(Set.of("A","B","C","D"),db.genes.get("EFO:0001"));
        assertEquals(Set.of("E"),db.genes.get("EFO:0002")); assertEquals(1,db.rejectedRows);
        assertEquals("union",db.metadata.getProperty("catalog_gene_source"));
        assertEquals(Set.of("A"),EnrichmentDatabase.load(catalog,"GWASCatalog",Map.of("--catalog-gene-source","reported")).genes.get("EFO:0001"));
    }
    @Test void ontologyUsesSafeEdgesAndRejectsCycles() throws Exception {
        Path obo=file("go.obo","format-version: 1.2\n[Term]\nid: GO:1\nname: Parent\nnamespace: biological_process\n[Term]\nid: GO:2\nname: Child\nnamespace: biological_process\nalt_id: GO:20\nis_a: GO:1 ! Parent\n[Term]\nid: GO:3\nname: Regulator\nnamespace: biological_process\nrelationship: regulates GO:1\n");
        Path gaf=file("go.gaf","DB\tP1\tA\t\tGO:20\tPMID:1\tIDA\t\tP\t\t\tprotein\ttaxon:9606\t20260913\tDB\n"
            +"DB\tP2\tB\tNOT\tGO:1\tPMID:1\tIDA\t\tP\t\t\tprotein\ttaxon:9606\t20260913\tDB\n"
            +"DB\tP3\tC\t\tGO:3\tPMID:1\tIDA\t\tP\t\t\tprotein\ttaxon:9606\t20260913\tDB\n");
        var db=EnrichmentDatabase.load(gaf,"GO:BP",Map.of("--ontology",obo.toString()));
        assertEquals(Set.of("A"),db.genes.get("GO:1")); assertEquals(Set.of("A"),db.genes.get("GO:2"));
        Files.writeString(obo,Files.readString(obo)+"is_a: GO:3\n");
        assertThrows(java.io.IOException.class,()->EnrichmentDatabase.load(gaf,"GO",Map.of("--ontology",obo.toString())));
    }
    @Test void localDownloadPackageIntegrityAndNamespace() throws Exception {
        Path gmt=file("source.gmt","T\tTerm\tA\tB\n"); Path target=dir.resolve("database");
        var options=Map.of("--source-file",gmt.toString(),"--gene-id-type","symbol");
        new EnrichmentDownloads(new PrintStream(new ByteArrayOutputStream()),100000).install("MSigDB:H",target,options);
        assertEquals(1,EnrichmentDatabase.load(target,"MSigDB:H",Map.of()).genes.size());
        assertThrows(java.io.IOException.class,()->EnrichmentDatabase.load(target,"MSigDB:GO:BP",Map.of()));
        assertThrows(java.io.IOException.class,()->EnrichmentDatabase.load(target,"MSigDB:H",Map.of("--gene-id-type","ensembl")));
        Files.writeString(target.resolve("data.gmt"),"corrupted");
        assertThrows(java.io.IOException.class,()->EnrichmentDatabase.load(target,"MSigDB:H",Map.of()));
    }
    @Test void gsamethCliUsesTheFullProbeUniverse() throws Exception {
        Path input=file("ewas.tsv","probe\tp\np1\t.001\np2\t.2\np3\t.5\np4\t.001\np5\t.5\np6\t.5\np7\t.5\np8\t.5\n");
        Path annotation=file("annot.tsv","probe\tgenes\np1\tA;B\np2\tA\np3\tC\np4\tD\np5\tE\np6\tF\np7\tG\np8\tH\n");
        Path sets=file("sets.gmt","T\tTerm\tA\tB\tD\nZ\tZero\tG\tH\n");
        ByteArrayOutputStream b=new ByteArrayOutputStream(); PrintStream s=new PrintStream(b);
        assertEquals(0,JLinAlgCli.run(new String[]{"--enrichment","Custom","--enrichment-method","gsameth","--enrichment-db",sets.toString(),"--input",input.toString(),"--input-id","probe","--selection","p < .05","--annot",annotation.toString(),"--annot-id","probe","--gene-col","genes","--min-set-size","1","--out",dir.resolve("ewas-out.tsv").toString()},s,s),b.toString());
        assertTrue(b.toString().contains("fractional weights"));
    }
    @Test void sourceSchemasSpeciesAndParquet() throws Exception {
        Path csv=file("membership.csv","term_id,term_name,gene_id\nT,\"Quoted, term\",A\n");
        Path imported=dir.resolve("csv-package");
        new EnrichmentDownloads(new PrintStream(new ByteArrayOutputStream()),100000).install("Custom",imported,Map.of("--source-file",csv.toString(),"--gene-id-type","symbol"));
        assertEquals(Set.of("A"),EnrichmentDatabase.load(imported,"Custom",Map.of()).genes.get("T"));
        Path reactome=file("reactome.tsv","ENSG1\tR-HSA-1\turl\tHuman pathway\tIEA\tHomo sapiens\nENSM1\tR-MMU-1\turl\tMouse pathway\tIEA\tMus musculus\n");
        assertEquals(Set.of("ENSG1"),EnrichmentDatabase.load(reactome,"Reactome",Map.of()).genes.get("R-HSA-1"));
        Path hpo=file("hpo.tsv","ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\n1\tA\tHP:1\tPhenotype\n");
        assertEquals(Set.of("1"),EnrichmentDatabase.load(hpo,"HPO",Map.of()).genes.get("HP:1"));
        Path monarch=file("monarch.tsv","subject\tobject\tobject_label\tsubject_taxon\nNCBIGene:1\tMONDO:1\tDisease\tNCBITaxon:9606\nMGI:2\tMONDO:1\tDisease\tNCBITaxon:10090\n");
        assertEquals(Set.of("NCBIGene:1"),EnrichmentDatabase.load(monarch,"Monarch:Mondo",Map.of()).genes.get("MONDO:1"));
        Path kegg=file("links.tsv","hsa:1\tpath:hsa00010\n");
        assertEquals(Set.of("1"),EnrichmentDatabase.load(kegg,"KEGG",Map.of("--db-format","kegg")).genes.get("hsa00010"));
        Path target=Files.createDirectory(dir.resolve("ot"));
        Path associations=Files.createDirectory(target.resolve("associations"));
        Path parquet=associations.resolve("part.parquet");
        try(var connection=java.sql.DriverManager.getConnection("jdbc:duckdb:"); var query=connection.createStatement()) {
            query.execute("COPY (SELECT * FROM (VALUES ('ENSG1','EFO:1',0.9),('ENSG2','EFO:1',0.2)) AS t(targetId,diseaseId,score)) TO '"+parquet.toString().replace("'","''")+"' (FORMAT PARQUET)");
        }
        java.util.Properties manifest=new java.util.Properties();
        manifest.setProperty("format_version","1"); manifest.setProperty("provider","OpenTargets"); manifest.setProperty("gene_id_type","ensembl");
        manifest.setProperty("data_format","opentargets"); manifest.setProperty("data_file","associations");
        manifest.setProperty("sha256.associations/part.parquet",EnrichmentDownloads.sha256(parquet));
        try(var writer=Files.newBufferedWriter(target.resolve("manifest.properties"))) { manifest.store(writer,"test"); }
        assertEquals(Set.of("ENSG1"),EnrichmentDatabase.load(target,"OpenTargets",Map.of()).genes.get("EFO:1"));
        assertEquals(Set.of("ENSG1","ENSG2"),EnrichmentDatabase.load(target,"OpenTargets",Map.of("--association-min-score","0.1")).genes.get("EFO:1"));
        Files.delete(parquet);
        try(var connection=java.sql.DriverManager.getConnection("jdbc:duckdb:"); var query=connection.createStatement()) {
            query.execute("COPY (SELECT 'ENSG3' AS targetId,'EFO:2' AS diseaseId,0.8 AS associationScore) TO '"+parquet.toString().replace("'","''")+"' (FORMAT PARQUET)");
        }
        manifest.setProperty("sha256.associations/part.parquet",EnrichmentDownloads.sha256(parquet));
        try(var writer=Files.newBufferedWriter(target.resolve("manifest.properties"))) { manifest.store(writer,"test"); }
        assertEquals(Set.of("ENSG3"),EnrichmentDatabase.load(target,"OpenTargets",Map.of()).genes.get("EFO:2"));
    }
}
