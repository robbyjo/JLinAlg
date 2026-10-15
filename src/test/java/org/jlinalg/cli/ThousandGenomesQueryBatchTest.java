package org.jlinalg.cli;

import static org.junit.jupiter.api.Assertions.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.jlinalg.cli.StudyGeneticEvidence.Variant;
import org.junit.jupiter.api.Test;

class ThousandGenomesQueryBatchTest {
    private static Variant variant(String id,long position,String effect,String other) {
        return new Variant(id,"1",position,effect,other,.1,.01,1e-20);
    }
    private static String row(long position,String id,String ref,String alt) {
        return "chr1\t"+position+"\t"+id+"\t"+ref+"\t"+alt+"\t.\tPASS\t.\tGT\t0|1";
    }
    private static ThousandGenomesGrch38Ld.RecordLines lines(List<String> rows) {
        var iterator=rows.iterator();return ()->iterator.hasNext()?iterator.next():null;
    }

    @Test void denseRequestsUseOneQueryAndRetainInputIndexWithoutRetainingVcfRows() throws Exception {
        List<Variant> requested=new ArrayList<>();List<String> records=new ArrayList<>();
        for(int i=1;i<=10000;i++){requested.add(0,variant("rs"+i,i*10L,"A","G"));records.add(row(i*10L,"rs"+i,"G","A"));}
        int[] callbacks={0},queries={0};
        int count=ThousandGenomesGrch38Ld.collectRecords(requested,10,(start,end)->{
            queries[0]++;assertEquals(10,start);assertEquals(100000,end);return lines(records);
        },(i,fields)->{assertEquals(requested.get(i).id(),fields[2]);callbacks[0]++;});
        assertEquals(1,count);assertEquals(1,queries[0]);assertEquals(10000,callbacks[0]);
    }

    @Test void sparsePositionsAreSplitAtTabixBinSizedGapAndMissingQueriesRemainAbsent() throws Exception {
        var requested=List.of(variant("rs1",1,"A","G"),variant("rs2",16385,"A","G"),variant("rs3",32770,"A","G"));
        List<List<Long>> windows=new ArrayList<>();
        int count=ThousandGenomesGrch38Ld.collectRecords(requested,10,(start,end)->{windows.add(List.of(start,end));return null;},(i,row)->fail());
        assertEquals(2,count);assertEquals(List.of(List.of(1L,16385L),List.of(32770L,32770L)),windows);
    }

    @Test void exactRequestedPositionsKeepPointQueryIdOrAlleleMatchingAndIgnoreOtherRecords() throws Exception {
        var requested=List.of(variant("rs2",200,"T","C"),variant("alias",100,"A","G"),variant("missing",150,"A","C"));
        var records=List.of(row(100,"rs1","G","A"),row(125,"noise","C","A"),
            row(150,"noise","N","A"),row(150,"noise","C","A,G"),row(200,"rs2","C","T"));
        List<Integer> indices=new ArrayList<>();
        assertEquals(1,ThousandGenomesGrch38Ld.collectRecords(requested,10,(a,b)->lines(records),(i,row)->indices.add(i)));
        assertEquals(List.of(1,0),indices);
    }

    @Test void multipleRecordsAtSamePositionRemainSeparateAndAmbiguousMatchesFail() throws Exception {
        var requested=List.of(variant("rsA",100,"A","G"),variant("rsC",100,"C","G"));
        var records=List.of(row(100,"rsA","G","A"),row(100,"rsC","G","C"));
        List<Integer> indices=new ArrayList<>();
        ThousandGenomesGrch38Ld.collectRecords(requested,10,(a,b)->lines(records),(i,row)->indices.add(i));
        assertEquals(List.of(0,1),indices);
        var error=assertThrows(IOException.class,()->ThousandGenomesGrch38Ld.collectRecords(
            requested,10,(a,b)->lines(List.of(records.get(0),records.get(0))),(i,row)->{}));
        assertTrue(error.getMessage().contains("Multiple GRCh38 LD records match rsA"));
    }
}
