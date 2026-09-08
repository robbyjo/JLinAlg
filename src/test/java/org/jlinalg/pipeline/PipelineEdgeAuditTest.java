package org.jlinalg.pipeline;
import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import org.jlinalg.association.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.ols.OlsOptions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
class PipelineEdgeAuditTest {
    @TempDir Path directory;
    @Test void streamedEventsPartitionInputAndStayInSourceOrder() throws Exception {
        Path path=directory.resolve("events.tsv");
        Files.writeString(path,"id\ts1\ts2\ts3\ts4\ts5\ts6\n"
            +"signal\t0\t0\t1\t1\t2\t2\nfiltered\t1\t1\t1\t1\t1\t1\n"
            +"collinear\t0\t1\t0\t1\t0\t1\n");
        var source=DelimitedVariantSource.open(path);List<String> events=new ArrayList<>();
        var summary=StreamingAssociationPipeline.fastOlsTo(source,source.metadata().sampleIds(),
            new double[]{1,1.2,2,2.1,3,3.2},new double[][]{{1,0},{1,1},{1,0},{1,1},{1,0},{1,1}},
            null,null,OlsOptions.defaults(),AssociationEngineOptions.defaults().withBackendPolicy(BackendPolicy.CPU)
                .withFailurePolicy(AssociationFailurePolicy.RECORD_NAN),
            new AssociationPipelineOptions(100,VariantFilterOptions.defaults()),new AssociationPipelineSink(){
                public void acceptEstimate(AssociationPipelineEstimate e){events.add("ok:"+e.variant().id());}
                public void acceptExcluded(VariantFilterResult e){events.add("filtered:"+e.variant().id());}
                public void acceptFailure(AssociationPipelineFailure e){events.add("failed:"+e.variantId());}
            });
        assertEquals(List.of("ok:signal","filtered:filtered","failed:collinear"),events);
        assertEquals(3,summary.sourceVariants());assertEquals(1,summary.testedVariants());
        assertEquals(1,summary.excludedVariants());assertEquals(1,summary.failures());
    }
    @Test void dosageVarianceDoesNotCancelAndInfinityIsNotMissing() {
        var v=new VariantRecord("v","",0,"","",new double[]{1-1e-8,1,1+1e-8},Double.NaN);
        assertEquals(1e-16,VariantStatistics.of(v).dosageVariance(),2e-24);
        assertTrue(VariantFilters.evaluate(v,VariantFilterOptions.defaults()).included());
        assertThrows(IllegalArgumentException.class,()->VariantStatistics.of(
            new VariantRecord("bad","",0,"","",new double[]{1,Double.POSITIVE_INFINITY},Double.NaN)));
    }
    @Test void dsWithoutGtRetainsDosageAndRequestedOrder() throws Exception {
        Path path=directory.resolve("ds.vcf");
        Files.writeString(path,"##fileformat=VCFv4.2\n##contig=<ID=1>\n"
            +"##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"Dosage\">\n"
            +"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tz\ta\n"
            +"1\t10\tv\tA\tG\t.\tPASS\t.\tDS\t0.25\t1.75\n");
        var source=VcfVariantSource.open(path);
        try(var reader=source.open(SampleAlignment.requireOrder(source.metadata().sampleIds(),List.of("a","z")))) {
            assertArrayEquals(new double[]{1.75,.25},reader.read(1).variants().get(0).dosages(),0);
        }
    }
    @Test void transformMustNotSilentlyDropSamplesOrMutateSource() throws Exception {
        Path path=directory.resolve("matrix.tsv");
        Files.writeString(path,"id\ts4\ts2\ts1\ts3\nf\t4\t2\t1\t3\n");
        var source=DelimitedMatrixSource.open(path);
        double[][] x={{1},{1},{1},{1}};
        assertThrows(IllegalArgumentException.class,()->StreamingOmicsAssociationPipeline.scanPredictors(
            source,List.of("s1","s2","s3","s4"),new double[]{1,3,2,5},x,
            values->Arrays.copyOf(values,3),OmicsMissingPolicy.ERROR,1,null,null,OlsOptions.defaults(),
            AssociationEngineOptions.defaults().withBackendPolicy(BackendPolicy.CPU)));
    }
    @Test void manuallyConstructedAlignmentRejectsAmbiguousRows() {
        assertThrows(IllegalArgumentException.class,()->new SampleAlignment(List.of("a","b"),new int[]{0,0},new int[]{0,1}));
        assertThrows(IllegalArgumentException.class,()->new SampleAlignment(List.of("a"),new int[]{-1},new int[]{0}));
        assertThrows(IllegalArgumentException.class,()->new SampleAlignment(List.of("a","a"),new int[]{0,1},new int[]{0,1}));
    }
}
