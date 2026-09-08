/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.genetics;

import java.nio.file.*;
import java.util.*;
import java.util.function.Supplier;
import org.jlinalg.coloc.*;
import org.jlinalg.susie.*;
import org.jlinalg.mr.*;
import org.jlinalg.compute.BackendPolicy;

/** Accuracy-gated warmed benchmark against independent base-R equations. */
public final class GeneticAuditBenchmark {
    private GeneticAuditBenchmark() { }
    private static double checksum;
    public static void main(String[] args) throws Exception {
        Properties ref=new Properties();
        try(var stream=Files.newInputStream(Path.of("src/test/resources/r-reference/genetic-audit-reference.properties"))){ref.load(stream);}
        double[] a=new double[64],b=new double[64],z=new double[32];
        for(int i=0;i<64;i++){a[i]=5*Math.sin((i+1)*.37)+.1*(i+1);b[i]=4*Math.cos((i+1)*.23)+.07*(i+1);}
        for(int i=0;i<32;i++)z[i]=5*Math.sin((i+1)*.37);
        var names=java.util.stream.IntStream.range(0,64).mapToObj(i->"s"+i).toList();
        var first=new ColocSusieInput(names,new double[][]{a});var second=new ColocSusieInput(names,new double[][]{b});
        Supplier<double[]> coloc=()->{var r=ColocSusie.analyze(first,second);return concat(r.signalPairs().get(0).hypothesisPosteriors(),r.sharedVariantPosterior());};
        double[][] ld=new double[32][32];for(int i=0;i<32;i++)ld[i][i]=1;
        var options=new SusieOptions(1,100,1e-8,.2,false,.95,0);
        Supplier<double[]> ser=()->{
            var r=Susie.fitSummary(z,ld,400,null,options,BackendPolicy.CPU);
            if(!r.converged())throw new AssertionError("benchmark fit did not converge");
            return concat(r.pip(),r.posteriorMean(),r.logBayesFactors());
        };
        double[][] dosages={{0,1,2},{0,Double.NaN,2},{1,2,0}};
        Supplier<double[]> grm=()->GenomicRelationshipMatrix.fromVariantDosages(dosages,List.of("a","b","c"),new GenomicRelationshipOptions(.01,.5),BackendPolicy.CPU).relationshipMatrix();
        Supplier<double[]> winner=()->new double[]{WinnerCurseCorrection.correct(40.01,1,40)};
        Supplier<double[]> phase=phaseSolver();
        var phaseRef=new Properties();
        try(var stream=Files.newInputStream(Path.of("src/test/resources/r-reference/genetic-final-fixes.properties"))){phaseRef.load(stream);}
        double phaseExpected=Double.parseDouble(phaseRef.getProperty("ld.slow.rSquared"));
        double phaseError=Math.abs(phase.get()[0]/phaseExpected-1);
        if(!Double.isFinite(phaseError)||phaseError>2e-10)throw new AssertionError("LD41 phase likelihood gate failed");
        System.out.printf(Locale.ROOT,"gate LD41 relative_error=%.9g%n",phaseError);
        check(ref,"coloc64",coloc.get(),3e-12);
        check(ref,"ser32",ser.get(),3e-12);
        check(ref,"grm",grm.get(),3e-12);
        double[] t={5.46,6,40.01,41},estimates=new double[4];
        for(int i=0;i<4;i++)estimates[i]=WinnerCurseCorrection.correct(t[i],1,i<2?5.45:40);
        check(ref,"winner",estimates,3e-11);
        System.out.println("accuracy gates passed: coloc64 H0-H4 + H4 SNPs, SER32 alpha/beta/logBF, GRM3, four conditional MLEs");
        timing("coloc64_Java",coloc,2000);
        timing("SER32_Java_IBSS_including_PSD_validation",ser,2000);
        timing("winner40_Java_score",winner,1000);
        timing("GRM3_Java_builder",grm,20000);
        timing("LD41_Java_stationary_core_reflection_included",phase,1000);
        System.out.println("Java "+System.getProperty("java.version")+" availableProcessors="+Runtime.getRuntime().availableProcessors());
    }
    private static Supplier<double[]> phaseSolver() throws ReflectiveOperationException {
        int[][] counts={{4,0,2},{3,21,4},{3,3,1}};
        byte[] a=new byte[41],b=new byte[41];boolean[] founders=new boolean[41];Arrays.fill(founders,true);
        int index=0;for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<counts[i][j];k++){a[index]=(byte)i;b[index++]=(byte)j;}
        // Benchmark the existing private core without changing the public API.
        // Reflection dispatch and genotype counting are included; BED I/O is not.
        var method=PlinkBedLdClumper.class.getDeclaredMethod("haplotypeRSquared",byte[].class,byte[].class,boolean[].class);
        method.setAccessible(true);
        return ()->{try{return new double[]{(double)method.invoke(null,a,b,founders)};}
            catch(ReflectiveOperationException e){throw new IllegalStateException("LD core benchmark failed",e);}};
    }
    private static double[] concat(double[]... arrays){int n=0;for(var a:arrays)n+=a.length;double[] r=new double[n];int i=0;for(var a:arrays){System.arraycopy(a,0,r,i,a.length);i+=a.length;}return r;}
    private static void check(Properties ref,String key,double[] actual,double tolerance){
        double max=0;
        for(int i=0;i<actual.length;i++){
            double expected=Double.parseDouble(ref.getProperty(key+"."+i));
            double error=Math.abs(actual[i]-expected)/Math.max(1,Math.abs(expected)); max=Math.max(max,error);
            if(!Double.isFinite(actual[i])||error>tolerance)throw new AssertionError(key+"["+i+"]: "+actual[i]+" vs "+expected);
        }
        System.out.printf(Locale.ROOT,"gate %s max_scaled_absolute_error=%.9g%n",key,max);
    }
    private static double consume(double[] values){double s=0;for(int i=0;i<values.length;i++)s+=(i+1)*values[i]*values[i];return s;}
    private static void timing(String label,Supplier<double[]> operation,int repeats){
        for(int i=0;i<1000;i++)checksum+=consume(operation.get());
        double[] times=new double[5];checksum=0;
        for(int batch=0;batch<5;batch++){
            long start=System.nanoTime();for(int i=0;i<repeats;i++)checksum+=consume(operation.get());
            times[batch]=(System.nanoTime()-start)*1e-6/repeats;
        }
        Arrays.sort(times);
        System.out.printf(Locale.ROOT,"%s median_ms=%.9g checksum=%.17g calls=%d%n",label,times[2],checksum,5*repeats);
    }
}
