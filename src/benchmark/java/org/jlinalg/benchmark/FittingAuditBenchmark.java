package org.jlinalg.benchmark;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;
import java.util.function.DoubleSupplier;
import java.lang.management.ManagementFactory;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.*;
import org.jlinalg.ols.*;
import org.jlinalg.gee.*;
import org.jlinalg.penalized.*;

/** Warm CPU audit benchmark; every timed estimate contributes to a checksum. */
public final class FittingAuditBenchmark {
    private FittingAuditBenchmark() { }
    public static void main(String[] args) throws IOException {
        int n=args.length==0?3000:Integer.parseInt(args[0]);
        int[] permutation=new int[60], folds=new int[60];
        for(int i=0;i<60;i++)permutation[i]=i;
        Random random=new Random(17);
        for(int i=59;i>0;i--){int j=random.nextInt(i+1);int tmp=permutation[i];permutation[i]=permutation[j];permutation[j]=tmp;}
        StringBuilder foldText=new StringBuilder();
        for(int i=0;i<60;i++)folds[permutation[i]]=i%5;
        for(int f:folds)foldText.append(f).append('\n');
        Path output=Path.of("build/fitting-audit-v030/folds.txt");
        Files.createDirectories(output.getParent());Files.writeString(output,foldText);
        double[][] x=new double[n][2], X=new double[n][3];
        double[] y=new double[n], w=new double[n], offset=new double[n], counts=new double[n];
        int[] id=new int[n], waves=new int[n];
        for(int i=0;i<n;i++) {
            x[i][0]=Math.sin(i*.37);x[i][1]=Math.cos(i*.19);
            X[i][0]=1;X[i][1]=x[i][0];X[i][2]=x[i][1];
            w[i]=1+i%5;offset[i]=.2*Math.sin(i*.11);
            y[i]=1+.7*x[i][0]-.4*x[i][1]+.3*Math.sin(i*1.7)+offset[i];
            counts[i]=(i*7+3)%11;id[i]=i/3;waves[i]=i%3;
        }
        GlmOptions glm=GlmOptions.builder().relativeTolerance(1e-10).maximumIterations(200).build();
        GeeOptions gee=GeeOptions.builder().correlation(GeeCorrelation.FIXED)
            .fixedAssociation(new double[][]{{1,.2,.2},{.2,1,.2},{.2,.2,1}})
            .fixedDispersion(1).build();
        ElasticNetOptions enet=ElasticNetOptions.builder().alpha(1).observationWeights(w).relativeTolerance(1e-14).build();
        bench("ols",n,()->sum(Ols.fit(y,X,w,offset,OlsOptions.defaults(),BackendPolicy.CPU).coefficients()));
        bench("poisson",n,()->{var f=Glm.fit(counts,X,GlmFamilies.poisson(),w,offset,glm,BackendPolicy.CPU);require(f.converged());return sum(f.coefficients());});
        bench("gee_fixed",n,()->{var f=Gee.fit(y,X,id,waves,GlmFamilies.gaussian(),w,offset,gee,BackendPolicy.CPU);require(f.converged());return sum(f.coefficients());});
        bench("lasso_path",n,()->{var p=PenalizedRegression.path(y,x,new double[]{.2,.08,.02},enet);double total=0;for(var f:p.fits()){require(f.converged());total+=f.intercept()+sum(f.coefficients());}return total;});
        bench("lasso_cv",n,()->{var p=PenalizedRegressionCrossValidation.fit(y,x,new double[]{.2,.08,.02},5,17,enet);return sum(p.meanSquaredErrors());});
    }
    private static double sum(double[] values){return Arrays.stream(values).sum();}
    private static void require(boolean valid){if(!valid)throw new IllegalStateException("benchmark fit did not converge");}
    private static void bench(String name,int n,DoubleSupplier task){
        int batch=20;
        for(int i=0;i<10;i++)require(Double.isFinite(task.getAsDouble()));
        var bean=(com.sun.management.ThreadMXBean)ManagementFactory.getThreadMXBean();
        long thread=Thread.currentThread().getId();
        long allocatedBefore=bean.getThreadAllocatedBytes(thread);
        double[] elapsed=new double[7];double checksum=0;
        for(int i=0;i<7;i++){long start=System.nanoTime();for(int j=0;j<batch;j++){double value=task.getAsDouble();require(Double.isFinite(value));checksum+=value;}elapsed[i]=(System.nanoTime()-start)/1e6/batch;}
        long allocated=bean.getThreadAllocatedBytes(thread)-allocatedBefore;
        Arrays.sort(elapsed);System.out.printf(Locale.ROOT,"%s n=%d median_ms=%.6f bytes_per_fit=%.1f checksum=%.17g%n",name,n,elapsed[3],(double)allocated/(7*batch),checksum);
    }
}
