/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.function.DoubleSupplier;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.glm.GlmFamilies;
import org.jlinalg.glmm.MultidimensionalGlmmQuadrature;
import org.jlinalg.glmm.MultidimensionalQuadratureOptions;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.regression.ProductKernelRegression;
import org.jlinalg.sem.SemDwls;
import org.jlinalg.sem.SemMixed;
import org.jlinalg.sem.SemModel;
import org.jlinalg.timeseries.ArimaRegression;
import org.jlinalg.timeseries.ArimaOrder;
import org.jlinalg.timeseries.ArimaSmoothing;

/** Small reproducible workloads; reference gates precede timed fits. */
public final class EstimatorExtensionsBenchmark {
    private EstimatorExtensionsBenchmark() { }
    private static double[][] data(String name) throws Exception {
        return Files.readAllLines(Path.of("src/test/resources/estimator-extensions/"+name)).stream()
            .map(s->Arrays.stream(s.split("\\s+")).mapToDouble(Double::parseDouble).toArray()).toArray(double[][]::new);
    }
    private static double[] vector(String name)throws Exception{return Arrays.stream(data(name)).flatMapToDouble(Arrays::stream).toArray();}
    private static void time(String name,DoubleSupplier run,double expected,double tolerance) {
        double value=run.getAsDouble(),error=Math.abs(value-expected);
        if(!Double.isFinite(error)||error>tolerance)throw new AssertionError(name+" reference error "+error);
        for(int i=0;i<3;i++)run.getAsDouble();
        int repeats=10;double checksum=0;long start=System.nanoTime();
        for(int i=0;i<repeats;i++)checksum+=run.getAsDouble();
        System.out.printf(java.util.Locale.ROOT,"%s,%d,%.3f,%.3g,%.9g%n",name,repeats,(System.nanoTime()-start)/1e6/repeats,error,checksum);
    }
    public static void main(String[] args)throws Exception {
        System.out.println("workload,repeats,milliseconds_per_run,reference_error,checksum");
        double[][] mixed=data("mixed-data.tsv");
        var model=SemModel.builder("x","y").meanStructure().intercept("mu","x",.3).variance("vx","x",1).fixedVariance("y",1).fixedIntercept("y",0).regression("b","y","x",.4).build();
        time("mixed_MAR_140",()->{var f=SemMixed.fit(mixed,new int[]{0,2},model);if(!f.informationAvailable())throw new AssertionError("mixed convergence");return f.logLikelihood();},vector("mixed-reference.tsv")[4],1e-7);
        var b=SemModel.builder("y1","y2","y3","y4").latent("F").fixedVariance("F",1);
        for(int i=1;i<=4;i++)b.loading("l"+i,"y"+i,"F",.6).fixedVariance("y"+i,1);
        var ordinal=b.build();var moments=new SemDwls.Moments(600,new int[]{3,3,3,3},vector("dwls-moments.tsv"),vector("dwls-gamma.tsv"));
        time("WLSMV_14_moments",()->{var f=SemDwls.fit(moments,ordinal);if(!f.converged())throw new AssertionError("DWLS convergence");return f.wlsmvStatistic();},vector("dwls-reference.tsv")[4],2e-5);
        double[][] source=data("arima-data.tsv"),x=new double[source.length][2];double[] y=new double[source.length];
        for(int i=0;i<y.length;i++){y[i]=source[i][0];x[i][0]=source[i][1];x[i][1]=source[i][2];}
        time("AR1_regression_80",()->{var f=ArimaRegression.fit(y,x,new ArimaOrder(1,0,0),null);if(!f.converged())throw new AssertionError("ARIMA convergence");return f.logLikelihood();},vector("arima-reference.tsv")[4],1e-7);
        time("diffuse_smoothing_80",()->ArimaSmoothing.smooth(y,new double[]{.4},new double[]{.2},new double[]{1,-1},1).states()[0][0],data("smooth-integrated.tsv")[0][0],2e-5);
        double[] bin={0,1,0,1,1,1,0,1};double[][] bx=new double[8][1];for(double[] row:bx)row[0]=1;
        var groups=List.of(RandomEffectTerm.randomIntercept("g",List.of("a","a","a","a","b","b","b","b")));
        double tensor=MultidimensionalGlmmQuadrature.evaluate(bin,bx,groups,null,GlmFamilies.binomial(),new double[]{0},new double[]{1},MultidimensionalQuadratureOptions.defaults(),BackendPolicy.CPU).logLikelihood();
        time("non_tensor_2D",()->{var f=MultidimensionalGlmmQuadrature.evaluate(bin,bx,groups,null,GlmFamilies.binomial(),new double[]{0},new double[]{1},MultidimensionalQuadratureOptions.nonTensorDefaults(),BackendPolicy.CPU);if(!f.converged())throw new AssertionError("QMC refinement");return f.logLikelihood();},tensor,2e-4);
        double[][] kx=new double[120][2];double[] ky=new double[120];
        for(int i=0;i<120;i++){kx[i]=new double[]{(i%60)/59.0,i/60};ky[i]=1+2*kx[i][0]+3*kx[i][1];}
        var types=new ProductKernelRegression.Type[]{ProductKernelRegression.Type.CONTINUOUS,ProductKernelRegression.Type.UNORDERED};
        time("mixed_kernel_120",()->ProductKernelRegression.predict(kx,ky,new double[]{.5,1},types,new double[]{.2,0}),5,1e-10);
    }
}
