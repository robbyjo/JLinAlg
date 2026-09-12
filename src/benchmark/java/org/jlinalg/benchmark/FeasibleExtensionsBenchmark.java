/* SPDX-License-Identifier: GPL-2.0-or-later */
package org.jlinalg.benchmark;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Properties;
import java.util.function.Supplier;
import org.jlinalg.regression.QuantileRegressionInference;
import org.jlinalg.sem.SemModel;
import org.jlinalg.sem.SemOrdinal;
import org.jlinalg.timeseries.ArimaOptions;
import org.jlinalg.timeseries.ArimaOrder;
import org.jlinalg.timeseries.DiffuseArima;

/** Full fit plus covariance, frozen inputs and accuracy gates outside timing.
 * No I/O, R startup, or reference generation is included in measured iterations. */
public final class FeasibleExtensionsBenchmark {
    private static volatile double checksum;
    private FeasibleExtensionsBenchmark() { }
    public static void main(String[] args) throws Exception {
        System.out.println("case,median_ms,max_absolute_covariance_error,checksum");
        Path root=Path.of("src/test/resources/feasible-extensions");
        List<String> rows=Files.readAllLines(Path.of("src/test/resources/exact-quantile/regular.csv"));
        double[] y=new double[rows.size()];double[][] x=new double[rows.size()][];
        for(int i=0;i<y.length;i++) {
            double[] row=Arrays.stream(rows.get(i).split(",")).mapToDouble(Double::parseDouble).toArray();
            y[i]=row[0];x[i]=Arrays.copyOfRange(row,1,row.length);
        }
        Properties qr=reference(root.resolve("quantile.properties"));double[] densities=numbers(qr,"densities");
        run("quantile-supplied-density-n"+y.length,()->QuantileRegressionInference.fitExact(y,x,.35,densities).coefficientCovariance(),numbers(qr,"covariance"),1e-10);
        run("quantile-iid-kernel-n"+y.length,()->QuantileRegressionInference.fitExactIidKernel(y,x,.35,.45).coefficientCovariance(),numbers(qr,"kernelCovariance"),2e-8);

        int[][] ordinal=Files.readAllLines(root.resolve("ordinal-missing.csv")).stream()
            .map(row->Arrays.stream(row.split(",")).mapToInt(Integer::parseInt).toArray()).toArray(int[][]::new);
        SemModel model=SemModel.builder("z1","z2","z3","z4").latent("f")
            .loading("a","z1","f",.8).loading("b","z2","f",.8).loading("c","z3","f",.8).loading("d","z4","f",.8)
            .fixedVariance("f",1).fixedVariance("z1",1).fixedVariance("z2",1).fixedVariance("z3",1).fixedVariance("z4",1).build();
        run("ordinal-missing-n"+ordinal.length,()->{
            var f=SemOrdinal.fitPairwiseMissing(ordinal,new int[]{3,3,3,3},model);
            if(!f.converged() || !f.informationAvailable())throw new AssertionError("ordinal inference unavailable");
            return f.parameterCovariance();
        },numbers(reference(root.resolve("ordinal-missing.properties")),"covariance"),4e-5);

        double[] series=Files.readAllLines(Path.of("src/test/resources/timeseries/drift-series.csv")).stream().skip(1).mapToDouble(Double::parseDouble).toArray();
        run("diffuse-ar1-drift-missing-n"+series.length,()->{
            var f=DiffuseArima.fit(series,new ArimaOrder(1,1,0),ArimaOptions.builder().includeDrift(true).build());
            if(!f.coefficientInferenceAvailable())throw new AssertionError("diffuse inference unavailable");
            return f.coefficientCovariance();
        },numbers(reference(root.resolve("drift.properties")),"covariance"),2e-5);
    }
    private static void run(String name,Supplier<double[]> workload,double[] reference,double tolerance) {
        double[] checked=workload.get();double error=0;
        if(checked.length!=reference.length)throw new AssertionError("covariance dimensions differ");
        for(int i=0;i<checked.length;i++)error=Math.max(error,Math.abs(checked[i]-reference[i]));
        if(!(error<=tolerance))throw new AssertionError(name+" covariance error="+error);
        int warmups=Integer.getInteger("jlinalg.benchmark.warmups",3),measurements=Integer.getInteger("jlinalg.benchmark.measurements",7);
        if(warmups<0 || measurements<1)throw new IllegalArgumentException("invalid benchmark repetitions");
        double[] times=new double[measurements];
        for(int i=-warmups;i<measurements;i++) {
            long start=System.nanoTime();double[] result=workload.get();long elapsed=System.nanoTime()-start;
            for(double v:result)checksum+=v;
            if(i>=0)times[i]=elapsed/1e6;
        }
        Arrays.sort(times);
        System.out.printf(Locale.ROOT,"%s,%.6f,%.9g,%.9g%n",name,times[measurements/2],error,checksum);
    }
    private static Properties reference(Path path) throws Exception {
        Properties p=new Properties();try(var in=Files.newInputStream(path)){p.load(in);}return p;
    }
    private static double[] numbers(Properties p,String key){return Arrays.stream(p.getProperty(key).split(",")).mapToDouble(Double::parseDouble).toArray();}
}
