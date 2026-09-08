package org.jlinalg.settest;
import java.util.Locale;
public class PipelineKernelAuditBenchmark {
 public static void main(String[] args) {
  System.out.println("method,rep,n,milliseconds,checksum");
  run("java-rank-two",20,new double[]{1,10},.16908111409965279,200);
  run("java-rank-three",100,new double[]{1,.5,.25},2.5052827321161237e-23,100);
 }
 private static void run(String method,double q,double[] lambda,double expected,int repeats) {
  for(int rep=-5;rep<9;rep++) {
   long start=System.nanoTime();double checksum=0;
   for(int i=0;i<repeats;i++)checksum+=QuadraticFormDistribution.survival(q,lambda).pValue();
   double ms=(System.nanoTime()-start)/(1e6*repeats);
   if(!Double.isFinite(checksum)||Math.abs(checksum/(repeats*expected)-1)>1e-10)throw new AssertionError("inaccurate kernel");
   System.out.printf(Locale.ROOT,"%s,%d,%d,%.9f,%.17g%n",method,rep,repeats,ms,checksum);
  }
 }
}
