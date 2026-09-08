package org.jlinalg.association;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.ols.OlsOptions;
/** Exact alternative-model OLS scans, including mean imputation and inference. */
public class PipelineOlsAuditBenchmark {
 public static void main(String[] args) {
  int n=400,m=128;double[] y=new double[n];double[][] x=new double[n][2],g=new double[n][m];
  for(int i=0;i<n;i++) {
   y[i]=Math.sin(.13*i)+.3*Math.cos(.07*i);x[i][0]=1;x[i][1]=Math.cos(.07*i);
   for(int j=0;j<m;j++)g[i][j]=(i+j)%97==0?Double.NaN:Math.sin(.023*(i+1)*(j+1))+.2*Math.cos(.031*(i+1)*(j+1));
  }
  var options=AssociationEngineOptions.defaults().withBackendPolicy(BackendPolicy.CPU);
  System.out.println("method,rep,n,markers,milliseconds,checksum");
  for(int rep=-30;rep<9;rep++) {
   long start=System.nanoTime();
   var result=FastOlsAssociation.scanPredictors(y,x,g,null,null,null,OlsOptions.defaults(),options);
   double ms=(System.nanoTime()-start)/1e6,checksum=0;
   if(!result.failures().isEmpty())throw new AssertionError("failed fits");
   for(int j=0;j<m;j++) {
    var e=result.estimate(j);
    if(!Double.isFinite(e.beta())||!Double.isFinite(e.standardError())||!Double.isFinite(e.pValue()))throw new AssertionError("nonfinite fit");
    checksum+=(j+1)*(e.beta()+e.standardError()+e.pValue());
   }
   if(!Double.isFinite(checksum)||Math.abs(checksum-8075.531253288807)>1e-8)throw new AssertionError("OLS checksum differs from frozen R reference: "+checksum);
   System.out.printf(Locale.ROOT,"java-ols,%d,%d,%d,%.9f,%.17g%n",rep,n,m,ms,checksum);
  }
 }
}
