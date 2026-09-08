package org.jlinalg.mediation;
import static org.junit.jupiter.api.Assertions.*;
import static org.jlinalg.nonlinear.RemainingAuditData.*;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mixed.RandomEffectTerm;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;
public class MediationAuditTest {
 @Test void extremeConfidenceLevelsKeepFiniteNormalIntervals() throws Exception {
  for(double[] row:data("mediation-tail")) {
   var f=MediationAnalysis.indirectEffect(1,2,.04,.09,row[0]);
   assertEquals(2+.5*row[1],f.confidenceUpper(),1e-14);
   assertEquals(2-.5*row[1],f.confidenceLower(),1e-14);
  }
 }
 @Test void sobelCalculationPreservesReciprocalUnits() {
  var f=MediationAnalysis.indirectEffect(1e160,1e-160,1e300,1e-300,.95);
  assertEquals(1,f.estimate(),1e-14);assertEquals(1e10,f.standardError(),1e-5);
  assertThrows(IllegalArgumentException.class,()->MediationAnalysis.indirectEffect(1,1,-1,1,.95));
  var zero=MediationAnalysis.indirectEffect(0,0,1,1,.95);assertTrue(Double.isNaN(zero.pValue()));
 }
 @Test void ordinaryAndMixedSobelMatchIndependentRPaths() throws Exception {
  double[][] d=data("mediation");double[] y=new double[d.length],t=y.clone(),m=y.clone();double[][] c=new double[d.length][1];List<String> groups=new ArrayList<>();
  for(int i=0;i<d.length;i++){y[i]=d[i][0];t[i]=d[i][1];m[i]=d[i][2];c[i][0]=d[i][3];groups.add("g"+(int)d[i][4]);}
  var o=MediationAnalysis.fit(y,t,m,c,org.jlinalg.ols.OlsOptions.defaults(),BackendPolicy.CPU);check("mediation",o.aPath(),o.bPath(),o.directEffect(),o.totalEffect(),o.indirectEffect());
  assertEquals(o.totalEffect().estimate(),o.directEffect().estimate()+o.indirectEffect().estimate(),1e-12);
  var f=MediationAnalysis.fitMixed(y,t,m,c,List.of(RandomEffectTerm.randomIntercept("g",groups)),RemlOptions.defaults(),.9,BackendPolicy.CPU);
  assertTrue(f.converged());check("mediation_mixed",f.aPath(),f.bPath(),f.directEffect(),f.totalEffect(),f.indirectEffect());
  double critical=jdistlib.T.quantile(.95,f.aPath().degreesOfFreedom(),true,false);
  assertEquals(f.aPath().estimate()+critical*f.aPath().standardError(),f.aPath().confidenceUpper(),1e-12);
  var extreme=MediationAnalysis.fitMixed(y,t,m,c,List.of(RandomEffectTerm.randomIntercept("g",groups)),RemlOptions.defaults(),Math.nextDown(1.0),BackendPolicy.CPU);
  double[] tail=data("mediation-tail")[2];
  for(var effect:List.of(extreme.aPath(),extreme.bPath(),extreme.directEffect(),extreme.totalEffect())) {
   double df=effect.degreesOfFreedom();
   assertTrue(df==Math.rint(df)&&df>=114&&df<=117||df==Double.POSITIVE_INFINITY,"unexpected reference df="+df);
   double q=tail[df==Double.POSITIVE_INFINITY?1:2+(int)df-114];
   assertEquals(effect.estimate()+q*effect.standardError(),effect.confidenceUpper(),1e-12);
  }
 }
 private static void check(String name,MediationEffect a,MediationEffect b,MediationEffect direct,MediationEffect total,MediationEffect indirect)throws Exception {
  double[] values={a.estimate(),b.estimate(),direct.estimate(),total.estimate(),indirect.estimate(),indirect.standardError()};
  for(int j=0;j<values.length;j++)assertEquals(ref(name,"effects",j),values[j],2e-6,name+" "+j);
 }
 public static void main(String[] args)throws Exception{run(MediationAuditTest.class);}
}
