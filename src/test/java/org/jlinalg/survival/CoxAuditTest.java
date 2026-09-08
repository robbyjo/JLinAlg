package org.jlinalg.survival;
import static org.junit.jupiter.api.Assertions.*;
import static org.jlinalg.nonlinear.RemainingAuditData.*;
import java.util.*;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.mixed.*;
import org.junit.jupiter.api.Test;
public class CoxAuditTest {
 @Test void disjointRiskSetsRecoverSmallRisksAfterLargeOffsetsLeave() {
  var survival=new CoxSurvivalData(new double[]{0,0,2,2},new double[]{1,1,3,3},new boolean[]{true,false,true,false},null);
  for(double[] offset:new double[][]{{40,40,0,0},{1000,1000,0,0},{0,0,1000,1000},{1e100,1e100,1e100,1e100}})for(CoxTies ties:CoxTies.values()) {
   var e=CoxPartialLikelihood.evaluate(survival,new double[]{0,1,0,1},1,new double[]{0},offset,ties);
   assertEquals(-2*Math.log(2),e.logLikelihood(),1e-12);assertEquals(-1,e.score()[0],1e-12);assertEquals(.5,e.information()[0],1e-12);
  }
  var baseline=CoxPartialLikelihood.baseline(survival,new double[]{0,1,0,1},1,new double[]{0},new double[]{40,40,0,0},CoxTies.BRESLOW);
  assertEquals(.5*Math.exp(-40),baseline.get(0).hazardIncrement(),1e-30);
  assertEquals(.5,baseline.get(1).hazardIncrement(),1e-14);
 }
 @Test void tiedStratifiedCountingFitsMatchSurvivalR() throws Exception {
  double[][] d=data("cox");double[] start=new double[d.length],stop=start.clone(),offset=start.clone();boolean[] event=new boolean[d.length];int[] strata=new int[d.length];double[][] x=new double[d.length][2];
  for(int i=0;i<d.length;i++){start[i]=d[i][0];stop[i]=d[i][1];event[i]=d[i][2]==1;strata[i]=(int)d[i][3];x[i][0]=d[i][4];x[i][1]=d[i][5];offset[i]=d[i][6];}
  for(CoxTies ties:CoxTies.values()) {
   var f=CoxRegression.fit(new CoxSurvivalData(start,stop,event,strata),x,offset,CoxOptions.defaults().withTies(ties),BackendPolicy.CPU);
   assertTrue(f.converged());String name="cox_"+ties.toString().toLowerCase(java.util.Locale.ROOT);
   for(int j=0;j<2;j++){assertEquals(ref(name,"beta",j),f.beta()[j],1e-7);assertEquals(ref(name,"se",j),f.standardErrors()[j],1e-8);}
   assertEquals(ref(name,"ll",0),f.logPartialLikelihood(),1e-9);
  }
 }
 @Test void splittingAnIntervalDoesNotChangeTheRiskSetLikelihood() {
  int n=60;double[] stop=new double[n];boolean[] event=new boolean[n];double[][] x=new double[n][1];
  double[] begin2=new double[2*n],stop2=new double[2*n];boolean[] event2=new boolean[2*n];double[][] x2=new double[2*n][1];
  for(int i=0;i<n;i++){stop[i]=i+1;event[i]=i%4!=0;x[i][0]=Math.sin(i*1.7);stop2[2*i]=stop[i]/2;begin2[2*i+1]=stop[i]/2;stop2[2*i+1]=stop[i];event2[2*i+1]=event[i];x2[2*i][0]=x2[2*i+1][0]=x[i][0];}
  var a=CoxRegression.fit(CoxSurvivalData.rightCensored(stop,event),x,null,CoxOptions.defaults(),BackendPolicy.CPU);
  var b=CoxRegression.fit(new CoxSurvivalData(begin2,stop2,event2,null),x2,null,CoxOptions.defaults(),BackendPolicy.CPU);
  assertTrue(a.converged()&&b.converged());assertEquals(a.beta()[0],b.beta()[0],1e-9);assertEquals(a.logPartialLikelihood(),b.logPartialLikelihood(),1e-10);assertEquals(a.standardErrors()[0],b.standardErrors()[0],1e-9);
 }
 @Test void sparseFrailtySurvivesIrrelevantEarlyCensorOffsets() {
  double[] time={1,2,3,4,5,6,7,8};boolean[] event={false,false,true,true,false,true,true,false};double[][] x={{0},{1},{-1},{2},{0},{-2},{1},{-1}};
  var survival=CoxSurvivalData.rightCensored(time,event);
  var term=RandomEffectTerm.randomIntercept("g",List.of("a","a","b","b","c","c","d","d"));
  var options=new CoxMixedOptions(CoxOptions.defaults(),new double[]{.2},2,1e-4,.199999,.200001);
  try(var prepared=SparseCoxMixedModel.prepare(survival,term,SparsePrecisionMatrix.identity(4),List.of(),options,BackendPolicy.CPU)) {
   var baseline=prepared.fitAtVariances(x,null,.2);assertTrue(baseline.converged());
   for(double shift:new double[]{40,1000}) {
    var f=prepared.fitAtVariances(x,new double[]{shift,shift,0,0,0,0,0,0},.2);
    assertTrue(f.converged());assertEquals(baseline.beta()[0],f.beta()[0],1e-8);
    assertArrayEquals(baseline.randomEffects("g").modes(),f.randomEffects("g").modes(),1e-8);
    assertEquals(baseline.standardErrors()[0],f.standardErrors()[0],1e-8);
   }
  }
 }
 @Test void fixedCoxInferenceAndConvergenceAreInvariantToPredictorUnits() {
  double[] time={1,2,3,4,5,6,7,8};boolean[] event={false,false,true,true,false,true,true,false};double[][] x={{0},{1},{-1},{2},{0},{-2},{1},{-1}};
  var survival=CoxSurvivalData.rightCensored(time,event);var base=CoxRegression.fit(survival,x,null,CoxOptions.defaults(),BackendPolicy.CPU);
  for(double scale:new double[]{1e-12,1e12}) {
   double[][] transformed=Arrays.stream(x).map(r->new double[]{r[0]*scale}).toArray(double[][]::new);
   var f=CoxRegression.fit(survival,transformed,null,CoxOptions.defaults(),BackendPolicy.CPU);
   assertTrue(f.converged());assertEquals(base.beta()[0],f.beta()[0]*scale,1e-10);assertEquals(base.standardErrors()[0],f.standardErrors()[0]*scale,1e-10);
   assertEquals(base.logPartialLikelihood(),f.logPartialLikelihood(),1e-10);
  }
 }
 @Test void gammaModesAndLaplaceNormalizationMatchIndependentRObjective() throws Exception {
  double[][] d=data("gamma_frailty");double[] time=new double[d.length];boolean[] event=new boolean[d.length];double[][] x=new double[d.length][1];List<String> groups=new ArrayList<>();
  for(int i=0;i<d.length;i++){time[i]=d[i][0];event[i]=d[i][1]==1;x[i][0]=d[i][2];groups.add("g"+(int)d[i][3]);}
  var f=CoxGammaFrailty.fitAtVariance(CoxSurvivalData.rightCensored(time,event),x,groups,null,.4,CoxGammaFrailtyOptions.defaults(),BackendPolicy.CPU);
  assertTrue(f.converged());assertEquals(ref("gamma_frailty","beta",0),f.beta()[0],2e-6);
  for(int j=0;j<4;j++)assertEquals(ref("gamma_frailty","mode",j),f.logFrailties()[j],2e-6);
  assertEquals(ref("gamma_frailty","partial",0),f.partialLogLikelihood(),1e-6);
  assertEquals(ref("gamma_frailty","laplace",0),f.laplaceLogLikelihood(),2e-6);
 }
 @Test void countingScoreAndInformationMatchFiniteDifferences() throws Exception {
  double[][] d=data("cox");double[] start=new double[d.length],stop=start.clone(),offset=start.clone(),x=new double[2*d.length];boolean[] event=new boolean[d.length];int[] strata=new int[d.length];
  for(int i=0;i<d.length;i++){start[i]=d[i][0];stop[i]=d[i][1];event[i]=d[i][2]==1;strata[i]=(int)d[i][3];x[2*i]=d[i][4];x[2*i+1]=d[i][5];offset[i]=d[i][6];}
  var s=new CoxSurvivalData(start,stop,event,strata);double[] beta={.2,-.1};double h=1e-5;
  for(CoxTies ties:CoxTies.values()) {
   var center=CoxPartialLikelihood.evaluate(s,x,2,beta,offset,ties);
   for(int j=0;j<2;j++) {
    double[] a=beta.clone(),b=beta.clone();a[j]+=h;b[j]-=h;
    var plus=CoxPartialLikelihood.evaluate(s,x,2,a,offset,ties);var minus=CoxPartialLikelihood.evaluate(s,x,2,b,offset,ties);
    assertEquals(center.score()[j],(plus.logLikelihood()-minus.logLikelihood())/(2*h),2e-7);
    for(int k=0;k<2;k++)assertEquals(center.information()[k*2+j],-(plus.score()[k]-minus.score()[k])/(2*h),2e-7);
   }
  }
 }
 public static void main(String[] args)throws Exception{run(CoxAuditTest.class);}
}
