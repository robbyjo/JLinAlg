"""Compare deterministic spatial kernels with a frozen Git baseline in isolated JVMs.

No production files/dependencies are replaced. Requires a JDK (java and javac).
Example: python tools/benchmark_cell_spatial.py --javac /path/to/javac
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import statistics
import subprocess

ROOT=Path(__file__).resolve().parents[1]
SOURCE="src/main/java/org/jlinalg/spatial/SpatialStatistics.java"
HARNESS=r'''
import java.util.*;
import org.jlinalg.spatial.SpatialStatistics;
public class CellSpatialTiming {
  static double[][] coordinates=new double[600][2];
  static String[] group=new String[600];
  static double[] values=new double[800];
  static String[] strata=new String[800];
  static List<SpatialStatistics.Edge> chain=new ArrayList<>();
  static {
    Random r=new Random(921);
    for(int i=0;i<600;i++){coordinates[i][0]=r.nextDouble();coordinates[i][1]=r.nextDouble();group[i]="section";}
    for(int i=0;i<800;i++){values[i]=r.nextGaussian();strata[i]="all";if(i>0)chain.add(new SpatialStatistics.Edge(i-1,i,1));}
  }
  static String run(String name) {
    if(name.equals("permutations")) {
      var result=SpatialStatistics.test(values,chain,strata,5000,17);
      return result.moran()+","+result.geary()+","+result.moranP()+","+result.gearyP();
    }
    var edges=SpatialStatistics.graph(coordinates,group,name,2,6);long checksum=0;
    for(var e:edges)checksum+=e.source()*600L+e.target();
    return edges.size()+","+checksum;
  }
  public static void main(String[] args) {
    for(String name:List.of("radius","knn","permutations")) {
      for(int warmup=0;warmup<3;warmup++)run(name);
      for(int repeat=0;repeat<7;repeat++) {
        long start=System.nanoTime();String check=run(name);double ms=(System.nanoTime()-start)/1e6;
        System.out.println(name+"\t"+ms+"\t"+check);
      }
    }
  }
}
'''

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--baseline",default="3a8b4f12eee31db935efbada15150bdbd94e9830")
    p.add_argument("--java",default="java")
    p.add_argument("--javac",default="javac")
    p.add_argument("--out")
    args=p.parse_args()
    output=ROOT/(args.out or "build/cell-spatial-benchmark-"+datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S"))
    output.mkdir(parents=True,exist_ok=False)
    frozen=subprocess.run(["git","show",args.baseline+":"+SOURCE],cwd=ROOT,check=True,capture_output=True,text=True).stdout
    report={"baseline":args.baseline,"workloads":{"radius":"600 nodes, complete radius graph","knn":"600 nodes, k=6, dense candidate set","permutations":"800 nodes, 799 edges, 5000 random-label permutations"},"timings":{}}
    checks={}
    for version,source in (("baseline",frozen),("revised",(ROOT/SOURCE).read_text(encoding="utf-8"))):
        work=output/version;work.mkdir()
        implementation=work/"SpatialStatistics.java";implementation.write_text(source,encoding="utf-8")
        harness=work/"CellSpatialTiming.java";harness.write_text(HARNESS,encoding="utf-8")
        subprocess.run([args.javac,"--release","17","-d",str(work),str(implementation),str(harness)],check=True,cwd=ROOT,timeout=60)
        run=subprocess.run([args.java,"-Xms256m","-Xmx1g","-cp",str(work),"CellSpatialTiming"],check=True,cwd=ROOT,capture_output=True,text=True,timeout=120)
        (work/"timings.tsv").write_text(run.stdout,encoding="utf-8")
        times={}
        for line in run.stdout.splitlines():
            name,ms,check=line.split("\t");times.setdefault(name,[]).append(float(ms));checks[version,name]=[float(x) for x in check.split(",")]
        report["timings"][version]={name:statistics.median(ms) for name,ms in times.items()}
    for name in report["workloads"]:
        for before,after in zip(checks["baseline",name],checks["revised",name]):
            assert abs(before-after)<=1e-12*max(1,abs(before)),(name,before,after)
    report["median_speedup"]={name:report["timings"]["baseline"][name]/report["timings"]["revised"][name] for name in report["workloads"]}
    report["java_version"]=subprocess.run([args.java,"-version"],capture_output=True,text=True,check=True).stderr.strip()
    (output/"benchmark.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2));print(output)

if __name__=="__main__":main()
