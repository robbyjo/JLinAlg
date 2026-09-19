"""Execute and verify all follow-up workflows; preserves artifacts under a fresh output directory.

Example: python tools/validate_followup.py --rscript /path/to/Rscript --r-library /path/to/library
"""
import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import math

ROOT = Path(__file__).resolve().parents[1]

def read(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--jar",default="build/cli/jlinalg-0.3.6.jar")
    parser.add_argument("--rscript",required=True)
    parser.add_argument("--r-library")
    parser.add_argument("--out")
    args = parser.parse_args()
    output = (ROOT / (args.out or "build/followup-check-"+datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S"))).resolve()
    output.mkdir(parents=True,exist_ok=False)
    example = ROOT/"examples/followup"
    jar=(ROOT/args.jar).resolve()
    def run(name, *arguments):
        dest=output/name
        command=["java","-jar",str(jar),"--no-config",*map(str,arguments),"--out",str(dest)]
        result=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,timeout=600)
        (output/(name+".console.txt")).write_text(result.stdout+result.stderr,encoding="utf-8")
        if result.returncode: raise RuntimeError(result.stdout+result.stderr)
        assert (dest/"manifest.yaml").is_file()
        assert Path(str(dest)+".config.yaml").is_file()
        return dest
    db=run("database","variant-db","--source-file",example/"annotations.tsv","--genome-build","GRCh38","--release","synthetic-1","--license","CC0")
    ann=run("annotation","variant-annotate","--input",example/"variants.tsv","--genome-build","GRCh38","--database",db)
    assert len(read(ann/"annotations.tsv"))==3
    cons=run("consequence","variant-consequence","--input",example/"variants.tsv","--genome-build","GRCh38","--engine","import","--annotations",example/"vep.tsv","--format","vep")
    assert read(cons/"consequences.tsv")[0]["consequence"]=="missense_variant"
    score=run("score","variant-score","--input",example/"evidence.tsv","--components",example/"components.tsv")
    assert math.isclose(float(read(score/"scores.tsv")[0]["priority_score"]),.825)
    ref=run("reference","network","--method","reference","--edges",example/"edges.tsv","--hits",example/"hits.tsv","--permutations",99)
    assert math.isclose(sum(float(x["propagation_score"]) for x in read(ref/"nodes.tsv")),1,abs_tol=1e-9)
    sparse=run("sparse","network","--method","sparse","--matrix",example/"discovery.tsv","--lambda",.15,"--bootstraps",20)
    assert read(sparse/"edges.tsv")
    assert all(0<=float(x["bootstrap_frequency"])<=1 for x in read(sparse/"edges.tsv"))
    diff=run("differential","network","--method","differential","--matrix",example/"discovery.tsv","--matrix-b",example/"replication.tsv")
    assert len(read(diff/"edges.tsv"))==24*23//2
    rargs=["--rscript",args.rscript]
    if args.r_library: rargs += ["--r-library",args.r_library]
    wgcna=run("wgcna","network","--method","wgcna","--matrix",example/"discovery.tsv","--power",6,"--min-module-size",5,"--pheno",example/"phenotype.tsv","--trait","trait","--reference-matrix",example/"replication.tsv","--permutations",10,*rargs)
    modules=read(wgcna/"modules.tsv")
    assert len({x["module"] for x in modules[:12]})==1
    assert len({x["module"] for x in modules[12:]})==1
    assert modules[0]["module"]!=modules[-1]["module"]
    assert max(abs(float(x["correlation"])) for x in read(wgcna/"module-trait.tsv"))>.8
    assert (wgcna/"preservation.tsv").is_file()
    reg=run("regulatory","network","--method","regulatory","--matrix",example/"discovery.tsv","--regulators",example/"regulators.tsv","--trees",100,*rargs)
    edges=read(reg/"edges.tsv")
    assert edges and all(x["source"] in {"G01","G13"} and x["sign"]=="unknown" for x in edges)
    signaling=run("signaling","network","--method","signaling","--edges",example/"signaling.tsv","--activities",example/"activities.tsv",*rargs)
    edges=read(signaling/"edges.tsv")
    assert edges and all("Perturbation" not in (x["source"],x["target"]) for x in edges)
    assert any(x["source"]=="I1" and x["target"]=="N1" and x["sign"]=="1" for x in edges)
    for directory in (wgcna,reg,signaling): assert (directory/"sessionInfo.txt").is_file()
    print("Validated all three variant operations and all six network methods:",output)

if __name__=="__main__": main()
