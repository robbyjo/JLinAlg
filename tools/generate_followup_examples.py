"""Regenerate synthetic omics-follow-up examples; no biological claims or downloads."""
from pathlib import Path
import csv
import random

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "examples" / "followup"
OUT.mkdir(parents=True, exist_ok=True)

def table(name, header, rows):
    with (OUT / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)

rng = random.Random(20260919)
genes = [f"G{i:02d}" for i in range(1, 25)]
for cohort in ("discovery", "replication"):
    rows, traits = [], []
    for i in range(48):
        first, second = rng.gauss(0, 1), rng.gauss(0, 1)
        sample = f"{cohort}_{i+1}"
        values = [(first if j < 12 else second) + .35*rng.gauss(0, 1) for j in range(24)]
        rows.append([sample] + [f"{x:.12g}" for x in values])
        traits.append([sample, f"{first+.5*rng.gauss(0,1):.12g}"])
    table(cohort+".tsv", ["sample_id"]+genes, rows)
    if cohort == "discovery": table("phenotype.tsv", ["sample_id", "trait"], traits)

table("edges.tsv", ["source","target","weight","type"],
      [[genes[i], genes[i+1], .8, "synthetic_physical"] for i in range(23)])
table("hits.tsv", ["gene"], [["G01"],["G02"],["G03"],["G13"],["UNMAPPED"]])
table("regulators.tsv", ["gene"], [["G01"],["G13"]])
table("signaling.tsv", ["source","target","sign"],
      [["I1","N1",1],["I2","N2",-1],["N1","M1",1],["N1","M2",1],["N2","M1",1],["N2","M2",-1]])
table("activities.tsv", ["gene","activity"], [["M1",1],["M2",1]])
keys=["genome_build","chrom","pos","ref","alt"]
variants=[["GRCh38","1",101,"A","G"],["GRCh38","1",202,"C","T"]]
table("variants.tsv",keys,variants)
table("annotations.tsv",keys+["gene","transcript","functional_score"],
      [variants[0]+["G01","TX1",.8],variants[0]+["G01","TX2",.4]])
table("evidence.tsv",["variant_id","gene","transcript","pip","functional_score","coloc_h4"],
      [["GRCh38:1:101:A:G","G01","TX1",.8,.8,.9],["GRCh38:1:101:A:G","G01","TX2",.8,.4,.9],["GRCh38:1:202:C:T","G02","NA",.2,"NA",.1]])
table("components.tsv",["column","weight","minimum","maximum","direction"],
      [["pip",2,0,1,"higher"],["functional_score",1,0,1,"higher"],["coloc_h4",1,0,1,"higher"]])
(OUT/"vep.tsv").write_text("##VEP=synthetic_fixture\n#Uploaded_variation\tGene\tFeature\tConsequence\tExtra\nGRCh38:1:101:A:G\tG01\tTX1\tmissense_variant\tIMPACT=MODERATE\n",encoding="utf-8")
(OUT/"jlinalg.yaml").write_text("""schema_version: 1
commands:
  network:
    method: sparse
    matrix: {path: discovery.tsv}
    lambda: 0.15
    rule: and
    seed: 19
  variant-annotate:
    input: {path: variants.tsv}
    genome-build: GRCh38
  variant-score:
    input: {path: evidence.tsv}
    components: {path: components.tsv}
""",encoding="utf-8")
