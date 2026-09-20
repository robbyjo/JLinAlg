"""Deterministic synthetic replicated RNA/tissue inputs. No biological validation claim."""
from pathlib import Path
import csv
import random

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "examples/cell-spatial"
OUT.mkdir(parents=True, exist_ok=True)
rng = random.Random(814)

def write(name, header, rows):
    with (OUT / name).open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)

genes = ["MT-GENE"] + [f"GENE{i:02}" for i in range(1, 31)]
write("features.tsv", ["feature_id", "mitochondrial"], [[g, str(i == 0).lower()] for i, g in enumerate(genes)])
samples, paired, cells, counts = [], [], [], []
for s in range(8):
    case = s >= 4
    sid = f"S{s+1}"
    row = [sid, f"D{s+1}", "case" if case else "control", f"B{s%2}", "cortex", "scRNA", "baseline", 30+3*s+(s%3)]
    samples.append(row)
    paired.append([sid, f"P{s%4+1}", row[2], row[3], *row[4:6], "post" if case else "pre", row[-1]])
    for section in range(2):
        for c in range(20):
            oid = f"{sid}_{section}_{c}"
            kind = "T" if c < 8+(s%4)+int(case)*2 else "B"
            x, y = c%5*10, c//5*10
            cells.append([oid, sid, kind, f"spec-{sid}", f"sec-{sid}-{section}", "tissue", "all", x, y, x, "false"])
            for j, gene in enumerate(genes):
                level = 2 if j == 0 else 5+(j*3)%20
                level += (7 if kind == "T" and j in (1, 2, 3) else 0)
                level += (10 if case and kind == "T" and j in (1, 4, 5) else 0)
                level += x*.2 if case and j == 6 else 0
                value = max(0, round(level*rng.uniform(.5, 1.5)+(s%4)))
                if value: counts.append([oid, gene, value])
write("samples.tsv", ["sample_id", "donor_id", "condition", "batch", "source", "assay", "visit", "age"], samples)
write("paired-samples.tsv", ["sample_id", "donor_id", "condition", "batch", "source", "assay", "visit", "age"], paired)
write("cells.tsv", ["obs_id", "sample_id", "cell_type", "specimen_id", "section_id", "compartment", "stratum", "x", "y", "distance", "exclude"], cells)
write("counts.tsv", ["obs_id", "feature_id", "count"], counts)
write("selected-features.tsv", ["feature_id"], [[g] for g in genes[1:7]])
write("gene-sets.tsv", ["gene_set", "feature_id"], [["synthetic_response",g] for g in genes[1:6]]+[["synthetic_background",g] for g in genes[10:20]])
print(f"Wrote {len(cells)} observations and {len(counts)} sparse counts to {OUT}")
