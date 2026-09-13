"""Independent 90-digit full-distribution Wallenius reference, standard library only.

Unlike the Java absorbing-tail calculation, propagate every possible count through
all draws, then sum the final distribution. No probability subtraction is used.
Run from the repository root; the output is consumed by offline Java tests.
"""
from decimal import Decimal, localcontext
from pathlib import Path

cases = [(100, 900, 100, 25, "0.5"), (100, 100, 100, 90, "0.1"),
         (200, 800, 100, 45, "2"), (30, 70, 95, 28, "3")]
rows = ["red\twhite\tdraws\tobserved\todds\tp"]
with localcontext() as context:
    context.prec = 90
    for red, white, draws, observed, text_odds in cases:
        odds = Decimal(text_odds)
        distribution = {0: Decimal(1)}
        for t in range(draws):
            following = {}
            for x, mass in distribution.items():
                r, w = red - x, white - t + x
                success = r * odds / (r * odds + w)
                if r:
                    following[x + 1] = following.get(x + 1, Decimal(0)) + mass * success
                if w:
                    following[x] = following.get(x, Decimal(0)) + mass * (1 - success)
            distribution = following
        assert abs(sum(distribution.values()) - 1) < Decimal("1e-85")
        tail = sum(p for x, p in distribution.items() if x >= observed)
        rows.append("\t".join(map(str, (red, white, draws, observed, text_odds, tail))))
Path("src/test/resources/enrichment/wallenius-precision.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
