#!/usr/bin/env python3
"""Generate JLinAlg scientific-citation documentation from docs/citations.json."""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
from collections import OrderedDict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "docs" / "citations.json"
MARKDOWN_OUTPUT = ROOT / "docs" / "CITATIONS.md"
HTML_OUTPUT = ROOT / "site" / "citations.html"
HTML_START = "<!-- SCIENTIFIC-CITATIONS:START -->"
HTML_END = "<!-- SCIENTIFIC-CITATIONS:END -->"
ROOT_DOCS = {
    "additive-models.md": "additive-models",
    "beta-regression.md": "beta-regression",
    "gee.md": "predictions-and-contrasts",
    "grm-cli.md": "grm",
    "gwas-twas-pipeline.md": "association",
    "latent-confounder-validation.md": "latent-confounders-and-batch",
    "loess.md": "additive-models",
    "mr-timeseries-susie-sem.md": "advanced-extensions",
    "rare-variant-meta-analysis.md": "rare-variant-meta-analysis",
}


def load_registry() -> tuple[list[dict], dict[str, list[str]], dict[str, str]]:
    data = json.loads(REGISTRY.read_text(encoding="utf-8"))
    refs = data["references"]
    pages = data["pages"]
    aliases = data.get("markdownAliases", {})
    by_id = {ref["id"]: ref for ref in refs}
    if len(by_id) != len(refs):
        raise ValueError("Citation IDs must be unique")
    required = {"id", "category", "authors", "year", "title", "venue", "url", "features"}
    for ref in refs:
        missing = required - ref.keys()
        if missing:
            raise ValueError(f"{ref.get('id', '<unknown>')} lacks {sorted(missing)}")
        if not str(ref["url"]).startswith("https://"):
            raise ValueError(f"{ref['id']} must use an HTTPS source URL")
    for page, ids in pages.items():
        if not ids:
            raise ValueError(f"{page} has no citations")
        unknown = set(ids) - by_id.keys()
        if unknown:
            raise ValueError(f"{page} references unknown citations: {sorted(unknown)}")
    return refs, pages, aliases


def short_author(ref: dict) -> str:
    authors = ref["authors"]
    if " et al." in authors:
        return authors.split(" et al.", 1)[0].split()[-1] + " et al."
    first = re.split(r",| and ", authors, maxsplit=1)[0].strip()
    return first.split()[-1]


def render_markdown(refs: list[dict]) -> str:
    groups: OrderedDict[str, list[dict]] = OrderedDict()
    for ref in refs:
        groups.setdefault(ref["category"], []).append(ref)
    lines = [
        "# Scientific citations",
        "",
        "This bibliography links JLinAlg features to the primary publications that introduced or established their statistical methods. Cite both the relevant method paper and JLinAlg when reporting an analysis. Inclusion here documents methodological provenance; it is not a claim that every implementation detail is identical to the cited software.",
        "",
        "The website provides a [searchable citation index](https://robbyjo.github.io/JLinAlg/citations.html). This file and the per-vignette citation panels are generated from [citations.json](citations.json).",
        "",
    ]
    for category, entries in groups.items():
        lines.extend([f"## {category}", ""])
        for ref in entries:
            features = ", ".join(ref["features"])
            lines.extend([
                f'<a id="{ref["id"]}"></a>',
                f"### {ref['authors']} ({ref['year']})",
                "",
                f"**{ref['title']}.** {ref['venue']}. [Primary source]({ref['url']})",
                "",
                f"JLinAlg features: {features}.",
                "",
            ])
    return "\n".join(lines).rstrip() + "\n"


def render_html(refs: list[dict]) -> str:
    cards = []
    for ref in refs:
        search = " ".join([
            ref["category"], ref["authors"], str(ref["year"]), ref["title"],
            ref["venue"], *ref["features"], ref["url"],
        ]).lower()
        tags = "".join(f'<span>{html.escape(feature)}</span>' for feature in ref["features"])
        cards.append(f"""        <article class="citation-entry reveal" id="{ref['id']}" data-citation-entry data-search="{html.escape(search, quote=True)}">
          <div class="citation-meta"><span class="badge">{html.escape(ref['category'])}</span><span>{ref['year']}</span></div>
          <h2>{html.escape(ref['title'])}</h2>
          <p class="citation-authors">{html.escape(ref['authors'])}</p>
          <p>{html.escape(ref['venue'])}</p>
          <div class="citation-tags">{tags}</div>
          <a class="text-link" href="{html.escape(ref['url'], quote=True)}">Open primary source ↗</a>
        </article>""")
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="description" content="Search the primary scientific sources behind JLinAlg statistical methods.">
  <title>Scientific citations · JLinAlg</title>
  <link rel="stylesheet" href="assets/styles.css">
  <script src="assets/app.js" defer></script>
</head>
<body>
  <site-header root="./"></site-header>
  <main id="main">
    <section class="page-hero citation-hero">
      <div class="container narrow">
        <span class="eyebrow">Method provenance</span>
        <h1>Scientific citations</h1>
        <p>Search the primary publications behind JLinAlg methods by author, year, feature, title, journal, or DOI.</p>
        <label class="citation-search">
          <span>Search {len(refs)} sources</span>
          <input type="search" data-citation-search placeholder="Try MR-Egger, REML, SuSiE, 2015, or a DOI" autocomplete="off">
        </label>
        <p class="citation-result-summary" aria-live="polite"><strong data-citation-count>{len(refs)}</strong> sources shown</p>
      </div>
    </section>
    <section class="section">
      <div class="container citation-list" data-citation-list>
{chr(10).join(cards)}
        <p class="citation-empty" data-citation-empty hidden>No citations match that search. Try a method, author, year, or DOI.</p>
      </div>
    </section>
  </main>
  <site-footer root="./"></site-footer>
</body>
</html>
"""


def citation_block_html(slug: str, ids: list[str], by_id: dict[str, dict]) -> str:
    items = []
    for ref_id in ids:
        ref = by_id[ref_id]
        label = f"{short_author(ref)} ({ref['year']}) — {ref['features'][0]}"
        items.append(f'<li><a href="../citations.html#{ref_id}">{html.escape(label)}</a></li>')
    return f"""{HTML_START}
<section class="citation-strip" aria-labelledby="scientific-citations-{slug}">
  <div class="container narrow">
    <span class="eyebrow">Scientific foundations</span>
    <h2 id="scientific-citations-{slug}">Primary method citations</h2>
    <p>Cite the relevant method publication as well as JLinAlg when reporting results from this workflow.</p>
    <ul class="citation-links">
      {''.join(items)}
    </ul>
    <a class="text-link" href="../citations.html">Search all scientific citations →</a>
  </div>
</section>
{HTML_END}"""


def citation_block_markdown(ids: list[str], by_id: dict[str, dict], central: str) -> str:
    lines = [
        HTML_START,
        "## Scientific citations",
        "",
        "These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.",
        "",
    ]
    for ref_id in ids:
        ref = by_id[ref_id]
        lines.append(f"- [{ref['authors']} ({ref['year']}) — {ref['title']}]({central}#{ref_id})")
    lines.extend([
        "",
        "[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).",
        HTML_END,
    ])
    return "\n".join(lines)


def with_marked_block(text: str, block: str, before: str | None = None) -> str:
    pattern = re.compile(re.escape(HTML_START) + r".*?" + re.escape(HTML_END), re.DOTALL)
    if pattern.search(text):
        return pattern.sub(block, text)
    if before:
        index = text.rfind(before)
        if index < 0:
            raise ValueError(f"Cannot find insertion marker {before!r}")
        return text[:index].rstrip() + "\n\n" + block + "\n" + text[index:]
    return text.rstrip() + "\n\n" + block + "\n"


def expected_files(refs: list[dict], pages: dict[str, list[str]], aliases: dict[str, str]) -> dict[Path, str]:
    by_id = {ref["id"]: ref for ref in refs}
    expected = {
        MARKDOWN_OUTPUT: render_markdown(refs),
        HTML_OUTPUT: render_html(refs),
    }

    site_pages = sorted((ROOT / "site" / "vignettes").glob("*.html"))
    for path in site_pages:
        if path.stem == "index":
            continue
        if path.stem not in pages:
            raise ValueError(f"No citation mapping for {path.relative_to(ROOT)}")
        original = path.read_text(encoding="utf-8")
        block = citation_block_html(path.stem, pages[path.stem], by_id)
        expected[path] = with_marked_block(original, block, "</main>")

    markdown_pages = sorted((ROOT / "docs" / "vignettes").glob("*.md"))
    for path in markdown_pages:
        if path.name.lower() == "readme.md":
            continue
        key = aliases.get(path.stem, path.stem)
        if key not in pages:
            raise ValueError(f"No citation mapping for {path.relative_to(ROOT)}")
        original = path.read_text(encoding="utf-8")
        block = citation_block_markdown(pages[key], by_id, "../CITATIONS.md")
        expected[path] = with_marked_block(original, block)

    for filename, key in ROOT_DOCS.items():
        path = ROOT / "docs" / filename
        if not path.exists():
            continue
        original = path.read_text(encoding="utf-8")
        block = citation_block_markdown(pages[key], by_id, "CITATIONS.md")
        expected[path] = with_marked_block(original, block)
    return expected


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="fail when generated citation content is stale")
    args = parser.parse_args()
    refs, pages, aliases = load_registry()
    expected = expected_files(refs, pages, aliases)
    stale = []
    for path, content in expected.items():
        current = path.read_text(encoding="utf-8") if path.exists() else None
        if current != content:
            stale.append(path)
            if not args.check:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(content, encoding="utf-8", newline="\n")
    if args.check and stale:
        print("Citation outputs are stale:", file=sys.stderr)
        for path in stale:
            print(f"  {path.relative_to(ROOT)}", file=sys.stderr)
        print("Run: python docs/render-citations.py", file=sys.stderr)
        return 1
    verb = "validated" if args.check else "generated"
    print(f"{verb} {len(refs)} references across {len(expected) - 2} cited documents")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
