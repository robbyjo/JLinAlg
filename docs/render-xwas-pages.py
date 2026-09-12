"""Regenerate the four xWAS pages and validation page. Requires Python Markdown.

Run from any directory: python docs/render-xwas-pages.py
The static output is checked in; website deployment needs no Python runtime.
"""
from pathlib import Path
import html
import re
import markdown

ROOT = Path(__file__).resolve().parents[1]
PAGES = {
    "ldsc": ("LDSC genetic architecture", "Estimate heritability and shared genetic architecture, then carry sampling uncertainty into the next model."),
    "predicted-omics": ("Genetically predicted TWAS and PWAS", "Connect molecular prediction weights to GWAS evidence with explicit allele, scale and LD alignment."),
    "genomic-factor": ("Shared genetic factors", "Model shared genetic variation and inspect common-factor SNP effects alongside heterogeneity."),
    "prediction-scores": ("Prediction scores", "Train portable scores and evaluate frozen predictions in an independent cohort."),
}

def render(source, destination, title, description, root):
    text = source.read_text(encoding="utf-8").split("\n", 1)[1]
    md = markdown.Markdown(extensions=["fenced_code", "tables", "toc"])
    body = md.convert(text)
    def links(match):
        url = match.group(1)
        if url.endswith("TODO.md"):
            url = "https://github.com/robbyjo/JLinAlg/blob/main/TODO.md"
        elif not url.startswith("http"):
            url = url.replace(".md", ".html")
        return 'href="' + url + '"'
    body = re.sub(r'href="([^"]+)"', links, body)
    body = re.sub(r'<pre><code(?: class="language-([^\"]+)")?>([\s\S]*?)</code></pre>',
        lambda m: '<div class="code-window"><div class="code-bar"><span>'
        + html.escape(m.group(1) or "Input / output") + '</span></div><pre><code>'
        + m.group(2) + '</code></pre></div>', body)
    toc = ''.join(f'<a href="#{item["id"]}">{html.escape(item["name"])}</a>' for item in md.toc_tokens)
    page = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="description" content="{html.escape(description)}"><title>{html.escape(title)} · JLinAlg</title>
<link rel="stylesheet" href="{root}assets/styles.css"><script src="{root}assets/app.js" defer></script></head>
<body><site-header root="{root}"></site-header><main id="main">
<section class="page-hero"><div class="container"><div class="breadcrumbs"><a href="{root}index.html">Home</a><span>/</span><a href="{root}vignettes/index.html">Vignettes</a><span>/</span><span>{html.escape(title)}</span></div>
<span class="eyebrow">CLI + Java · Source build</span><h1>{html.escape(title)}</h1><p>{html.escape(description)}</p></div></section>
<div class="container page-layout"><article class="prose">{body}</article>
<aside class="toc"><strong>On this page</strong>{toc}</aside></div>
</main><site-footer root="{root}"></site-footer></body></html>
'''
    destination.write_text(page, encoding="utf-8")

for slug, (title, description) in PAGES.items():
    render(ROOT / f"docs/vignettes/{slug}.md", ROOT / f"site/vignettes/{slug}.html", title, description, "../")
render(ROOT / "docs/xwas-followup-validation.md", ROOT / "site/xwas-followup-validation.html",
       "xWAS follow-up validation", "Independent numerical fixtures, reproducible checks and explicit estimator boundaries.", "./")
