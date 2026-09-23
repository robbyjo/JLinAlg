"""Regenerate selected Markdown-backed vignettes and xWAS validation. Requires Python Markdown.

Run from any directory: python docs/render-xwas-pages.py
The static output is checked in; website deployment needs no Python runtime.
"""
from pathlib import Path
import html
import re
import markdown

ROOT = Path(__file__).resolve().parents[1]
PAGES = {
    "single-cell": ("Replicated single-cell RNA analysis", "Sparse QC, pseudobulk, donor-aware state and abundance, functional scores and exploratory representations.", "CLI + Java/R · Source build"),
    "spatial-analysis": ("Physical spatial RNA analysis", "Tissue graphs, Moran and Geary tests, neighborhoods, sample-aware comparisons and distance gradients.", "CLI + Java · Source build"),
    "project-configuration": ("Project configuration", "Project-over-local YAML precedence and reproducible CLI runs.", "CLI · Source build"),
    "auditable-omics-pipeline": ("Auditable omics research pipeline", "A minimal study contract, project status, amendments, provenance, and future LLM-assisted planning.", "Roadmap design · Fixed pipeline first"),
    "variant-followup": ("Variant annotation and scoring", "Local snapshots, POST/API lookup, VEP/ANNOVAR adapters and transparent evidence integration.", "CLI · Source build"),
    "network-followup": ("Networks and candidate regulators", "Reference networks, WGCNA, sparse and differential associations, GENIE3 and CARNIVAL.", "CLI + Java/R · Source build"),
    "censored-ordinal": ("Censored and ordinal regression", "Fit Tobit, parametric survival-time and ordered logit/probit likelihoods with validated covariance.", "CLI + Java · Source build"),
    "sampling-models": ("Rare-event and survey inference", "Separate coefficient and prevalence corrections, and model sampling weights, strata and PSUs explicitly.", "CLI + Java · Source build"),
    "acat-rare-variants": ("ACAT rare-variant tests", "Run ACAT-V and canonical six-component ACAT-O from pooled scores or participant-level null models.", "CLI + Java · v0.3.6"),
    "enrichment": ("Gene-set enrichment", "Define analysis-specific backgrounds, select defensible tests, and interpret overlapping gene, disease and phenotype sets.", "CLI + Java · v0.3.6"),
    "ldsc": ("LDSC genetic architecture", "Estimate heritability and shared genetic architecture, then carry sampling uncertainty into the next model.", "CLI + Java · v0.3.6"),
    "predicted-omics": ("Genetically predicted TWAS and PWAS", "Connect molecular prediction weights to GWAS evidence with explicit allele, scale and LD alignment.", "CLI + Java · v0.3.6"),
    "genomic-factor": ("Shared genetic factors", "Model shared genetic variation and inspect common-factor SNP effects alongside heterogeneity.", "CLI + Java · v0.3.6"),
    "prediction-scores": ("Prediction scores", "Train portable scores and evaluate frozen predictions in an independent cohort.", "CLI + Java · v0.3.6"),
    "predictions-and-contrasts": ("Predictions, scenarios, and marginal effects", "Report expected responses, standardized scenario contrasts, risk ratios, and marginal effects with covariance-aware uncertainty.", "CLI + Java · v0.3.6"),
    "instrumental-variable-regression": ("Individual-level instrumental-variable regression", "Fit linear 2SLS models with robust or clustered inference, instrument-strength diagnostics, and explicit identification limits.", "CLI + Java · v0.3.6"),
    "conditional-score-conditioning": ("Summary-only conditional score analysis", "Import compatible cohort score blocks, condition by Schur complement, and pool auditable aggregate-data results.", "CLI + Java · v0.3.6"),
    "latent-confounders-and-batch": ("Latent confounders and batch effects", "Estimate unknown sample factors with PCA, SVA, AutoSVA, or PEER, and adjust known batches with ComBat.", "CLI + Java · v0.3.6"),
    "differential-regions-testing-imputation": ("Differential, regional, multiple-testing, and imputation workflows", "Run moderated continuous or count differential analysis, spatial EWAS regions, prespecified adaptive or hierarchical testing, and uncertainty-aware multiple imputation.", "CLI + Java · v0.3.6"),
}

def render(source, destination, title, description, root, eyebrow="CLI + Java · v0.3.6"):
    text = source.read_text(encoding="utf-8").split("\n", 1)[1]
    md = markdown.Markdown(extensions=["fenced_code", "tables", "toc"])
    body = md.convert(text)
    def links(match):
        url = match.group(1)
        if url.endswith("TODO.md"):
            url = "https://github.com/robbyjo/JLinAlg/blob/main/TODO.md"
        elif not url.startswith("http"):
            url = re.sub(r"(?i)CITATIONS\.md", "citations.html", url)
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
<span class="eyebrow">{html.escape(eyebrow)}</span><h1>{html.escape(title)}</h1><p>{html.escape(description)}</p></div></section>
<div class="container page-layout"><article class="prose">{body}</article>
<aside class="toc"><strong>On this page</strong>{toc}</aside></div>
</main><site-footer root="{root}"></site-footer></body></html>
'''
    destination.write_text(page, encoding="utf-8")

for slug, (title, description, eyebrow) in PAGES.items():
    render(ROOT / f"docs/vignettes/{slug}.md", ROOT / f"site/vignettes/{slug}.html", title, description, "../", eyebrow)
render(ROOT / "docs/xwas-followup-validation.md", ROOT / "site/xwas-followup-validation.html",
       "xWAS follow-up validation", "Independent numerical fixtures, reproducible checks and explicit estimator boundaries.", "./")
render(ROOT / "docs/latent-confounder-validation.md", ROOT / "site/latent-confounder-validation.html",
       "Latent-confounder validation and source audit", "Pinned upstream sources, numerical reference gates, performance design, and explicit implementation boundaries.", "./",
       "Source audit · Reproducible validation")
render(ROOT / "docs/regression-inference-validation.md", ROOT / "site/regression-inference-validation.html",
       "Regression inference validation", "Independent R fixtures, numerical checks, workload timings and estimator boundaries.", "./",
       "Source build · Reproducible validation")
render(ROOT / "docs/inference-workflows-validation.md", ROOT / "site/inference-workflows-validation.html",
       "Modern inference workflow validation", "Frozen Bioconductor comparisons, independent formula checks, reproducible commands, and explicit estimator boundaries.", "./",
       "Source audit · Reproducible validation")
render(ROOT / "docs/audit-fixes-validation.md", ROOT / "site/audit-fixes-validation.html",
       "Statistical audit fixes", "Regression checks, bootstrap imputation contracts and measured smoothing performance.", "./",
       "Source audit · Reproducible validation")
render(ROOT / "docs/cell-spatial-validation.md", ROOT / "site/cell-spatial-validation.html",
       "Single-cell and spatial validation", "Independent RNA and tissue workflow fixtures, scientific units and explicit remaining roadmap.", "./",
       "Source build · Reproducible validation")
