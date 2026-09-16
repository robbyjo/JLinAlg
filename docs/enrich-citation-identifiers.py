#!/usr/bin/env python3
"""Enrich docs/citations.json with PMID and PMCID values from NCBI."""

from __future__ import annotations

import argparse
import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

REGISTRY = Path(__file__).with_name("citations.json")
ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TOOL = "JLinAlgCitationIndexer"


def post(url: str, fields: dict[str, str]) -> bytes:
    request = urllib.request.Request(
        url,
        data=urllib.parse.urlencode(fields).encode("ascii"),
        headers={"User-Agent": f"{TOOL}/1.0"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read()


def doi_from_url(url: str) -> str | None:
    prefix = "https://doi.org/"
    return urllib.parse.unquote(url[len(prefix):]) if url.startswith(prefix) else None


def pubmed_ids(dois: list[str], email: str | None) -> list[str]:
    ids: set[str] = set()
    for offset in range(0, len(dois), 25):
        chunk = dois[offset:offset + 25]
        fields = {
            "db": "pubmed",
            "term": " OR ".join(f'"{doi}"[AID]' for doi in chunk),
            "retmax": "200",
            "retmode": "json",
            "tool": TOOL,
        }
        if email:
            fields["email"] = email
        result = json.loads(post(ESEARCH, fields))
        ids.update(result["esearchresult"]["idlist"])
        time.sleep(0.4)
    return sorted(ids, key=int)


def identifier_map(pmids: list[str], email: str | None) -> dict[str, dict[str, str]]:
    mapped: dict[str, dict[str, str]] = {}
    for offset in range(0, len(pmids), 100):
        fields = {
            "db": "pubmed",
            "id": ",".join(pmids[offset:offset + 100]),
            "retmode": "xml",
            "tool": TOOL,
        }
        if email:
            fields["email"] = email
        root = ET.fromstring(post(EFETCH, fields))
        for article in root.findall(".//PubmedArticle"):
            values = {
                node.attrib.get("IdType", "").lower(): (node.text or "").strip()
                for node in article.findall("./PubmedData/ArticleIdList/ArticleId")
            }
            doi = values.get("doi", "").lower()
            if doi and values.get("pubmed"):
                mapped[doi] = {
                    "pmid": values["pubmed"],
                    **({"pmcid": values["pmc"]} if values.get("pmc") else {}),
                }
        time.sleep(0.4)
    return mapped


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--email", help="contact email sent to NCBI E-utilities")
    args = parser.parse_args()

    data = json.loads(REGISTRY.read_text(encoding="utf-8"))
    dois = [doi for ref in data["references"] if (doi := doi_from_url(ref["url"]))]
    mapped = identifier_map(pubmed_ids(dois, args.email), args.email)
    enriched = 0
    pmcids = 0
    for ref in data["references"]:
        ref.pop("pmid", None)
        ref.pop("pmcid", None)
        doi = doi_from_url(ref["url"])
        identifiers = mapped.get(doi.lower(), {}) if doi else {}
        if identifiers:
            ref.update(identifiers)
            enriched += 1
            pmcids += int("pmcid" in identifiers)

    REGISTRY.write_text(
        json.dumps(data, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(f"Added PMID to {enriched} references and PMCID to {pmcids} references")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
