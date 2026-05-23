"""
=============================================================================
  tier1_citation_audit.py
=============================================================================

  Single-file analysis pipeline for:

  Miller RC, Wrightson T. The post-ChatGPT citation fabrication surge is
  not universal across medical specialties. The Lancet (Correspondence,
  submitted 2026).
  Zenodo DOI: 10.5281/zenodo.20350450

  In reply to: Topaz M, Roguin N, Gupta P, et al. Fabricated citations: an
  audit across 2.5 million biomedical papers. Lancet 2026; 407: 1779-81.


  WHAT THIS FILE DOES
  -------------------
  Audits citation integrity across the seven fully open-access, PMC-indexed
  radiation-oncology journals, applying the Topaz et al. analysis design
  (Categories 1, 2, and 3 strong fabrication candidates).

  Three pipeline stages:
      fetch    PMC esearch + efetch + JATS parsing  -> data/<key>_refs_raw.csv
      verify   DOI -> PMID -> title-search cascade   -> data/<key>_verified.csv
      analyze  Topaz-style analysis                   -> analysis_report.md


  REQUIRED CREDENTIALS
  --------------------
  No credentials are stored in this file. The two required values are read
  from environment variables at startup:

      NCBI_EMAIL    YOUR institutional or personal email address.
                    Required. NCBI may throttle or block anonymous traffic.
                    Example:  jane.smith@hospital.org

      NCBI_API_KEY  YOUR free NCBI API key. Optional but recommended -
                    raises the rate limit from 3 to 10 req/s, cutting a
                    full seven-journal run from ~60 h to ~20 h.
                    Get one in 1 minute at:
                    https://account.ncbi.nlm.nih.gov/settings/

  Set them in your shell BEFORE running this script:

      export NCBI_EMAIL="your.email@institution.edu"
      export NCBI_API_KEY="your_api_key_here"        # optional

  Or put those two lines in a config.sh file (gitignored) and source it:

      source config.sh
      python tier1_citation_audit.py

  This pipeline does NOT require an ORCID login, password, or token.
  ORCID identifiers appear only in the manuscript (they are public IDs,
  not authentication credentials) and are not used anywhere in this code.


  REQUIREMENTS
  ------------
  Python 3.9+ and:
      pip install requests tenacity lxml rapidfuzz pandas scipy tqdm matplotlib


  USAGE
  -----
      python tier1_citation_audit.py                              # all 7 journals, all stages
      python tier1_citation_audit.py --journal adro               # one journal
      python tier1_citation_audit.py --stage fetch                # fetch only
      python tier1_citation_audit.py --stage verify               # re-verify (no re-fetch)
      python tier1_citation_audit.py --stage analyze              # analysis only
      python tier1_citation_audit.py --n 25                       # pilot: 25 articles per journal


  LICENSE
  -------
  MIT
=============================================================================
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import random
import sys
import time
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Iterator, Optional

import numpy as np
import pandas as pd
import requests
from lxml import etree
from rapidfuzz import fuzz
from scipy.optimize import brentq
from scipy.stats import fisher_exact, norm
from tenacity import retry, stop_after_attempt, wait_exponential
from tqdm import tqdm


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
log = logging.getLogger("audit")

ROOT = Path(__file__).parent
DATA = ROOT / "data"


# =============================================================================
# SECTION 1 - NCBI E-UTILITIES ACCESS
# -----------------------------------------------------------------------------
# Credentials are loaded from environment variables only. Never put your
# email or API key directly in this file.
# =============================================================================

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

NCBI_EMAIL: str = os.environ.get("NCBI_EMAIL", "").strip()
NCBI_API_KEY: str = os.environ.get("NCBI_API_KEY", "").strip()

# Throttle interval: 10 req/s with an API key, 3 req/s without (a small
# margin is added so bursts never trip NCBI's limiter).
_MIN_INTERVAL = 0.11 if NCBI_API_KEY else 0.34
_last_call = [0.0]


def environment_ok() -> tuple[bool, list[str]]:
    """Return (ok, warnings). ok is False only for fatal misconfiguration."""
    warnings: list[str] = []
    ok = True
    if not NCBI_EMAIL:
        warnings.append(
            "NCBI_EMAIL is not set.\n"
            "  Set it in your shell: export NCBI_EMAIL='your.email@institution.edu'\n"
            "  NCBI may throttle or block anonymous traffic."
        )
        ok = False
    if not NCBI_API_KEY:
        warnings.append(
            "NCBI_API_KEY is not set. Pipeline will run at 3 req/s (~3x slower).\n"
            "  Get a free key at: https://account.ncbi.nlm.nih.gov/settings/"
        )
    return ok, warnings


def _throttle() -> None:
    elapsed = time.time() - _last_call[0]
    if elapsed < _MIN_INTERVAL:
        time.sleep(_MIN_INTERVAL - elapsed)
    _last_call[0] = time.time()


def eutils_params(tool: str, extra: dict) -> dict:
    """Build an E-utilities parameter dict with auth fields attached."""
    p: dict = {"tool": tool, "email": NCBI_EMAIL or "anonymous@example.org"}
    if NCBI_API_KEY:
        p["api_key"] = NCBI_API_KEY
    p.update(extra)
    return p


@retry(stop=stop_after_attempt(4), wait=wait_exponential(multiplier=1, min=2, max=30))
def eutils_get(endpoint: str, query: dict) -> requests.Response:
    """GET an E-utilities endpoint (e.g. 'esearch.fcgi') with throttle + retry."""
    _throttle()
    r = requests.get(f"{EUTILS}/{endpoint}", params=query, timeout=60)
    r.raise_for_status()
    return r


# =============================================================================
# SECTION 2 - JOURNAL CONFIGURATION
# -----------------------------------------------------------------------------
# The seven Tier 1 open-access radiation-oncology journals indexed in PMC,
# with verified NLM IDs and MEDLINE abbreviations.
# =============================================================================

JOURNALS_CONFIG = [
    {"key": "adro",   "name": "Advances in Radiation Oncology",
     "medline": "Adv Radiat Oncol",                "issn_e": "2452-1094",
     "nlm_id": "101677247", "doi_prefix": "10.1016/j.adro",
     "publisher": "ASTRO / Elsevier",              "first_year": 2016},

    {"key": "ctro",   "name": "Clinical and Translational Radiation Oncology",
     "medline": "Clin Transl Radiat Oncol",        "issn_e": "2405-6308",
     "nlm_id": "101713416", "doi_prefix": "10.1016/j.ctro",
     "publisher": "ESTRO / Elsevier",              "first_year": 2016},

    {"key": "radonc", "name": "Radiation Oncology",
     "medline": "Radiat Oncol",                    "issn_e": "1748-717X",
     "nlm_id": "101265111", "doi_prefix": "10.1186/s13014",
     "publisher": "BMC / Springer Nature",         "first_year": 2006},

    {"key": "phro",   "name": "Physics and Imaging in Radiation Oncology",
     "medline": "Phys Imaging Radiat Oncol",       "issn_e": "2405-6316",
     "nlm_id": "101704276", "doi_prefix": "10.1016/j.phro",
     "publisher": "ESTRO / Elsevier",              "first_year": 2017},

    {"key": "tipsro", "name": "Technical Innovations and Patient Support in Radiation Oncology",
     "medline": "Tech Innov Patient Support Radiat Oncol", "issn_e": "2405-6324",
     "nlm_id": "101762366", "doi_prefix": "10.1016/j.tipsro",
     "publisher": "ESTRO / Elsevier",              "first_year": 2017},

    {"key": "jrr",    "name": "Journal of Radiation Research",
     "medline": "J Radiat Res",                    "issn_e": "1349-9157",
     "nlm_id": "0376611",   "doi_prefix": "10.1093/jrr",
     "publisher": "Oxford UP / JRRS + JASTRO",     "first_year": 2000},

    {"key": "rpor",   "name": "Reports of Practical Oncology and Radiotherapy",
     "medline": "Rep Pract Oncol Radiother",       "issn_e": "2083-4640",
     "nlm_id": "100885761", "doi_prefix": "10.1016/j.rpor",
     "publisher": "Polish Society / Via Medica",   "first_year": 2009},
]


@dataclass
class JournalConfig:
    key: str
    name: str
    medline: str
    issn_e: str
    nlm_id: str
    doi_prefix: str
    publisher: str
    first_year: int

    @property
    def journal_term(self) -> str:
        return f'"{self.medline}"[Journal]'

    @property
    def tool_name(self) -> str:
        return f"tier1-citation-audit-{self.key}"


def load_journals() -> dict[str, JournalConfig]:
    return {j["key"]: JournalConfig(**j) for j in JOURNALS_CONFIG}


# =============================================================================
# SECTION 3 - REFERENCE RECORD AND PMID HYGIENE
# =============================================================================

@dataclass
class Reference:
    journal_key: str
    article_pmcid: str
    ref_id: str
    publication_type: Optional[str] = None
    first_author_surname: Optional[str] = None
    title: Optional[str] = None
    source: Optional[str] = None
    year: Optional[int] = None
    volume: Optional[str] = None
    issue: Optional[str] = None
    fpage: Optional[str] = None
    lpage: Optional[str] = None
    doi: Optional[str] = None
    pmid: Optional[str] = None
    pmcid: Optional[str] = None


RAW_FIELDS = [f.name for f in fields(Reference)]


def clean_pmid(raw) -> Optional[str]:
    """Reduce a PMID to bare digits. Guards against '25705639.0' and stray
    whitespace so PubMed lookups never fail on a malformed identifier."""
    if raw is None:
        return None
    s = str(raw).strip()
    if s in ("", "nan", "None", "<NA>"):
        return None
    if s.endswith(".0"):
        s = s[:-2]
    return s if s.isdigit() else None


# =============================================================================
# SECTION 4 - PMC FETCH AND JATS PARSING
# =============================================================================

def search_journal_pmcids(journal: JournalConfig, retmax: int = 100000) -> list[str]:
    """Return every PMCID in the journal, paginating through esearch."""
    pmcids: list[str] = []
    retstart = 0
    batch = 500
    while True:
        resp = eutils_get(
            "esearch.fcgi",
            eutils_params(journal.tool_name, {
                "db": "pmc",
                "term": journal.journal_term,
                "retmax": batch,
                "retstart": retstart,
                "retmode": "json",
            }),
        ).json()
        result = resp.get("esearchresult", {})
        ids = result.get("idlist", [])
        if not ids:
            break
        pmcids.extend(ids)
        total = int(result.get("count", 0))
        retstart += batch
        if retstart >= total or retstart >= retmax:
            break
    return [i if i.startswith("PMC") else f"PMC{i}" for i in pmcids]


def fetch_article_xml(pmcid: str, journal: JournalConfig) -> bytes:
    """Download the full JATS XML for one article."""
    bare = pmcid.replace("PMC", "")
    resp = eutils_get(
        "efetch.fcgi",
        eutils_params(journal.tool_name, {"db": "pmc", "id": bare, "rettype": "xml"}),
    )
    return resp.content


def _text(elem) -> Optional[str]:
    if elem is None:
        return None
    s = "".join(elem.itertext()).strip()
    return s or None


def parse_references(journal_key: str, pmcid: str, xml_bytes: bytes) -> list[Reference]:
    """Parse JATS XML into Reference records (one per <ref>)."""
    parser = etree.XMLParser(recover=True, huge_tree=True)
    try:
        root = etree.fromstring(xml_bytes, parser=parser)
    except etree.XMLSyntaxError as e:
        log.warning("XML parse failed for %s: %s", pmcid, e)
        return []
    if root is None:
        return []

    refs: list[Reference] = []
    for ref in root.iter("ref"):
        cit = ref.find(".//element-citation") or ref.find(".//mixed-citation")
        if cit is None:
            continue

        r = Reference(
            journal_key=journal_key,
            article_pmcid=pmcid,
            ref_id=ref.get("id") or "",
        )
        r.publication_type = cit.get("publication-type")

        authors: list[str] = []
        for name in cit.findall(".//person-group[@person-group-type='author']/name"):
            surname = _text(name.find("surname"))
            if surname:
                authors.append(surname)
        if not authors:
            for name in cit.findall(".//name"):
                surname = _text(name.find("surname"))
                if surname:
                    authors.append(surname)
        if authors:
            r.first_author_surname = authors[0].split()[0]

        r.title = _text(cit.find("article-title")) or _text(cit.find("chapter-title"))
        r.source = _text(cit.find("source"))
        year_text = _text(cit.find("year"))
        if year_text and year_text.isdigit():
            r.year = int(year_text)
        r.volume = _text(cit.find("volume"))
        r.issue = _text(cit.find("issue"))
        r.fpage = _text(cit.find("fpage"))
        r.lpage = _text(cit.find("lpage"))

        for pid in cit.findall("pub-id"):
            kind = pid.get("pub-id-type", "")
            val = (_text(pid) or "").strip()
            if kind == "doi" and val:
                r.doi = val.lower()
            elif kind == "pmid":
                r.pmid = clean_pmid(val)
            elif kind == "pmcid" and val:
                r.pmcid = val

        refs.append(r)
    return refs


def collect_journal_references(
    journal: JournalConfig,
    pmcids: list[str],
) -> Iterator[Reference]:
    """Fetch and parse the given articles, yielding Reference records.
    Failures on a single article are logged and skipped, not fatal."""
    for pmcid in pmcids:
        try:
            xml = fetch_article_xml(pmcid, journal)
            yield from parse_references(journal.key, pmcid, xml)
        except Exception as e:  # noqa: BLE001
            log.warning("fetch/parse failed for %s (%s): %s", pmcid, journal.key, e)


# =============================================================================
# SECTION 5 - REFERENCE VERIFICATION CASCADE
# -----------------------------------------------------------------------------
# Three fixes baked in:
#   FIX 1  cascade order: DOI -> PMID -> title search (v1 skipped PMID
#          after a DOI failure)
#   FIX 2  PMIDs sanitised to bare digits before every PubMed call
#   FIX 3  title similarity >= 85 confirms the work alone (aligns with
#          Topaz Category 2: only sim < 50 counts as substantially different)
# =============================================================================

CROSSREF = "https://api.crossref.org"
USER_AGENT = f"tier1-citation-audit/3.0 (mailto:{NCBI_EMAIL or 'anonymous@example.org'})"

TITLE_STRONG = 85     # title similarity at/above this confirms the work alone
TITLE_THRESHOLD = 80  # minimum title similarity when corroboration is also required
YEAR_TOLERANCE = 1    # permitted |cited year - record year| in corroborated path


@dataclass
class Verdict:
    status: str
    method: str
    matched_doi: Optional[str] = None
    matched_pmid: Optional[str] = None
    matched_title: Optional[str] = None
    matched_year: Optional[int] = None
    title_similarity: Optional[float] = None
    notes: Optional[str] = None


VERDICT_FIELDS = [
    "verify_status", "verify_method", "matched_doi", "matched_pmid",
    "matched_title", "matched_year", "title_similarity", "notes",
]


def verdict_to_row(v: Verdict) -> dict:
    return {
        "verify_status": v.status,
        "verify_method": v.method,
        "matched_doi": v.matched_doi,
        "matched_pmid": v.matched_pmid,
        "matched_title": v.matched_title,
        "matched_year": v.matched_year,
        "title_similarity": v.title_similarity,
        "notes": v.notes,
    }


def _names_match(claimed_surname: Optional[str], candidate_authors: list[str]) -> bool:
    if not claimed_surname or not candidate_authors:
        return True
    cs = claimed_surname.lower()
    return any(cs in a.lower() for a in candidate_authors)


def _metadata_match(
    ref: Reference,
    cand_title: str,
    cand_year: Optional[int],
    cand_authors: list[str],
) -> tuple[bool, float]:
    """Return (is_match, title_similarity).
    FIX 3: a strong title match alone confirms the work."""
    sim = fuzz.token_set_ratio((ref.title or "").lower(), cand_title.lower())
    if ref.title and sim >= TITLE_STRONG:
        return (True, sim)
    title_ok = sim >= TITLE_THRESHOLD if ref.title else True
    year_ok = (
        ref.year is None
        or cand_year is None
        or abs(int(cand_year) - int(ref.year)) <= YEAR_TOLERANCE
    )
    author_ok = (
        _names_match(ref.first_author_surname, cand_authors)
        if ref.first_author_surname else True
    )
    return (title_ok and year_ok and author_ok, sim)


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=20))
def _crossref_get(path: str, query: Optional[dict] = None) -> dict:
    r = requests.get(
        f"{CROSSREF}{path}", params=query,
        headers={"User-Agent": USER_AGENT}, timeout=30,
    )
    if r.status_code == 404:
        return {}
    r.raise_for_status()
    return r.json()


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=20))
def _pubmed_summary(pmid: str) -> dict:
    time.sleep(0.05)  # gentle extra spacing on top of the eutils throttle
    r = eutils_get(
        "esummary.fcgi",
        eutils_params("tier1-citation-audit-verify",
                      {"db": "pubmed", "id": pmid, "retmode": "json"}),
    )
    return r.json()


def verify_doi(ref: Reference) -> Optional[Verdict]:
    if not ref.doi:
        return None
    data = _crossref_get(f"/works/{ref.doi}")
    if not data or "message" not in data:
        return Verdict(
            status="not_found", method="doi",
            notes=f"DOI {ref.doi} did not resolve at Crossref",
        )
    msg = data["message"]
    cand_title = (msg.get("title") or [""])[0]
    parts = msg.get("issued", {}).get("date-parts") or [[]]
    cand_year = parts[0][0] if parts and parts[0] else None
    cand_authors = [a.get("family", "") for a in msg.get("author", [])]
    ok, sim = _metadata_match(ref, cand_title, cand_year, cand_authors)
    return Verdict(
        status="verified" if ok else "metadata_mismatch",
        method="doi",
        matched_doi=msg.get("DOI"),
        matched_title=cand_title,
        matched_year=cand_year,
        title_similarity=sim,
        notes=None if ok else "DOI resolved but cited metadata does not match record",
    )


def verify_pmid(ref: Reference) -> Optional[Verdict]:
    pmid = clean_pmid(ref.pmid)
    if not pmid:
        return None
    resp = _pubmed_summary(pmid)
    rec = resp.get("result", {}).get(pmid)
    if not rec or rec.get("error"):
        return Verdict(
            status="not_found", method="pmid",
            notes=f"PMID {pmid} not found in PubMed",
        )
    cand_title = rec.get("title", "")
    pubdate = rec.get("pubdate", "")
    cand_year = int(pubdate[:4]) if pubdate[:4].isdigit() else None
    cand_authors = [a.get("name", "").split()[0] for a in rec.get("authors", [])]
    ok, sim = _metadata_match(ref, cand_title, cand_year, cand_authors)
    return Verdict(
        status="verified" if ok else "metadata_mismatch",
        method="pmid",
        matched_pmid=pmid,
        matched_title=cand_title,
        matched_year=cand_year,
        title_similarity=sim,
        notes=None if ok else "PMID resolved but cited metadata does not match record",
    )


def verify_by_search(ref: Reference) -> Verdict:
    """Last resort: Crossref bibliographic search on title + author + year."""
    if not ref.title or not ref.first_author_surname:
        return Verdict(
            status="no_identifiers", method="none",
            notes="insufficient metadata for search-based verification",
        )
    query: dict = {
        "query.bibliographic": ref.title,
        "query.author": ref.first_author_surname,
        "rows": 5,
    }
    if ref.year:
        query["filter"] = f"from-pub-date:{ref.year - 1},until-pub-date:{ref.year + 1}"
    data = _crossref_get("/works", query=query)
    items = data.get("message", {}).get("items", [])
    best, best_sim = None, 0.0
    for it in items:
        cand_title = (it.get("title") or [""])[0]
        parts = it.get("issued", {}).get("date-parts") or [[]]
        cand_year = parts[0][0] if parts and parts[0] else None
        cand_authors = [a.get("family", "") for a in it.get("author", [])]
        ok, sim = _metadata_match(ref, cand_title, cand_year, cand_authors)
        if ok and sim > best_sim:
            best, best_sim = it, sim
    if best is None:
        return Verdict(
            status="not_found", method="title_search",
            notes="no Crossref result matched the cited metadata",
        )
    parts = best.get("issued", {}).get("date-parts") or [[None]]
    return Verdict(
        status="verified",
        method="title_search",
        matched_doi=best.get("DOI"),
        matched_title=(best.get("title") or [""])[0],
        matched_year=parts[0][0] if parts and parts[0] else None,
        title_similarity=best_sim,
    )


def verify(ref: Reference) -> Verdict:
    """FIX 1: DOI -> PMID -> title search, in that order."""
    try:
        v_doi = verify_doi(ref)
        if v_doi is not None and v_doi.status == "verified":
            return v_doi

        v_pmid = verify_pmid(ref)
        if v_pmid is not None and v_pmid.status == "verified":
            if v_doi is not None:
                note = f"DOI ({v_doi.method}={v_doi.status}) overridden by PMID"
                v_pmid.notes = f"{v_pmid.notes}; {note}" if v_pmid.notes else note
            return v_pmid

        v_search = verify_by_search(ref)
        if v_search.status == "verified":
            failed = [v for v in (v_doi, v_pmid) if v is not None]
            if failed:
                tags = "; ".join(f"{v.method}={v.status}" for v in failed)
                base = v_search.notes + "; " if v_search.notes else ""
                v_search.notes = f"{base}identifier(s) failed: {tags}"
            return v_search

        for v in (v_doi, v_pmid):
            if v is not None:
                return v
        return v_search
    except Exception as e:  # noqa: BLE001
        log.warning("verify error for %s/%s: %s", ref.article_pmcid, ref.ref_id, e)
        return Verdict(status="error", method="none", notes=str(e))


# =============================================================================
# SECTION 6 - TOPAZ-STYLE STATISTICAL ANALYSIS
# -----------------------------------------------------------------------------
# Methodological decisions baked in:
#   - restricted to journal-type references
#   - article publication year proxied by max reference year (clipped 1990-2026)
#   - Pre/post-ChatGPT boundary: proxy year <= 2022 vs >= 2023
#   - Strong fabrication candidate = Topaz Cat 1, 2, or 3:
#       Cat 1  not_found with DOI present
#       Cat 2  metadata_mismatch with title_similarity < 50
#       Cat 3  not_found with neither DOI nor PMID
#   - PRIMARY analysis is article level; reference-level is sensitivity
#   - >=2018 sensitivity is PER REFERENCE, never per article (per-article form
#     is confounded by recent articles carrying more recent references)
# =============================================================================

def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def diff_ci(p1, n1, p2, n2, z: float = 1.96) -> tuple[float, float]:
    se = np.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    d = p2 - p1
    return (d - z * se, d + z * se)


def logit_year(years, outcomes, iters: int = 100):
    """Single-predictor logistic regression via IRLS. Predictor standardised
    internally. Returns (OR_per_year, ci_lo, ci_hi, p_value)."""
    x = np.asarray(years, float)
    y = np.asarray(outcomes, float)
    sd = x.std()
    if sd == 0 or len(np.unique(y)) < 2:
        return (float("nan"),) * 4
    X = np.column_stack([np.ones_like(x), (x - x.mean()) / sd])
    b = np.zeros(2)
    for _ in range(iters):
        eta = X @ b
        mu = 1 / (1 + np.exp(-eta))
        W = np.clip(mu * (1 - mu), 1e-9, None)
        z = eta + (y - mu) / W
        XtW = X.T * W
        b_new = np.linalg.solve(XtW @ X, XtW @ z)
        if np.max(np.abs(b_new - b)) < 1e-10:
            b = b_new
            break
        b = b_new
    eta = X @ b
    mu = 1 / (1 + np.exp(-eta))
    W = np.clip(mu * (1 - mu), 1e-9, None)
    cov = np.linalg.inv((X.T * W) @ X)
    slope = b[1] / sd
    se = np.sqrt(cov[1, 1]) / sd
    p = 2 * (1 - norm.cdf(abs(slope / se)))
    return (np.exp(slope), np.exp(slope - 1.96 * se), np.exp(slope + 1.96 * se), p)


def min_detectable_rate(p0: float, n1: int, n2: int, power: float = 0.80) -> float:
    """Smallest post-period rate detectable at the given power, alpha=0.05."""
    za = norm.ppf(0.975)
    zb_target = power

    def gap(p1):
        pbar = (p0 * n1 + p1 * n2) / (n1 + n2)
        se0 = np.sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n2))
        se1 = np.sqrt(p0 * (1 - p0) / n1 + p1 * (1 - p1) / n2)
        return norm.cdf((abs(p1 - p0) - za * se0) / se1) - zb_target

    try:
        return brentq(gap, p0 + 1e-6, 0.999)
    except ValueError:
        return float("nan")


def load_references(paths: list[Path]) -> pd.DataFrame:
    """Load verified CSV(s), concatenate, restrict to journal-type references,
    and attach Topaz category flags."""
    frames = []
    for p in paths:
        df = pd.read_csv(p, low_memory=False, dtype={"pmid": "Int64"})
        if "journal_key" not in df.columns:
            df["journal_key"] = p.stem.replace("_verified", "")
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)

    df["has_doi"] = df["doi"].notna()
    df["has_pmid"] = df["pmid"].notna()

    jr = df[df["publication_type"] == "journal"].copy()
    jr["cat1"] = (jr["verify_status"] == "not_found") & jr["has_doi"]
    jr["cat2"] = (jr["verify_status"] == "metadata_mismatch") & (jr["title_similarity"] < 50)
    jr["cat3"] = (jr["verify_status"] == "not_found") & (~jr["has_doi"]) & (~jr["has_pmid"])
    jr["candidate"] = jr["cat1"] | jr["cat2"] | jr["cat3"]
    jr["year_clip"] = jr["year"].where((jr["year"] >= 1990) & (jr["year"] <= 2026))
    return jr


def article_table(jr: pd.DataFrame) -> pd.DataFrame:
    g = jr.groupby("article_pmcid")
    art = pd.DataFrame({
        "journal_key": g["journal_key"].first(),
        "proxy_year": g["year_clip"].max(),
        "n_refs": g.size(),
        "has_candidate": g["candidate"].any(),
    }).dropna(subset=["proxy_year"])
    art["proxy_year"] = art["proxy_year"].astype(int)
    return art


class Report:
    def __init__(self) -> None:
        self.lines: list[str] = []

    def w(self, line: str = "") -> None:
        print(line)
        self.lines.append(line)

    def text(self) -> str:
        return "\n".join(self.lines)


def _window_block(rep: Report, jr: pd.DataFrame, art: pd.DataFrame,
                  cutoff: int, label: str) -> None:
    a = art[art["proxy_year"] >= cutoff].copy()
    a["era"] = np.where(a["proxy_year"] <= 2022, "pre", "post")
    pre, post = a[a["era"] == "pre"], a[a["era"] == "post"]
    if len(pre) == 0 or len(post) == 0:
        rep.w(f"  {label}: insufficient data in one era; skipped")
        rep.w("")
        return

    # Article-level primary
    kp, np_, kq, nq = (pre["has_candidate"].sum(), len(pre),
                       post["has_candidate"].sum(), len(post))
    pp, pq = kp / np_, kq / nq
    OR, pv = fisher_exact([[kp, np_ - kp], [kq, nq - kq]])
    dlo, dhi = diff_ci(pp, np_, pq, nq)
    lop, hip = wilson(kp, np_)
    loq, hiq = wilson(kq, nq)

    # Reference-level sensitivity
    j = jr.dropna(subset=["year_clip"]).merge(
        a[["era"]], left_on="article_pmcid", right_index=True, how="inner")
    rpre, rpost = j[j["era"] == "pre"]["candidate"], j[j["era"] == "post"]["candidate"]
    rkp, rnp, rkq, rnq = rpre.sum(), len(rpre), rpost.sum(), len(rpost)
    rOR, rpv = fisher_exact([[rkp, rnp - rkp], [rkq, rnq - rkq]])

    # Logistic trend
    lOR, llo, lhi, lp = logit_year(a["proxy_year"], a["has_candidate"].astype(int))

    # Power
    mdr = min_detectable_rate(pp, np_, nq)
    fold = mdr / pp if pp > 0 else float("nan")

    rep.w(f"  {label}")
    rep.w(f"    PRIMARY (article level)")
    rep.w(f"      Pre  : {kp:>6}/{np_:<7} = {100*pp:6.2f}%  (95% CI {100*lop:.2f}-{100*hip:.2f})")
    rep.w(f"      Post : {kq:>6}/{nq:<7} = {100*pq:6.2f}%  (95% CI {100*loq:.2f}-{100*hiq:.2f})")
    rep.w(f"      Diff : {100*(pq-pp):+.2f} pp  (95% CI {100*dlo:+.2f} to {100*dhi:+.2f})")
    rep.w(f"      Fisher exact: OR = {OR:.3f},  p = {pv:.4f}")
    rep.w(f"    SENSITIVITY (reference level)")
    rep.w(f"      Pre {100*rkp/rnp:.3f}%   Post {100*rkq/rnq:.3f}%   OR = {rOR:.3f},  p = {rpv:.4f}")
    rep.w(f"    SENSITIVITY (logistic, continuous year)")
    rep.w(f"      OR/year = {lOR:.4f}  (95% CI {llo:.4f}-{lhi:.4f}),  p = {lp:.4f}")
    rep.w(f"    POWER: min detectable post-rate {100*mdr:.2f}%  =  {fold:.2f}x fold-change")
    rep.w("")


def _sensitivity_2018(rep: Report, jr: pd.DataFrame, art: pd.DataFrame) -> None:
    """Per-reference candidate rate among references dated >= 2018."""
    a = art.copy()
    a["era"] = np.where(a["proxy_year"] <= 2022, "pre", "post")
    j18 = jr[jr["year_clip"] >= 2018].merge(
        a[["era"]], left_on="article_pmcid", right_index=True, how="inner")
    pre, post = j18[j18["era"] == "pre"]["candidate"], j18[j18["era"] == "post"]["candidate"]
    if len(pre) == 0 or len(post) == 0:
        rep.w("  References >= 2018: insufficient data; skipped")
        rep.w("")
        return
    kp, np_, kq, nq = pre.sum(), len(pre), post.sum(), len(post)
    OR, pv = fisher_exact([[kp, np_ - kp], [kq, nq - kq]])
    rep.w("  SENSITIVITY -- references dated >= 2018 (per reference, full corpus)")
    rep.w(f"    Pre  refs >=2018: {kp:>5}/{np_:<7} = {100*kp/np_:.3f}%")
    rep.w(f"    Post refs >=2018: {kq:>5}/{nq:<7} = {100*kq/nq:.3f}%")
    rep.w(f"    Fisher exact: OR = {OR:.3f},  p = {pv:.4f}")
    rep.w("")


def run_analysis(paths: list[Path]) -> str:
    jr = load_references(paths)
    art = article_table(jr)
    rep = Report()

    n_journals = jr["journal_key"].nunique()
    scope = f"{n_journals} journals pooled" if n_journals > 1 else str(jr["journal_key"].iloc[0])

    rep.w("=" * 76)
    rep.w(f"  TOPAZ-STYLE CITATION-INTEGRITY ANALYSIS  --  {scope}")
    rep.w("=" * 76)
    rep.w(f"  Journal-type references : {len(jr):,}")
    rep.w(f"  Articles                : {jr['article_pmcid'].nunique():,}")
    rep.w("")

    if n_journals > 1:
        rep.w("  Per-journal contribution:")
        rep.w(f"  {'journal':<10}{'articles':>10}{'refs':>11}{'candidates':>12}{'cand %':>9}")
        for jk, g in jr.groupby("journal_key"):
            na = g["article_pmcid"].nunique()
            nc = int(g["candidate"].sum())
            rep.w(f"  {jk:<10}{na:>10,}{len(g):>11,}{nc:>12,}{100*nc/len(g):>8.3f}%")
        rep.w("")

    rep.w("  Strong fabrication candidates (Topaz Category 1/2/3):")
    rep.w(f"    Category 1 (DOI unresolved)          : {int(jr['cat1'].sum()):>6,}")
    rep.w(f"    Category 2 (DOI resolves, title<50)  : {int(jr['cat2'].sum()):>6,}")
    rep.w(f"    Category 3 (no identifier, no match) : {int(jr['cat3'].sum()):>6,}")
    rep.w(f"    Any strong candidate                 : {int(jr['candidate'].sum()):>6,}  "
          f"({100*jr['candidate'].mean():.3f}% of references)")
    rep.w("")

    rep.w("=" * 76)
    rep.w("  PRE / POST-CHATGPT COMPARISON  (boundary: proxy year 2022 | 2023)")
    rep.w("=" * 76)
    _window_block(rep, jr, art, 0, "Window A -- full corpus")
    _window_block(rep, jr, art, 2017, "Window B -- articles after 2016 (proxy year >= 2017)")
    _window_block(rep, jr, art, 2019, "Window C -- balanced 4+4 (proxy year >= 2019)")

    rep.w("=" * 76)
    rep.w("  ADDITIONAL SENSITIVITY")
    rep.w("=" * 76)
    _sensitivity_2018(rep, jr, art)

    rep.w("=" * 76)
    rep.w("  ARTICLE-LEVEL CANDIDATE RATE BY PROXY YEAR")
    rep.w("=" * 76)
    for y in range(2000, 2027):
        sub = art[art["proxy_year"] == y]
        if len(sub) < 10:
            continue
        k, n = int(sub["has_candidate"].sum()), len(sub)
        lo, hi = wilson(k, n)
        era = "post" if y >= 2023 else "pre "
        bar = "#" * int(round(k / n * 45))
        rep.w(f"  {y}  n={n:>6}  {100*k/n:>5.1f}%  [{100*lo:5.1f}-{100*hi:5.1f}]  {era}  {bar}")
    rep.w("")

    return rep.text()


# =============================================================================
# SECTION 7 - ORCHESTRATOR: FETCH / VERIFY / ANALYZE STAGES
# =============================================================================

VERIFIED_FIELDS = RAW_FIELDS + VERDICT_FIELDS


def reference_from_row(row: dict) -> Reference:
    """Rebuild a Reference from a raw-CSV row, with safe type coercion."""
    def s(v):
        if v is None:
            return None
        v = str(v).strip()
        return v or None

    def i(v):
        v = s(v)
        if v is None:
            return None
        try:
            return int(float(v))
        except ValueError:
            return None

    return Reference(
        journal_key=str(row["journal_key"]),
        article_pmcid=str(row["article_pmcid"]),
        ref_id=str(row.get("ref_id") or ""),
        publication_type=s(row.get("publication_type")),
        first_author_surname=s(row.get("first_author_surname")),
        title=s(row.get("title")),
        source=s(row.get("source")),
        year=i(row.get("year")),
        volume=s(row.get("volume")),
        issue=s(row.get("issue")),
        fpage=s(row.get("fpage")),
        lpage=s(row.get("lpage")),
        doi=s(row.get("doi")),
        pmid=clean_pmid(row.get("pmid")),
        pmcid=s(row.get("pmcid")),
    )


def _seen_articles(csv_path: Path) -> set[str]:
    if not csv_path.exists():
        return set()
    with csv_path.open(encoding="utf-8") as f:
        return {row["article_pmcid"] for row in csv.DictReader(f)}


def stage_fetch(journal: JournalConfig, n_articles: int | None, seed: int) -> Path:
    raw_path = DATA / f"{journal.key}_refs_raw.csv"
    log.info("[fetch] %s (%s)", journal.name, journal.key)

    pmcids = search_journal_pmcids(journal)
    log.info("[fetch] %s: %d articles in PMC", journal.key, len(pmcids))

    done = _seen_articles(raw_path)
    pmcids = [p for p in pmcids if p not in done]
    if done:
        log.info("[fetch] %s: resuming, %d already fetched", journal.key, len(done))

    if n_articles is not None:
        random.Random(seed).shuffle(pmcids)
        pmcids = pmcids[:n_articles]
        log.info("[fetch] %s: pilot mode, %d articles", journal.key, len(pmcids))

    DATA.mkdir(parents=True, exist_ok=True)
    is_new = not raw_path.exists()
    with raw_path.open("a", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=RAW_FIELDS, extrasaction="ignore")
        if is_new:
            writer.writeheader()
        current = None
        for ref in tqdm(collect_journal_references(journal, pmcids),
                        desc=f"fetch {journal.key}", unit="ref"):
            writer.writerow(asdict(ref))
            if ref.article_pmcid != current:
                current = ref.article_pmcid
                fh.flush()
    log.info("[fetch] %s: wrote %s", journal.key, raw_path)
    return raw_path


def _verified_keys(csv_path: Path) -> set[tuple[str, str]]:
    """Return (article_pmcid, ref_id) pairs already verified with a final
    (non-error) status, so a re-run retries only errors and new rows."""
    if not csv_path.exists():
        return set()
    keys: set[tuple[str, str]] = set()
    with csv_path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row.get("verify_status") and row["verify_status"] != "error":
                keys.add((row["article_pmcid"], row["ref_id"]))
    return keys


def stage_verify(journal: JournalConfig) -> Path:
    raw_path = DATA / f"{journal.key}_refs_raw.csv"
    out_path = DATA / f"{journal.key}_verified.csv"
    if not raw_path.exists():
        log.warning("[verify] %s: no raw file (%s); run fetch first", journal.key, raw_path)
        return out_path

    with raw_path.open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    log.info("[verify] %s: %d references", journal.key, len(rows))

    done = _verified_keys(out_path)
    if done:
        log.info("[verify] %s: resuming, %d already verified", journal.key, len(done))

    is_new = not out_path.exists()
    with out_path.open("a", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=VERIFIED_FIELDS, extrasaction="ignore")
        if is_new:
            writer.writeheader()
        for row in tqdm(rows, desc=f"verify {journal.key}", unit="ref"):
            key = (row["article_pmcid"], row.get("ref_id") or "")
            if key in done:
                continue
            ref = reference_from_row(row)
            v = verify(ref)
            writer.writerow({**asdict(ref), **verdict_to_row(v)})
            fh.flush()
    log.info("[verify] %s: wrote %s", journal.key, out_path)
    return out_path


def stage_analyze(journals: list[JournalConfig]) -> None:
    verified = [DATA / f"{j.key}_verified.csv" for j in journals]
    present = [p for p in verified if p.exists()]
    missing = [p for p in verified if not p.exists()]
    if missing:
        log.warning("[analyze] missing verified files: %s",
                    ", ".join(p.name for p in missing))
    if not present:
        log.error("[analyze] no verified CSVs found; run fetch + verify first")
        return

    # Build the pooled CSV.
    pooled = DATA / "pooled_verified.csv"
    with pooled.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=VERIFIED_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for p in present:
            with p.open(encoding="utf-8") as fh:
                for row in csv.DictReader(fh):
                    writer.writerow(row)
    log.info("[analyze] pooled CSV: %s", pooled)

    text = run_analysis(present)
    report = ROOT / "analysis_report.md"
    report.write_text("```\n" + text + "\n```\n", encoding="utf-8")
    log.info("[analyze] report: %s", report)


# =============================================================================
# SECTION 8 - COMMAND LINE INTERFACE
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Tier 1 radiation-oncology citation-integrity audit "
                    "(single-file consolidated pipeline)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python tier1_citation_audit.py                       # all 7, all stages\n"
            "  python tier1_citation_audit.py --journal adro        # one journal\n"
            "  python tier1_citation_audit.py --stage fetch         # fetch only\n"
            "  python tier1_citation_audit.py --stage verify        # re-verify\n"
            "  python tier1_citation_audit.py --stage analyze       # analysis only\n"
            "  python tier1_citation_audit.py --n 25                # pilot mode"
        ),
    )
    p.add_argument("--journal", default="all",
                   help="journal key (adro|ctro|radonc|phro|tipsro|jrr|rpor) "
                        "or 'all' (default: all)")
    p.add_argument("--stage", default="all",
                   choices=["fetch", "verify", "analyze", "all"],
                   help="pipeline stage to run (default: all)")
    p.add_argument("--n", type=int, default=None, metavar="N",
                   help="pilot mode: only N random articles per journal")
    p.add_argument("--seed", type=int, default=42, help="pilot sampling seed")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    journals = load_journals()

    ok, warnings = environment_ok()
    for w in warnings:
        log.warning(w)
    if not ok and args.stage in ("fetch", "verify", "all"):
        log.warning("Continuing without NCBI_EMAIL set -- expect throttling.")

    if args.journal == "all":
        selected = list(journals.values())
    elif args.journal in journals:
        selected = [journals[args.journal]]
    else:
        log.error("unknown journal '%s'; valid: %s, all",
                  args.journal, ", ".join(journals))
        return 1

    if args.stage in ("fetch", "all"):
        for j in selected:
            stage_fetch(j, args.n, args.seed)
    if args.stage in ("verify", "all"):
        for j in selected:
            stage_verify(j)
    if args.stage in ("analyze", "all"):
        stage_analyze(selected)

    log.info("done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
