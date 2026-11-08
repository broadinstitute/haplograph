#!/usr/bin/env python3
"""Compare cohort HLA frequencies with Allele Frequency Net (AFND).

Reads the TSV written by ``hla_allele_frequency.py`` and computes the Pearson
R² between each ancestry group and the matching public frequencies for the six
superpopulations AFR, AMR, EAS, EUR, SAS, and MID.

The cohort frequency used here is the within-population haplotype frequency,

    local_af = count(allele, group) / haplotypes of that gene in that group

which is the quantity AFND stores as ``geneFrequency``. The input column
``freq`` is each allele's share of every haplotype of that gene pooled across
groups; the correlation uses the within-population frequency above instead.

Public frequencies are the sample-size-weighted mean of AFND gold-standard
populations (each gold population contributes proportionally to how many
individuals it typed). ``--weight equal`` gives every population the same
weight instead. Alleles are compared at a shared field resolution (default 2,
``A*01:01``). A population typed at coarser resolution is left out of that
locus. Higher-resolution alleles are rolled up: a study typed at 4 fields
contributes the sum of its children to each 2-field allele. Alleles below
0.1% in both the cohort and AFND are left out of the correlation
(``--min-af``, default 0.001) so registry-only singletons do not dominate R².

A gene/population cell is dropped when it has fewer than ``--min-ref-n`` typed
individuals (default 500), or fewer than ``--min-ref-pops`` populations
(default 2) without reaching ``LARGE_REF_N``. AFND's gold data for the
non-classical genes is one small study per superpopulation — HLA-E in EUR is a
single series of 82 — and because those alleles sit at 40-85% frequency they
are the highest-leverage points in the scatter. Left in, they set the
regression line for every panel.

Every population also has to report frequencies that sum into
``[MIN_FREQ_MASS, MAX_FREQ_MASS]``. The upper bound matters: AFND's
``geneFrequency`` column is not always an allele frequency, and one study
stored there on a percentage scale can outweigh a registry a thousand times
its size.

``r2`` measures correlation, not agreement: a cohort uniformly twice AFND
still scores 1.0. Read it next to ``slope``, which a matching spectrum puts
near 1, and prefer the per-gene rows over the pooled ``(all)`` row — pooling
genes lets the spread *between* genes stand in for the within-gene agreement
that actually tests the calls.

Two presets select which AFND populations back each superpopulation. Both
filter the populations an AFND ethnicity label pools in but an ancestry-
labeled cohort does not: "Asian" mixes Han Chinese and Malay cohorts into the
Indian subcontinent, and "Caucasoid" covers the Volga Turkic and North
Caucasus series that AFND files under Russia. See the filters on ``_REGISTRY``.

registry (default)
    Ethnicity labels, worldwide. The large gold registries dominate the
    weighted mean: European Caucasian and DKMS for EUR, African American for
    AFR, NMDP Hispanic for AMR, Oriental registries for EAS, South Asian Indian
    for SAS, Arab registries for MID. This is the usual public catalog for an
    ancestry-labeled biobank cohort.
region
    The same comparison restricted to a geographical region, which is what
    ``geographicalRegion=Europe`` means for EUR. East Asia keeps the Oriental
    label so the large Russian series inside "North-East Asia" stay out of EAS.
    AMR keeps Hispanic and Mestizo populations in the Americas, including the
    NMDP Hispanic registries.

The populations table records every gold population that entered a mean and
the fraction of the weight it carried, so a registry-dominated estimate is
visible in the output rather than hidden inside the average.

Notebook
--------
    import sys; sys.path.append("scripts")
    from hla_afnd_r2 import AFNDComparison

    fit = AFNDComparison.from_table("hla_allele_frequency.tsv")
    # fit = AFNDComparison.from_table(af.frequency_matrix())
    fit.summary          # R² per gene and superpopulation
    fit.comparison.head()
    fit.plot()           # matplotlib Figure, renders inline

The same call is written out in ``notebooks/hla_afnd_r2.ipynb``.

CLI
---
    python scripts/hla_afnd_r2.py --freq hla_allele_frequency.tsv --out hla_afnd_r2.tsv

Writes the per-allele comparison, a per-gene R² summary, the population
weights, a scatter plot, and a cache of the AFND responses next to ``--out``.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = [
    "AFNDComparison",
    "AFNDError",
    "PRESETS",
    "SUPERPOPS",
    "aggregate_afnd",
    "compare",
    "correlation_tables",
    "fetch_afnd",
    "load_local_frequencies",
    "plot_comparison",
]

SEARCH_URL = "https://www.allelefrequencies.net/api/hla/allele-frequencies/search"
USER_AGENT = "haplograph-hla-afnd-r2/1.0 (research; allele frequency comparison)"

# Order matches hla_allele_frequency.SUPERPOP_ORDER, without OTH / Unknown.
SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS", "MID"]

CLASSICAL_GENES = [
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "HLA-DRB1",
    "HLA-DQA1",
    "HLA-DQB1",
    "HLA-DPA1",
    "HLA-DPB1",
]

# A gold population's allele frequencies have to sum into this band at one
# resolution before it is allowed into the mean. The lower bound drops studies
# that mix 1-field and 2-field calls or report only part of the locus. The
# upper bound drops studies whose geneFrequency column is not an allele
# frequency at all: "Czech Republic pop 3" stores DQA1 and DQB1 phenotype
# frequencies as percentages there, summing to 185, and weighted in at face
# value it alone pushes EUR DQA1*05:05 to a frequency of 1.98.
MIN_FREQ_MASS = 0.8
MAX_FREQ_MASS = 1.1

# A cell backed by a single population is normally a lone small study, but one
# large registry is a better reference than several small series. Above this
# many typed individuals the population-count floor does not apply: the MID
# DPB1 cell is one NMDP registry of 28,927 donors and is the best public data
# for that locus and group.
LARGE_REF_N = 5_000

# AFND's ethnicity labels are anthropological, not genetic, so some of them
# pool populations an ancestry-labeled biobank cohort would never group
# together. A preset query may therefore carry three post-filters applied to
# the populations the search returns:
#
#   countries  ISO codes to keep.
#   include    population names re-admitted whatever their country, for the
#              diaspora registries AFND files under their country of residence
#              ("USA NMDP South Asian Indian" is in the USA).
#   exclude    population names dropped whatever their country.
#
# The populations table records everything that survived, so a filter's effect
# is visible in the output rather than hidden in a constant here.

_SOUTH_ASIA = {"IND", "PAK", "BGD", "LKA", "NPL", "BTN", "MDV", "AFG"}
_SOUTH_ASIAN_NAMES = r"south asian|indian|sri lank|pakistan|banglades|nepal"
# "West Indian" and "American Indian" are not South Asian.
_NOT_SOUTH_ASIAN_NAMES = r"west indian|american indian|amerind"

# Russia is deliberately absent: AFND files every Russian series under the
# North-East Asia region, and the Caucasoid ones are Volga Turkic or North
# Caucasus cohorts (Tatar, Bashkir, Kabardin, ...), not European ancestry.
_EUROPE = {
    "ALB", "AUT", "BEL", "BGR", "BIH", "BLR", "CHE", "CYP", "CZE", "DEU",
    "DNK", "ENG", "ESP", "EST", "FIN", "FRA", "GBR", "GRC", "HRV", "HUN",
    "IRL", "ISL", "ITA", "LTU", "LUX", "LVA", "MDA", "MKD", "MLT", "NLD",
    "NOR", "POL", "PRT", "ROU", "SCG", "SRB", "SVK", "SVN", "SWE", "UKR",
    "WLS",
}
# European-ancestry registries and cord-blood banks outside Europe.
_EUROPEAN_DIASPORA = {"USA", "CAN", "AUS", "NZL"}
# ... and the labels that mark a non-European series inside those countries.
_NON_EUROPEAN_NAMES = (
    r"hispanic|mestizo|jewish|arab|african|black|asian|amerind|native|"
    r"chinese|japanese|korean|filipino|vietnam|hawaii|pacific|middle east"
)

# Ethnicity worldwide. Sample-size weighting then follows the large registries.
# AMR is admixed American (Hispanic / Mestizo), not Amerindian.
# MID is Arab. Jewish series are a different, internally diverse collection
# and are not pooled into MID.
_REGISTRY = {
    "AFR": [{"ethnicity": "Black"}],
    "AMR": [{"ethnicity": "Hispanic"}, {"ethnicity": "Mestizo"}],
    # Siberian cohorts are filed under Russia and labelled Oriental; they are
    # not what an EAS cohort label means.
    "EAS": [{"ethnicity": "Oriental", "exclude": r"^russia"}],
    "EUR": [
        {
            "ethnicity": "Caucasoid",
            "countries": _EUROPE | _EUROPEAN_DIASPORA,
            "exclude": _NON_EUROPEAN_NAMES,
        }
    ],
    # "Asian" is not South Asian: it pools Han Chinese and Malay cohorts with
    # the Indian subcontinent. Without this filter the SAS reference for the DP
    # loci comes out ~90% China Zhejiang Han.
    "SAS": [
        {
            "ethnicity": "Asian",
            "countries": _SOUTH_ASIA,
            "include": _SOUTH_ASIAN_NAMES,
            "exclude": _NOT_SOUTH_ASIAN_NAMES,
        }
    ],
    "MID": [{"ethnicity": "Arab"}],
}

# Geographical regions, with an ethnicity guard where the region label mixes
# ancestries (Russian cohorts are tagged North-East Asia; Jewish cohorts
# dominate Western Asia by sample size).
_REGION = {
    "AFR": [{"geographicalRegion": "Sub-Saharan Africa", "ethnicity": "Black"}],
    "AMR": [
        {"geographicalRegion": "North America", "ethnicity": "Hispanic"},
        {"geographicalRegion": "North America", "ethnicity": "Mestizo"},
        {"geographicalRegion": "South and Central America", "ethnicity": "Hispanic"},
        {"geographicalRegion": "South and Central America", "ethnicity": "Mestizo"},
    ],
    "EAS": [
        {
            "geographicalRegion": "North-East Asia",
            "ethnicity": "Oriental",
            "exclude": r"^russia",
        },
        {"geographicalRegion": "South-East Asia", "ethnicity": "Oriental"},
    ],
    "EUR": [{"geographicalRegion": "Europe", "ethnicity": "Caucasoid"}],
    "SAS": [{"geographicalRegion": "South Asia"}],
    "MID": [
        {"geographicalRegion": "Western Asia", "ethnicity": "Arab"},
        {"geographicalRegion": "North Africa", "ethnicity": "Arab"},
    ],
}

PRESETS = {"registry": _REGISTRY, "region": _REGION}

STANDARDS = {"gold": "goldOnly", "gold-silver": "goldAndSilver", "all": ""}

GENE_COLORS = {
    "HLA-A": "#2a78d6",
    "HLA-B": "#eb6834",
    "HLA-C": "#1baf7a",
    "HLA-DRB1": "#eda100",
    "HLA-DQA1": "#e87ba4",
    "HLA-DQB1": "#6250d6",
    "HLA-DPA1": "#008300",
    "HLA-DPB1": "#e34948",
}
# Standard superpopulation colors. MID is the teal listed first for that group.
GROUP_COLORS = {
    "AFR": "#FF7F00",
    "AMR": "#ECB01F",
    "EAS": "#778500",
    "EUR": "#377EB8",
    "SAS": "#C44CFD",
    "MID": "#4A90E2",
    "OCE": "#E41A1C",
}
OTHER_COLOR = "#898781"
SURFACE = "#fcfcfb"
INK = "k"
INK_MUTED = "grey"
GRID = "#e1e0d9"

_GROUP_RE = re.compile(r"[0-9][GP]$", re.IGNORECASE)
_EXPR_RE = re.compile(r"(?<=\d)[NLSQAC]$", re.IGNORECASE)
_NEW_RE = re.compile(r"new", re.IGNORECASE)

_CACHE_COLUMNS = [
    "locus",
    "geographical_region",
    "query_ethnicity",
    "population_standard",
    "population_id",
    "population",
    "country",
    "ethnicity",
    "sample_size",
    "allele",
    "gene_frequency",
    "pheno_frequency",
]
_QUERY_KEYS = ["locus", "geographical_region", "query_ethnicity", "population_standard"]


class AFNDError(ValueError):
    """Raised for a bad frequency table or an AFND response we cannot use."""


# --- alleles --------------------------------------------------------------


def canonical_gene(gene: str) -> str:
    """``hla-drb1`` / ``HLA-DRB1`` -> ``HLA-DRB1``."""
    gene = str(gene).strip()
    if gene.upper().startswith("HLA-"):
        return "HLA-" + gene.split("-", 1)[1].upper()
    return gene


def locus_of(gene: str) -> str:
    """``HLA-DRB1`` -> ``DRB1``, the AFND locus code."""
    gene = canonical_gene(gene)
    return gene.split("-", 1)[1] if gene.startswith("HLA-") else gene


def allele_fields(allele: str) -> list[str] | None:
    """Colon-separated fields of an allele, or None for G/P groups and junk.

    A trailing expression letter (N, L, S, A, Q, C) is removed from the last
    field. Ambiguous group codes (``A*02:01G``) are rejected so they are not
    treated as a real allele.
    """
    if not isinstance(allele, str) or "*" not in allele:
        return None
    name = allele.split("*", 1)[1].strip()
    name = _NEW_RE.sub("", name)
    if not name or _GROUP_RE.search(name.replace(":", "")):
        return None
    fields = [part for part in name.split(":") if part]
    if not fields:
        return None
    fields[-1] = _EXPR_RE.sub("", fields[-1])
    fields = [part for part in fields if part]
    return fields or None


def canonical_allele(gene: str, allele: str, resolution: int) -> str | None:
    """Trim ``allele`` to ``resolution`` fields under ``gene``.

    Returns None when the allele has fewer fields than requested: a 1-field
    call cannot be expanded to ``A*01:01``.
    """
    fields = allele_fields(allele)
    if not fields or len(fields) < resolution:
        return None
    return f"{canonical_gene(gene)}*{':'.join(fields[:resolution])}"


def _field_count(allele: str) -> int:
    fields = allele_fields(allele)
    return len(fields) if fields else 0


# --- cohort frequencies ---------------------------------------------------


def _is_wide_matrix(frame: pd.DataFrame) -> bool:
    """True for ``frequency_matrix()`` / ``top_alleles()``, false for the tidy table."""
    columns = {str(col).strip() for col in frame.columns}
    if {"Gene", "Allele", "Group", "count"} <= columns:
        return False
    if isinstance(frame.index, pd.MultiIndex):
        return True
    if frame.index.name in ("Allele", "Gene"):
        return True
    return "Allele" in columns and "Group" not in columns and "count" not in columns


def _from_wide_matrix(frame: pd.DataFrame) -> pd.DataFrame:
    """Melt an allele × population matrix into Gene, Allele, Group, local_af.

    ``frequency_matrix()`` stores each allele's share of every haplotype of
    that gene (its ``freq`` column). Dividing a group's column by the sum of
    that column, within the gene, recovers the within-population frequency.
    ``mode='composition'`` is the population mix of each allele and is rejected,
    because that matrix cannot be turned back into an allele frequency.
    """
    wide = frame.copy()
    if isinstance(wide.index, pd.MultiIndex) or wide.index.name is not None:
        wide = wide.reset_index()
    wide.columns = [str(col).strip() for col in wide.columns]
    if "Allele" not in wide.columns:
        for candidate in ("index", "level_0", "level_1"):
            if candidate in wide.columns and wide[candidate].astype(str).str.contains("*", regex=False).any():
                wide = wide.rename(columns={candidate: "Allele"})
                break
    if "Allele" not in wide.columns:
        raise AFNDError(
            "frequency matrix has no Allele index. Pass af.frequency_matrix() or af.table."
        )
    if "Gene" not in wide.columns:
        # top_alleles() is one gene; the gene is the prefix of HLA-A*01:01.
        for candidate in ("level_0",):
            if candidate in wide.columns:
                wide = wide.rename(columns={candidate: "Gene"})
        if "Gene" not in wide.columns:
            wide["Gene"] = wide["Allele"].astype(str).map(lambda allele: allele.split("*", 1)[0])

    value_vars = []
    for col in list(wide.columns):
        if col in ("Gene", "Allele"):
            continue
        numeric = pd.to_numeric(wide[col], errors="coerce")
        if numeric.notna().any():
            wide[col] = numeric.fillna(0.0)
            value_vars.append(col)
    if not value_vars:
        raise AFNDError("frequency matrix has no population columns")

    long = wide.melt(id_vars=["Gene", "Allele"], value_vars=value_vars, var_name="Group", value_name="local_af")
    long["local_af"] = pd.to_numeric(long["local_af"], errors="coerce").fillna(0.0)
    mass = long.groupby("Gene")["local_af"].sum()
    n_alleles = long.groupby("Gene")["Allele"].nunique()
    ratio = (mass / n_alleles.reindex(mass.index)).median()
    if len(mass) and float((mass > 1.5).mean()) > 0.5 and ratio > 0.8:
        raise AFNDError(
            "this matrix is the within-allele population mix from "
            "frequency_matrix(mode='composition'). Pass af.frequency_matrix() "
            "with the default mode, or af.table."
        )
    return long[["Gene", "Allele", "Group", "local_af"]]


def load_local_frequencies(
    freq: str | Path | pd.DataFrame,
    resolution: int = 2,
    groups: list[str] | None = None,
) -> pd.DataFrame:
    """Within-population allele frequencies from a cohort frequency table.

    Accepts the tidy table from ``hla_allele_frequency.py`` (columns Gene,
    Allele, Group, count), ``AlleleFrequencies.table``, or the wide matrix
    from ``af.frequency_matrix()`` / ``af.top_alleles()``. A group's column in
    that matrix is each allele's share of the whole gene; it is divided by
    the column sum so the result is the frequency inside the group.

    Groups other than the six superpopulations are dropped. Alleles coarser
    than ``resolution`` are dropped. Alleles that collapse to the same name
    have their frequencies added. The matrix has no haplotype counts, so
    ``local_count`` and ``local_n`` are left empty.
    """
    if resolution < 1:
        raise AFNDError(f"resolution must be >= 1, got {resolution}")
    frame = _read_table(freq) if not isinstance(freq, pd.DataFrame) else freq.copy()
    if _is_wide_matrix(frame):
        long = _from_wide_matrix(frame)
        has_counts = False
    else:
        frame.columns = [str(col).strip() for col in frame.columns]
        missing = {"Gene", "Allele", "Group", "count"} - set(frame.columns)
        if missing:
            raise AFNDError(
                "frequency table is missing column(s): "
                + ", ".join(sorted(missing))
                + ". Pass af.table, af.frequency_matrix(), or the TSV from hla_allele_frequency.py."
            )
        long = frame
        has_counts = True

    long["Gene"] = long["Gene"].map(canonical_gene)
    long["Group"] = long["Group"].astype(str).str.strip().str.upper()
    if groups is not None:
        groups = [str(group).strip().upper() for group in groups]
    keep = set(groups) if groups else set(SUPERPOPS)
    unknown = sorted(set(long["Group"]) - set(SUPERPOPS) - {"", "NAN"})
    long = long[long["Group"].isin(keep)].copy()
    if long.empty:
        raise AFNDError(
            "no rows left for the six superpopulations "
            + ", ".join(SUPERPOPS)
            + (f" (table groups included {', '.join(unknown)})" if unknown else "")
        )

    trimmed = [
        canonical_allele(gene, allele, resolution)
        for gene, allele in zip(long["Gene"], long["Allele"])
    ]
    dropped = int(sum(allele is None for allele in trimmed))
    if dropped:
        print(
            f"note: dropped {dropped} rows coarser than {resolution} fields or tagged as allele groups",
            file=sys.stderr,
        )
    long["Allele"] = trimmed
    long = long.dropna(subset=["Allele"])
    if long.empty:
        raise AFNDError(f"no alleles left at {resolution}-field resolution")

    if has_counts:
        long["count"] = pd.to_numeric(long["count"], errors="coerce")
        long = long.dropna(subset=["count"])
        long = long[long["count"] > 0]
        long = long.groupby(["Gene", "Allele", "Group"], as_index=False)["count"].sum()
        totals = long.groupby(["Gene", "Group"])["count"].transform("sum")
        long["local_af"] = long["count"] / totals
        long["local_n"] = totals.astype(int)
        long = long.rename(columns={"count": "local_count"})
    else:
        long["local_af"] = pd.to_numeric(long["local_af"], errors="coerce").fillna(0.0)
        long = long[long["local_af"] > 0]
        long = long.groupby(["Gene", "Allele", "Group"], as_index=False)["local_af"].sum()
        totals = long.groupby(["Gene", "Group"])["local_af"].transform("sum")
        long = long[totals > 0].copy()
        long["local_af"] = long["local_af"] / long.groupby(["Gene", "Group"])["local_af"].transform("sum")
        long["local_count"] = pd.NA
        long["local_n"] = pd.NA
    return long[["Gene", "Allele", "Group", "local_af", "local_count", "local_n"]]


def _read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path).expanduser()
    if not path.exists():
        raise AFNDError(f"frequency table not found: {path}")
    return pd.read_csv(path, sep=None, engine="python", encoding_errors="replace")


# --- AFND download --------------------------------------------------------


def _get_json(params: dict, timeout: float = 180, attempts: int = 4) -> dict:
    query = {key: value for key, value in params.items() if value not in (None, "")}
    url = SEARCH_URL + "?" + urllib.parse.urlencode(query)
    delay = 1.0
    last = "no response"
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(
                url,
                headers={"Accept": "application/json", "User-Agent": USER_AGENT},
            )
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.load(response)
        except urllib.error.HTTPError as exc:
            body = exc.read().decode("utf-8", errors="replace")[:300]
            last = f"HTTP {exc.code}: {body}"
            if exc.code not in (429, 500, 502, 503, 504) or attempt == attempts - 1:
                raise AFNDError(f"AFND request failed ({last})") from exc
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last = str(exc)
            if attempt == attempts - 1:
                raise AFNDError(f"AFND request failed ({last})") from exc
        time.sleep(delay)
        delay *= 2
    raise AFNDError(f"AFND request failed ({last})")


def _fetch_query(locus: str, region: str, ethnicity: str, standard: str) -> list[dict]:
    """Every positive-frequency row for one locus / region / ethnicity filter.

    The search endpoint returns ``totalElements`` but ignores page offsets, so
    the call asks for the whole result in one shot. A shortfall of a row or two
    that a larger page does not fill is a record the server counts and does not
    serialize (Oriental DPB1 is 698 of 699); those rows are kept. A large
    shortfall is a truncated page and still raises.
    """
    params = {
        "locus": locus,
        "geographicalRegion": region,
        "ethnicity": ethnicity,
        "populationStandard": standard,
        "showFrequencies": "positives",
        "size": 1,
    }
    probe = _get_json(params)
    total = int(probe.get("totalElements") or 0)
    if total == 0:
        return []
    if total > 200_000:
        raise AFNDError(f"AFND reports {total} rows for {locus}; refusing that download")
    if total == 1:
        return list(probe.get("content") or [])
    params["size"] = total
    data = _get_json(params)
    content = list(data.get("content") or [])
    reported = int(data.get("totalElements") or total)
    if len(content) < reported:
        # A second, larger page distinguishes a missing record from a cap.
        params["size"] = reported + max(50, reported // 10)
        retry = _get_json(params)
        retry_rows = list(retry.get("content") or [])
        if len(retry_rows) > len(content):
            content = retry_rows
            reported = int(retry.get("totalElements") or reported)
    missing = reported - len(content)
    if missing <= 0:
        return content
    label = " ".join(part for part in (locus, region, ethnicity) if part)
    if missing <= max(2, reported // 100):
        print(
            f"note: AFND returned {len(content)} of {reported} rows for {label}; "
            "the missing row is not in the response",
            file=sys.stderr,
        )
        return content
    raise AFNDError(
        f"AFND returned {len(content)} of {reported} rows for locus {locus}. "
        "The result was truncated and there is no working page offset."
    )


def _query_frame(content: list[dict], locus: str, region: str, ethnicity: str, standard: str) -> pd.DataFrame:
    standard_key = standard or "all"
    if not content:
        row = {column: "" for column in _CACHE_COLUMNS}
        row.update(
            {
                "locus": locus,
                "geographical_region": region,
                "query_ethnicity": ethnicity,
                "population_standard": standard_key,
            }
        )
        return pd.DataFrame([row])
    records = []
    for item in content:
        frequency = item.get("geneFrequency")
        if frequency is None:
            frequency = item.get("calculatedFrequency")
        # AFND returns geneFrequency 0.0 -- not null -- for studies that only
        # published a phenotype frequency, and showFrequencies=positives keeps
        # those rows because the phenotype frequency is positive. Left as 0.0
        # they would enter the weighted mean as a real observation of absence,
        # so they are recorded as missing instead.
        pheno = item.get("phenoFrequency")
        try:
            pheno_value = float(pheno)
        except (TypeError, ValueError):
            pheno_value = 0.0
        if pheno_value > 0 and (frequency is None or float(frequency) == 0.0):
            frequency = None
        records.append(
            {
                "locus": locus,
                "geographical_region": region,
                "query_ethnicity": ethnicity,
                "population_standard": standard_key,
                "population_id": "" if item.get("populationId") is None else str(item.get("populationId")),
                "population": item.get("population") or "",
                "country": item.get("countryCode") or "",
                "ethnicity": item.get("ethnicity") or "",
                "sample_size": item.get("sampleSize"),
                "allele": item.get("allele") or "",
                "gene_frequency": frequency,
                "pheno_frequency": pheno,
            }
        )
    return pd.DataFrame(records, columns=_CACHE_COLUMNS)


def _load_cache(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=_CACHE_COLUMNS)
    frame = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        encoding="utf-8",
        encoding_errors="replace",
    )
    if "pheno_frequency" not in frame.columns:
        # Written before phenotype-only rows were recorded as missing rather
        # than as a gene frequency of 0.0. Those zeros are still in this file
        # and would still bias the mean down, so the cache has to be refetched.
        print(
            f"warning: {path} predates the phenotype-frequency fix and still holds "
            "geneFrequency 0.0 for studies that only published phenotype "
            "frequencies. Delete it or pass --refresh / refresh=True.",
            file=sys.stderr,
        )
    for column in _CACHE_COLUMNS:
        if column not in frame.columns:
            frame[column] = ""
    return frame[_CACHE_COLUMNS]


def _write_cache(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    out = frame[_CACHE_COLUMNS].copy()
    out["gene_frequency"] = pd.to_numeric(out["gene_frequency"], errors="coerce")
    out["pheno_frequency"] = pd.to_numeric(out["pheno_frequency"], errors="coerce")
    out["sample_size"] = pd.to_numeric(out["sample_size"], errors="coerce")
    out.to_csv(path, sep="\t", index=False)


def _planned_queries(
    preset: str,
    superpops: list[str],
    loci: list[str],
    standard: str,
) -> list[tuple[str, str, str, str, str, dict]]:
    if preset not in PRESETS:
        raise AFNDError(f"preset must be one of {', '.join(PRESETS)}, got {preset}")
    standard_param = STANDARDS.get(standard)
    if standard_param is None:
        raise AFNDError(f"standard must be one of {', '.join(STANDARDS)}, got {standard}")
    jobs = []
    for superpop in superpops:
        for query in PRESETS[preset][superpop]:
            region = query.get("geographicalRegion", "")
            ethnicity = query.get("ethnicity", "")
            filters = {
                key: query[key]
                for key in ("countries", "include", "exclude")
                if key in query
            }
            for locus in loci:
                jobs.append((superpop, locus, region, ethnicity, standard_param, filters))
    return jobs


def fetch_afnd(
    genes: list[str],
    superpops: list[str],
    *,
    preset: str = "registry",
    standard: str = "gold",
    cache: str | Path | None = None,
    refresh: bool = False,
    verbose: bool = True,
) -> pd.DataFrame:
    """Download (or reload) the AFND rows needed for these genes and groups.

    The cache is the raw search result, keyed by locus, region, ethnicity, and
    population standard, so a second run does not repeat the download. The
    returned frame is restricted to ``superpops`` and carries a ``superpop``
    column. Marker rows for an empty search have a blank allele.
    """
    loci = []
    for gene in genes:
        locus = locus_of(gene)
        if locus not in loci:
            loci.append(locus)
    jobs = _planned_queries(preset, superpops, loci, standard)
    cache_path = Path(cache).expanduser() if cache else None
    cached = _load_cache(cache_path) if cache_path else pd.DataFrame(columns=_CACHE_COLUMNS)
    if refresh and not cached.empty:
        drop = pd.DataFrame([job[:5] for job in jobs], columns=["superpop", *_QUERY_KEYS])
        # standard stored as "all" when the parameter is omitted
        drop["population_standard"] = drop["population_standard"].replace({"": "all"})
        cached = cached.merge(drop[_QUERY_KEYS], on=_QUERY_KEYS, how="left", indicator=True)
        cached = cached[cached["_merge"] == "left_only"].drop(columns="_merge")

    wrote_cache = False
    for superpop, locus, region, ethnicity, standard_param, _filters in jobs:
        standard_key = standard_param or "all"
        have = cached
        if not have.empty:
            have = have[
                (have["locus"] == locus)
                & (have["geographical_region"] == region)
                & (have["query_ethnicity"] == ethnicity)
                & (have["population_standard"] == standard_key)
            ]
        if not have.empty:
            continue
        label = " ".join(
            part
            for part in (region or None, ethnicity or None, standard_key)
            if part
        )
        if verbose:
            print(f"AFND {superpop} {locus} ({label})", file=sys.stderr)
        content = _fetch_query(locus, region, ethnicity, standard_param)
        if verbose:
            print(f"  {len(content)} rows", file=sys.stderr)
        cached = pd.concat(
            [cached, _query_frame(content, locus, region, ethnicity, standard_param)],
            ignore_index=True,
        )
        cached = cached.drop_duplicates(_QUERY_KEYS + ["population_id", "allele"])
        if cache_path:
            _write_cache(cache_path, cached)
            wrote_cache = True

    if wrote_cache and verbose:
        print(f"cached {cache_path}", file=sys.stderr)

    if cached.empty:
        return pd.DataFrame(columns=_CACHE_COLUMNS + ["superpop"])
    return _attach_superpop(cached, jobs, verbose=verbose)


def _apply_population_filter(part: pd.DataFrame, filters: dict, label: str, verbose: bool) -> pd.DataFrame:
    """Drop the populations an AFND ethnicity label pools in but the cohort does not.

    ``countries`` keeps those ISO codes, ``include`` re-admits a population by
    name whatever its country (diaspora registries), and ``exclude`` drops one
    by name. Marker rows for an empty search carry no population and survive.
    """
    if not filters or part.empty:
        return part
    name = part["population"].astype(str).str.lower()
    keep = pd.Series(True, index=part.index)
    countries = filters.get("countries")
    if countries:
        keep &= part["country"].astype(str).str.upper().isin(set(countries))
        include = filters.get("include")
        if include:
            keep |= name.str.contains(include, regex=True, na=False)
    exclude = filters.get("exclude")
    if exclude:
        keep &= ~name.str.contains(exclude, regex=True, na=False)
    keep |= part["allele"].astype(str).str.len().eq(0)
    if verbose:
        dropped = sorted(set(part.loc[~keep, "population"]))
        if dropped:
            shown = ", ".join(dropped[:4])
            more = ", ..." if len(dropped) > 4 else ""
            print(
                f"note: {label}: left out {len(dropped)} population(s) the AFND "
                f"label pools in ({shown}{more})",
                file=sys.stderr,
            )
    return part[keep]


def _attach_superpop(cached: pd.DataFrame, jobs: list[tuple], verbose: bool = False) -> pd.DataFrame:
    pieces = []
    for superpop, locus, region, ethnicity, standard_param, filters in jobs:
        standard_key = standard_param or "all"
        part = cached[
            (cached["locus"] == locus)
            & (cached["geographical_region"] == region)
            & (cached["query_ethnicity"] == ethnicity)
            & (cached["population_standard"] == standard_key)
        ].copy()
        if part.empty:
            continue
        part["superpop"] = superpop
        part = _apply_population_filter(part, filters, f"{superpop} {locus}", verbose)
        if part.empty:
            continue
        pieces.append(part)
    if not pieces:
        return pd.DataFrame(columns=_CACHE_COLUMNS + ["superpop"])
    frame = pd.concat(pieces, ignore_index=True)
    frame["gene_frequency"] = pd.to_numeric(frame["gene_frequency"], errors="coerce")
    frame["sample_size"] = pd.to_numeric(frame["sample_size"], errors="coerce")
    if "pheno_frequency" in frame.columns:
        frame["pheno_frequency"] = pd.to_numeric(frame["pheno_frequency"], errors="coerce")
    return frame


# --- weighted reference frequencies ---------------------------------------


def _collapse_population(
    alleles: list[str], frequencies: list[float], resolution: int
) -> tuple[dict[str, float] | None, str]:
    """One population's allele frequencies rolled up to ``resolution`` fields.

    Uses only the modal field-count of that population, so a 2-field study and
    a 4-field study are not added together allele by allele.

    Returns ``(frequencies, reason)``. ``frequencies`` is None when the
    population cannot be used, and ``reason`` says why: ``"resolution"`` when
    it is typed coarser than requested, ``"mass_low"`` when the kept
    frequencies sum to less than ``MIN_FREQ_MASS``, and ``"mass_high"`` when
    they sum to more than ``MAX_FREQ_MASS``, which means the column is not an
    allele frequency.
    """
    mass: dict[int, float] = defaultdict(float)
    rows: list[tuple[int, str, float]] = []
    for allele, frequency in zip(alleles, frequencies):
        if frequency is None or not np.isfinite(frequency) or frequency < 0:
            continue
        count = _field_count(allele)
        if count == 0:
            continue
        mass[count] += float(frequency)
        rows.append((count, allele, float(frequency)))
    if not rows:
        return None, "mass_low"
    modal = max(mass, key=mass.get)
    if modal < resolution:
        return None, "resolution"
    if mass[modal] < MIN_FREQ_MASS:
        return None, "mass_low"
    if mass[modal] > MAX_FREQ_MASS:
        return None, "mass_high"
    collapsed: dict[str, float] = defaultdict(float)
    for count, allele, frequency in rows:
        if count != modal:
            continue
        # The gene prefix is replaced by the caller; keep the raw fields here.
        fields = allele_fields(allele)
        if not fields:
            continue
        name = ":".join(fields[:resolution])
        collapsed[name] += frequency
    total = sum(collapsed.values())
    if total < MIN_FREQ_MASS:
        return None, "mass_low"
    if total > MAX_FREQ_MASS:
        return None, "mass_high"
    return dict(collapsed), ""


def aggregate_afnd(
    raw: pd.DataFrame,
    *,
    resolution: int = 2,
    weight: str = "sample",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sample-size-weighted mean allele frequency per superpopulation and gene.

    ``weight="sample"`` weights each gold population by its sample size.
    ``weight="equal"`` weights each population by 1. Alleles absent from a
    population that typed the locus enter the mean as zero, so a private allele
    of a small series does not keep that series' frequency.

    Returns ``(reference, populations)``. ``reference`` has one row per
    superpopulation, gene, and allele with a positive mean. ``populations`` has
    one row per population that passed the resolution filter, with the fraction
    of the weight it received.
    """
    if weight not in ("sample", "equal"):
        raise AFNDError("weight must be 'sample' or 'equal'")
    if raw.empty:
        empty_ref = pd.DataFrame(columns=["Gene", "Allele", "Group", "afnd_af", "afnd_n", "afnd_populations"])
        empty_pop = pd.DataFrame(
            columns=[
                "Group", "Gene", "population_id", "population", "country", "ethnicity",
                "sample_size", "weight", "n_alleles", "n_pheno_only", "freq_sum",
            ]
        )
        return empty_ref, empty_pop

    usable = raw.dropna(subset=["allele"]).copy()
    usable = usable[usable["allele"].astype(str).str.len() > 0]
    usable["gene_frequency"] = pd.to_numeric(usable["gene_frequency"], errors="coerce")
    # Alleles this population published only as a phenotype frequency. They are
    # missing, not zero, so they never reach the mean; the count is reported so
    # a thinly-typed study is visible in the populations table.
    if "pheno_frequency" in usable.columns:
        pheno = pd.to_numeric(usable["pheno_frequency"], errors="coerce")
    else:
        pheno = pd.Series(np.nan, index=usable.index)
    pheno_only = (
        usable[usable["gene_frequency"].isna() & (pheno > 0)]
        .groupby(["superpop", "locus", "population_id"])
        .size()
    )
    usable = usable.dropna(subset=["gene_frequency", "sample_size"])
    usable = usable[usable["sample_size"] > 0]
    usable = usable.drop_duplicates(["superpop", "locus", "population_id", "allele"])

    ref_rows = []
    pop_rows = []
    n_seen = 0
    if usable.empty:
        grouped = []
    else:
        grouped = usable.groupby(["superpop", "locus", "population_id"], sort=False)
        n_seen = grouped.ngroups

    # Collect collapsed populations first, then average within (superpop, locus).
    buckets: dict[tuple[str, str], list[dict]] = defaultdict(list)
    rejected: dict[str, list[str]] = defaultdict(list)
    for (superpop, locus, population_id), part in grouped:
        collapsed, reason = _collapse_population(
            part["allele"].tolist(),
            part["gene_frequency"].tolist(),
            resolution,
        )
        if collapsed is None:
            rejected[reason].append(f"{superpop} {locus} {part['population'].iloc[0]}")
            continue
        gene = f"HLA-{locus}"
        sample_size = float(part["sample_size"].iloc[0])
        buckets[(superpop, gene)].append(
            {
                "n_pheno_only": int(pheno_only.get((superpop, locus, population_id), 0)),
                "population_id": population_id,
                "population": part["population"].iloc[0],
                "country": part["country"].iloc[0],
                "ethnicity": part["ethnicity"].iloc[0],
                "sample_size": sample_size,
                "frequencies": collapsed,
            }
        )

    coarse = len(rejected["resolution"]) + len(rejected["mass_low"])
    if n_seen and coarse:
        print(
            f"note: left out {coarse} of {n_seen} populations typed coarser than "
            f"{resolution} fields or reporting only part of the locus",
            file=sys.stderr,
        )
    # Named rather than counted: a frequency column that is not a frequency is
    # an AFND data problem worth seeing, and it moves the mean a long way.
    if rejected["mass_high"]:
        print(
            "note: left out "
            + str(len(rejected["mass_high"]))
            + " population(s) whose frequencies sum above "
            + f"{MAX_FREQ_MASS} and so are not allele frequencies: "
            + "; ".join(sorted(rejected["mass_high"])),
            file=sys.stderr,
        )

    for (superpop, gene), populations in buckets.items():
        if weight == "equal":
            weights = np.ones(len(populations), dtype=float)
        else:
            weights = np.array([item["sample_size"] for item in populations], dtype=float)
        total_weight = float(weights.sum())
        alleles = set()
        for item in populations:
            alleles.update(item["frequencies"])
        for allele_fields_name in alleles:
            numerator = 0.0
            for item, item_weight in zip(populations, weights):
                numerator += item_weight * item["frequencies"].get(allele_fields_name, 0.0)
            frequency = numerator / total_weight
            if frequency <= 0:
                continue
            ref_rows.append(
                {
                    "Gene": gene,
                    "Allele": f"{gene}*{allele_fields_name}",
                    "Group": superpop,
                    "afnd_af": frequency,
                    "afnd_n": int(round(sum(item["sample_size"] for item in populations))),
                    "afnd_populations": len(populations),
                }
            )
        for item, item_weight in zip(populations, weights):
            pop_rows.append(
                {
                    "Group": superpop,
                    "Gene": gene,
                    "population_id": item["population_id"],
                    "population": item["population"],
                    "country": item["country"],
                    "ethnicity": item["ethnicity"],
                    "sample_size": int(round(item["sample_size"])),
                    "weight": float(item_weight / total_weight),
                    "n_alleles": len(item["frequencies"]),
                    "n_pheno_only": item["n_pheno_only"],
                    "freq_sum": float(sum(item["frequencies"].values())),
                }
            )

    reference = pd.DataFrame(ref_rows)
    populations_out = pd.DataFrame(pop_rows)
    if reference.empty:
        reference = pd.DataFrame(columns=["Gene", "Allele", "Group", "afnd_af", "afnd_n", "afnd_populations"])
    if populations_out.empty:
        populations_out = pd.DataFrame(
            columns=[
                "Group", "Gene", "population_id", "population", "country", "ethnicity",
                "sample_size", "weight", "n_alleles", "n_pheno_only", "freq_sum",
            ]
        )
    else:
        populations_out = populations_out.sort_values(
            ["Group", "Gene", "weight"], ascending=[True, True, False]
        )
    # Invariant: a weighted mean of allele frequencies is itself an allele
    # frequency distribution, so every cell sums to 1. A cell that does not is
    # carrying a row AFND stores on some other scale.
    if not reference.empty:
        cell_sum = reference.groupby(["Group", "Gene"])["afnd_af"].sum()
        off = cell_sum[(cell_sum < 0.85) | (cell_sum > 1.1)]
        if len(off):
            print(
                "warning: reference frequencies do not sum to 1 for "
                + ", ".join(f"{gene} {group} ({value:.2f})" for (group, gene), value in off.items())
                + "; the AFND rows behind those cells are suspect",
                file=sys.stderr,
            )
    return reference, populations_out


# --- correlation ----------------------------------------------------------


def _pearson(x, y) -> tuple[float, float, float]:
    """Pearson r, R², and OLS slope of y on x. NaN when fewer than 3 points."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 3:
        return np.nan, np.nan, np.nan
    x_center = x - x.mean()
    y_center = y - y.mean()
    ss_x = float(np.dot(x_center, x_center))
    ss_y = float(np.dot(y_center, y_center))
    if ss_x == 0 or ss_y == 0:
        return np.nan, np.nan, np.nan
    cross = float(np.dot(x_center, y_center))
    r = cross / np.sqrt(ss_x * ss_y)
    return r, r * r, cross / ss_x


def _summary_row(gene: str, group: str, sub: pd.DataFrame, populations: pd.DataFrame) -> dict:
    r, r2, slope = _pearson(sub["afnd_af"], sub["local_af"])
    shared = sub[(sub["local_af"] > 0) & (sub["afnd_af"] > 0)]
    r_shared, r2_shared, _ = _pearson(shared["afnd_af"], shared["local_af"])
    if gene == "(all)":
        block = populations[populations["Group"] == group]
        # A population's share of the whole reference: its per-gene weights
        # summed over every gene in this group, not averaged over only the
        # genes it typed. Averaging over the latter hands the sole source of
        # one rarely-typed gene a weight of 1.0, so a single small study
        # outranks a registry carrying most of eight genes.
        n_genes = int(block["Gene"].nunique()) if not block.empty else 0
        pops = block.groupby("population_id", as_index=False).agg(
            population=("population", "first"), weight=("weight", "sum")
        )
        if n_genes:
            pops["weight"] = pops["weight"] / n_genes
    else:
        pops = populations[(populations["Group"] == group) & (populations["Gene"] == gene)]
    top_population = ""
    top_weight = np.nan
    if not pops.empty and pops["weight"].notna().any():
        top = pops.loc[pops["weight"].idxmax()]
        top_population = str(top["population"])
        top_weight = float(top["weight"])
    local_n = sub["local_n"].dropna()
    afnd_n = sub["afnd_n"].dropna()
    afnd_pops = sub["afnd_populations"].dropna()
    pooled = gene == "(all)"
    return {
        "Gene": gene,
        "Group": group,
        "n_alleles": int(len(sub)),
        "n_shared": int(len(shared)),
        "n_local_only": int(((sub["local_af"] > 0) & (sub["afnd_af"] <= 0)).sum()),
        "n_ref_only": int(((sub["local_af"] <= 0) & (sub["afnd_af"] > 0)).sum()),
        "n_haplotypes": pd.NA if pooled or local_n.empty else int(local_n.iloc[0]),
        "n_individuals": pd.NA if pooled or afnd_n.empty else int(afnd_n.iloc[0]),
        "n_populations": pd.NA if pooled or afnd_pops.empty else int(afnd_pops.iloc[0]),
        "top_population": top_population,
        "top_weight": top_weight,
        "r": r,
        "r2": r2,
        "slope": slope,
        "r_shared": r_shared,
        "r2_shared": r2_shared,
    }


def correlation_tables(
    local: pd.DataFrame,
    reference: pd.DataFrame,
    populations: pd.DataFrame,
    *,
    min_af: float = 0.001,
    min_ref_n: int = 500,
    min_ref_pops: int = 2,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Outer-join cohort and AFND frequencies and compute Pearson R².

    An allele present on only one side is filled with zero on the other, so a
    cohort allele missing from AFND counts against the correlation. Alleles
    below ``min_af`` on both sides are left out (default 0.001), which keeps a
    long tail of registry-only singletons from dominating the correlation.
    ``r2`` is computed on that union. ``r2_shared`` is the same correlation
    restricted to alleles with a positive frequency on both sides. ``slope`` is
    the OLS slope of the cohort frequency on the AFND frequency; a matching
    spectrum has slope near 1.

    The ``(all)`` rows pool every gene in a superpopulation into one correlation.
    """
    if local.empty:
        raise AFNDError("local frequency table is empty")
    covered = set(zip(reference["Gene"], reference["Group"])) if not reference.empty else set()
    if covered and (min_ref_n > 0 or min_ref_pops > 1):
        support = reference.groupby(["Gene", "Group"])[["afnd_n", "afnd_populations"]].first()
        thin = support[
            (support["afnd_n"] < min_ref_n)
            | (
                (support["afnd_populations"] < min_ref_pops)
                & (support["afnd_n"] < LARGE_REF_N)
            )
        ]
        if len(thin):
            print(
                "note: AFND support too thin for "
                + ", ".join(
                    f"{gene} {group} ({int(row.afnd_populations)} population(s), n={int(row.afnd_n)})"
                    for (gene, group), row in thin.iterrows()
                )
                + "; those cells are omitted (--min-ref-n / --min-ref-pops)",
                file=sys.stderr,
            )
            covered -= set(thin.index)
            if not covered:
                raise AFNDError(
                    "every gene/population cell fell below --min-ref-n / --min-ref-pops; "
                    "lower those thresholds or compare genes AFND types deeply"
                )
    wanted = local[["Gene", "Group"]].drop_duplicates()
    missing = [
        f"{gene} {group}"
        for gene, group in wanted.itertuples(index=False)
        if (gene, group) not in covered
    ]
    if missing:
        print(
            "note: no AFND frequencies for " + ", ".join(missing) + "; those cells are omitted",
            file=sys.stderr,
        )
    merged = local.merge(reference, on=["Gene", "Allele", "Group"], how="outer")
    merged = merged.merge(wanted, on=["Gene", "Group"], how="inner")
    merged = merged[merged.apply(lambda row: (row["Gene"], row["Group"]) in covered, axis=1)]
    if merged.empty:
        raise AFNDError("no gene overlapped the cohort table and AFND")

    n_map = local.groupby(["Gene", "Group"])["local_n"].first()
    ref_n = reference.groupby(["Gene", "Group"])[["afnd_n", "afnd_populations"]].first() if not reference.empty else pd.DataFrame()
    merged["local_af"] = merged["local_af"].fillna(0.0)
    merged["afnd_af"] = merged["afnd_af"].fillna(0.0)
    # A frequency matrix has no haplotype counts. Fill a count only for alleles
    # that are absent from the cohort; leave the matrix rows blank.
    merged.loc[merged["local_af"].eq(0), "local_count"] = merged.loc[merged["local_af"].eq(0), "local_count"].fillna(0)
    merged["local_count"] = pd.to_numeric(merged["local_count"], errors="coerce").astype("Int64")
    merged["local_n"] = [
        pd.NA if (gene, group) not in n_map.index or pd.isna(n_map.loc[(gene, group)]) else int(n_map.loc[(gene, group)])
        for gene, group in zip(merged["Gene"], merged["Group"])
    ]
    if not ref_n.empty:
        merged["afnd_n"] = [
            int(ref_n.loc[(gene, group), "afnd_n"]) if (gene, group) in ref_n.index else pd.NA
            for gene, group in zip(merged["Gene"], merged["Group"])
        ]
        merged["afnd_populations"] = [
            int(ref_n.loc[(gene, group), "afnd_populations"]) if (gene, group) in ref_n.index else pd.NA
            for gene, group in zip(merged["Gene"], merged["Group"])
        ]
    merged.loc[merged["local_af"] < 1e-12, "local_af"] = 0.0
    merged.loc[merged["afnd_af"] < 1e-12, "afnd_af"] = 0.0
    if min_af > 0:
        merged = merged[(merged["local_af"] >= min_af) | (merged["afnd_af"] >= min_af)]
    merged = merged[(merged["local_af"] > 0) | (merged["afnd_af"] > 0)]

    gene_rank = {gene: index for index, gene in enumerate(CLASSICAL_GENES)}
    group_rank = {group: index for index, group in enumerate(SUPERPOPS)}
    merged["_gene"] = merged["Gene"].map(lambda gene: gene_rank.get(gene, 100))
    merged["_group"] = merged["Group"].map(lambda group: group_rank.get(group, 100))
    merged = merged.sort_values(
        ["_group", "_gene", "local_af", "afnd_af"],
        ascending=[True, True, False, False],
    ).drop(columns=["_gene", "_group"])
    comparison = merged[
        [
            "Gene", "Allele", "Group", "local_af", "afnd_af",
            "local_count", "local_n", "afnd_n", "afnd_populations",
        ]
    ].reset_index(drop=True)

    rows = []
    for group in [group for group in SUPERPOPS if group in set(comparison["Group"])]:
        block = comparison[comparison["Group"] == group]
        rows.append(_summary_row("(all)", group, block, populations))
        for gene in CLASSICAL_GENES + sorted(set(block["Gene"]) - set(CLASSICAL_GENES)):
            sub = block[block["Gene"] == gene]
            if sub.empty:
                continue
            rows.append(_summary_row(gene, group, sub, populations))
    summary = pd.DataFrame(rows)
    return comparison, summary


# --- plot -----------------------------------------------------------------


def _p_math(p: float) -> str:
    """Pearson p-value as mathtext, in scientific notation below 0.001."""
    if not np.isfinite(p):
        return r"\mathrm{n/a}"
    if p <= 0:
        return r"< 10^{-300}"
    if p < 1e-3:
        exp = int(np.floor(np.log10(p)))
        mant = p / 10.0**exp
        return f"{mant:.1f} \\times 10^{{{exp}}}"
    return f"{p:.3f}"


def _reg_label(group: str, sub: pd.DataFrame, stats) -> str:
    """``y = kx + b``, R², and the two-sided Pearson p for one superpopulation."""
    x = np.asarray(sub["afnd_af"], dtype=float)
    y = np.asarray(sub["local_af"], dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 3 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return group
    slope, intercept = (float(value) for value in np.polyfit(x, y, 1))
    r, p = stats.pearsonr(x, y)
    sign = "+" if intercept >= 0 else "-"
    return (
        f"{group}   $y = {slope:.2f}x {sign} {abs(intercept):.3f}$,  "
        f"$R^2 = {r * r:.2f}$,  $p = {_p_math(float(p))}$"
    )


def plot_regplot(
    comparison: pd.DataFrame,
    summary: pd.DataFrame,
    path: str | Path | None = None,
    *,
    figsize: tuple[float, float] | None = None,
):
    """One axes, with an ``sns.regplot`` line for every superpopulation.

    Alleles from every gene are pooled within a superpopulation, which is the
    same set of points as the ``(all)`` row of the summary. The line is the
    OLS fit of cohort frequency on AFND frequency, with a 95% confidence
    band. The legend, inside the axes, gives that line as ``y = kx + b``
    with its R² and two-sided Pearson p-value. ``figsize`` is in inches
    and defaults to (7.6, 5.6).
    """
    del summary  # the line is fit to `comparison`
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.lines import Line2D
    from scipy import stats

    groups = [group for group in SUPERPOPS if group in set(comparison["Group"])]
    extra = sorted(set(comparison["Group"]) - set(groups))
    groups.extend(extra)
    if not groups:
        raise AFNDError("nothing to plot")

    hi = float(max(comparison["local_af"].max(), comparison["afnd_af"].max(), 0.01)) * 1.08
    fig, ax = plt.subplots(figsize=figsize or (7.6, 5.6))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    ax.plot([0, hi], [0, hi], color="#c6c5bd", linewidth=1.0, zorder=1)
    ax.set_xlim(0, hi)
    ax.set_ylim(0, hi)

    handles = []
    for group in groups:
        sub = comparison[comparison["Group"] == group]
        color = GROUP_COLORS.get(group, OTHER_COLOR)
        if len(sub) >= 3 and sub["afnd_af"].nunique() > 1 and sub["local_af"].nunique() > 1:
            scatter_kws = {"s": 22, "alpha": 0.45, "linewidths": 0, "zorder": 3}
            sns.regplot(
                data=sub,
                x="afnd_af",
                y="local_af",
                ax=ax,
                color=color,
                ci=95,
                truncate=False,
                scatter_kws=scatter_kws,
                line_kws={"linewidth": 2.0, "zorder": 4},
            )
        else:
            ax.scatter(
                sub["afnd_af"],
                sub["local_af"],
                s=22,
                color=color,
                alpha=0.45,
                linewidths=0,
                zorder=3,
            )
        handles.append(Line2D([0], [0], color=color, lw=2, label=_reg_label(group, sub, stats)))

    # regplot resets the limits to the points; put the shared diagonal back.
    ax.set_xlim(0, hi)
    ax.set_ylim(0, hi)
    ax.legend(
        handles=handles,
        loc="upper left",
        frameon=False,
        fontsize=8,
        labelcolor="black",
        handlelength=1.4,
        borderaxespad=0.4,
    )
    ax.set_xlabel("AFND allele frequency", fontsize=9, color="black")
    ax.set_ylabel("Cohort allele frequency", fontsize=9, color="black")
    ax.tick_params(labelsize=8, colors="black", length=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")
    ax.xaxis.set_major_formatter(lambda value, _: f"{value:.0%}")
    ax.yaxis.set_major_formatter(lambda value, _: f"{value:.0%}")
    fig.tight_layout()
    if path is not None:
        path = Path(path).expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    return fig


def plot_comparison(
    comparison: pd.DataFrame,
    summary: pd.DataFrame,
    path: str | Path | None = None,
    *,
    regplot: bool = False,
    figsize: tuple[float, float] | None = None,
):
    """Scatter of cohort vs AFND frequency, one panel per superpopulation.

    ``regplot=True`` draws every superpopulation on one axes instead, each
    with an ``sns.regplot`` line and 95% confidence band. Returns the
    matplotlib Figure so a notebook can render it. Pass ``path`` to also
    write a PNG. ``figsize`` is in inches; the faceted default is
    (4, 3.6) per panel, and the regplot default is (7.6, 5.6).
    """
    if regplot:
        return plot_regplot(comparison, summary, path, figsize=figsize)

    import matplotlib.pyplot as plt

    groups = [group for group in SUPERPOPS if group in set(comparison["Group"])]
    if not groups:
        raise AFNDError("nothing to plot")
    ncols = 3
    nrows = int(np.ceil(len(groups) / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=figsize or (4.0 * ncols, 3.6 * nrows), squeeze=False
    )
    fig.patch.set_facecolor(SURFACE)

    for axis, group in zip(axes.ravel(), groups):
        sub = comparison[comparison["Group"] == group]
        for gene, part in sub.groupby("Gene", sort=False):
            color = GENE_COLORS.get(gene, OTHER_COLOR)
            axis.scatter(
                part["afnd_af"],
                part["local_af"],
                s=18,
                color=color,
                alpha=0.9,
                linewidths=0,
                zorder=3,
                label=gene,
            )
        hi = float(max(sub["local_af"].max(), sub["afnd_af"].max(), 0.01)) * 1.08
        axis.plot([0, hi], [0, hi], color="#c6c5bd", linewidth=1.0, zorder=1)
        axis.set_xlim(0, hi)
        axis.set_ylim(0, hi)
        axis.set_title(group, loc="left", fontsize=11, color="black", pad=8)
        axis.tick_params(labelsize=8, colors="black", length=0)
        for side in ("top", "right"):
            axis.spines[side].set_visible(False)
        axis.spines["left"].set_color("black")
        axis.spines["bottom"].set_color("black")
        axis.set_facecolor(SURFACE)
        axis.xaxis.set_major_formatter(lambda value, _: f"{value:.0%}")
        axis.yaxis.set_major_formatter(lambda value, _: f"{value:.0%}")

    for axis in axes.ravel()[len(groups):]:
        axis.set_visible(False)

    fig.supxlabel("AFND allele frequency", fontsize=9, color="black")
    fig.supylabel("Cohort allele frequency", fontsize=9, color="black")
    handles, labels = [], []
    seen = set()
    for axis in axes.ravel():
        for handle, label in zip(*axis.get_legend_handles_labels()):
            if label not in seen:
                handles.append(handle)
                labels.append(label)
                seen.add(label)
    # Stable gene order.
    order = {gene: index for index, gene in enumerate(CLASSICAL_GENES)}
    paired = sorted(zip(labels, handles), key=lambda item: order.get(item[0], 100))
    if paired:
        labels, handles = zip(*paired)
        fig.legend(
            handles,
            labels,
            loc="upper left",
            bbox_to_anchor=(0.01, 1.0),
            ncol=min(len(labels), 8),
            frameon=False,
            fontsize=8,
            handlelength=0.8,
            handletextpad=0.4,
        )
    fig.tight_layout(rect=(0.015, 0.02, 1, 0.93))
    if path is not None:
        path = Path(path).expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=200, facecolor=SURFACE)
    return fig


# --- driver ---------------------------------------------------------------


def compare(
    freq: str | Path | pd.DataFrame,
    *,
    preset: str = "registry",
    standard: str = "gold",
    resolution: int = 2,
    weight: str = "sample",
    min_af: float = 0.001,
    min_ref_n: int = 500,
    min_ref_pops: int = 2,
    groups: list[str] | None = None,
    cache: str | Path | None = None,
    refresh: bool = False,
    verbose: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Cohort-vs-AFND comparison, R² summary, and population weights.

    ``freq`` is a path or DataFrame from ``hla_allele_frequency.py``. The
    download is cached at ``cache`` when that path is set.
    """
    local = load_local_frequencies(freq, resolution=resolution, groups=groups)
    present = [group for group in SUPERPOPS if group in set(local["Group"])]
    skipped = sorted(gene for gene in set(local["Gene"]) if not str(gene).startswith("HLA-"))
    if skipped:
        print(
            "note: AFND HLA search has no locus for " + ", ".join(skipped) + "; those genes are omitted",
            file=sys.stderr,
        )
        local = local[local["Gene"].str.startswith("HLA-")].copy()
        if local.empty:
            raise AFNDError("no HLA genes left to compare")
    genes = [gene for gene in CLASSICAL_GENES if gene in set(local["Gene"])]
    genes += sorted(set(local["Gene"]) - set(genes))
    raw = fetch_afnd(
        genes,
        present,
        preset=preset,
        standard=standard,
        cache=cache,
        refresh=refresh,
        verbose=verbose,
    )
    reference, populations = aggregate_afnd(raw, resolution=resolution, weight=weight)
    comparison, summary = correlation_tables(
        local,
        reference,
        populations,
        min_af=min_af,
        min_ref_n=min_ref_n,
        min_ref_pops=min_ref_pops,
    )
    return comparison, summary, populations


class AFNDComparison:
    """Cohort frequencies against Allele Frequency Net, for a notebook.

    ``fit.summary`` is the R² table, ``fit.comparison`` the per-allele
    frequencies, and ``fit.populations`` the weight each gold population
    contributed. ``fit.plot()`` returns a matplotlib Figure.
    """

    def __init__(self, comparison: pd.DataFrame, summary: pd.DataFrame, populations: pd.DataFrame):
        self.comparison = comparison
        self.summary = summary
        self.populations = populations

    @classmethod
    def from_table(cls, freq: str | Path | pd.DataFrame, **kwargs) -> "AFNDComparison":
        """Build the comparison from a frequency TSV or DataFrame.

        ``freq`` is the TSV ``hla_allele_frequency.py`` writes,
        ``AlleleFrequencies.table``, or ``af.frequency_matrix()``. Keyword
        arguments are the ones ``compare`` takes: preset, standard, resolution,
        weight, min_af, min_ref_n, min_ref_pops, groups, cache, refresh,
        verbose.
        """
        comparison, summary, populations = compare(freq, **kwargs)
        return cls(comparison, summary, populations)

    def plot(
        self,
        path: str | Path | None = None,
        *,
        regplot: bool = False,
        figsize: tuple[float, float] | None = None,
    ):
        """Scatter of cohort vs AFND frequency. Pass ``path`` to also save a PNG.

        ``regplot=True`` puts every superpopulation on one axes, with an
        ``sns.regplot`` line, a 95% confidence band, and a legend giving
        ``y = kx + b``, R², and the Pearson p-value. ``figsize`` is
        ``(width, height)`` in inches.
        """
        return plot_comparison(
            self.comparison, self.summary, path, regplot=regplot, figsize=figsize
        )

    def save(
        self,
        out: str | Path,
        *,
        regplot: bool = False,
        figsize: tuple[float, float] | None = None,
    ) -> None:
        """Write the comparison, the R² summary, the population weights, and the plot.

        Paths sit next to ``out``, the same way the command line does:
        ``hla_afnd_r2.tsv``, ``hla_afnd_r2_summary.tsv``,
        ``hla_afnd_r2_populations.tsv``, ``hla_afnd_r2.png``.
        """
        out = Path(out).expanduser()
        stem = out.with_suffix("")
        _write_table(self.comparison, out)
        _write_table(self.summary, Path(str(stem) + "_summary.tsv"))
        _write_table(self.populations, Path(str(stem) + "_populations.tsv"))
        self.plot(stem.with_suffix(".png"), regplot=regplot, figsize=figsize)


def _write_table(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, float_format="%.6g")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--freq", required=True, type=Path, help="TSV from hla_allele_frequency.py")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("hla_afnd_r2.tsv"),
        help="per-allele comparison TSV (summary, populations, plot, and cache are named from this)",
    )
    parser.add_argument("--preset", choices=sorted(PRESETS), default="registry")
    parser.add_argument("--standard", choices=sorted(STANDARDS), default="gold")
    parser.add_argument("--resolution", type=int, default=2, help="allele fields to compare (default 2)")
    parser.add_argument("--weight", choices=["sample", "equal"], default="sample")
    parser.add_argument(
        "--min-af",
        type=float,
        default=0.001,
        help="drop alleles below this frequency on both sides (default 0.001; 0 keeps every allele)",
    )
    parser.add_argument(
        "--min-ref-n",
        type=int,
        default=500,
        help="drop a gene/population cell whose AFND reference totals fewer individuals (default 500; 0 disables)",
    )
    parser.add_argument(
        "--min-ref-pops",
        type=int,
        default=2,
        help=f"drop a gene/population cell backed by fewer AFND populations, unless it has "
        f"{LARGE_REF_N}+ typed individuals (default 2; 1 disables)",
    )
    parser.add_argument(
        "--groups",
        nargs="+",
        help="superpopulations to compare (default: whichever of the six are in the table)",
    )
    parser.add_argument("--cache", type=Path, help="AFND cache TSV (default: next to --out)")
    parser.add_argument("--refresh", action="store_true", help="redownload queries already in the cache")
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument(
        "--regplot",
        action="store_true",
        help="draw every superpopulation on one axes with an sns.regplot line",
    )
    parser.add_argument(
        "--figsize",
        nargs=2,
        type=float,
        metavar=("WIDTH", "HEIGHT"),
        help="figure size in inches (default 7.6 5.6 for --regplot, else one 4 by 3.6 panel)",
    )
    parser.add_argument("--summary-out", type=Path, help="R² summary TSV")
    parser.add_argument("--populations-out", type=Path, help="population-weight TSV")
    parser.add_argument("--plot", type=Path, help="scatter plot path")
    args = parser.parse_args(argv)

    out = args.out.expanduser()
    stem = out.with_suffix("")
    summary_path = args.summary_out or Path(str(stem) + "_summary.tsv")
    populations_path = args.populations_out or Path(str(stem) + "_populations.tsv")
    plot_path = args.plot or stem.with_suffix(".png")
    cache_path = args.cache or stem.parent / "hla_afnd_cache.tsv"
    groups = [group.upper() for group in args.groups] if args.groups else None
    if groups:
        bad = [group for group in groups if group not in SUPERPOPS]
        if bad:
            print(f"error: unknown group(s) {', '.join(bad)}; expected {', '.join(SUPERPOPS)}", file=sys.stderr)
            return 1

    try:
        comparison, summary, populations = compare(
            args.freq,
            preset=args.preset,
            standard=args.standard,
            resolution=args.resolution,
            weight=args.weight,
            min_af=args.min_af,
            min_ref_n=args.min_ref_n,
            min_ref_pops=args.min_ref_pops,
            groups=groups,
            cache=cache_path,
            refresh=args.refresh,
        )
    except AFNDError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    _write_table(comparison, out)
    _write_table(summary, summary_path)
    _write_table(populations, populations_path)
    if not args.no_plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        figure = plot_comparison(
            comparison,
            summary,
            plot_path,
            regplot=args.regplot,
            figsize=tuple(args.figsize) if args.figsize else None,
        )
        plt.close(figure)

    shown = summary.copy()
    for column in ("r", "r2", "slope", "r_shared", "r2_shared", "top_weight"):
        shown[column] = shown[column].map(lambda value: "" if pd.isna(value) else f"{value:.4f}")
    print(shown.to_string(index=False))
    print(f"\nwrote {out}", file=sys.stderr)
    print(f"wrote {summary_path}", file=sys.stderr)
    print(f"wrote {populations_path}", file=sys.stderr)
    if not args.no_plot:
        print(f"wrote {plot_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
