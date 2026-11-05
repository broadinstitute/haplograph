#!/usr/bin/env python3
"""Classical HLA allele frequencies from Immuannot calls, stacked by population.

Importable from a notebook, and runnable as a CLI.

Notebook use
------------
    import sys; sys.path.append("scripts")
    from hla_allele_frequency import AlleleFrequencies

    af = AlleleFrequencies.from_files("hprc.immuannot.csv", "metadata.csv")
    af.table.head()                 # tidy per-(gene, allele, group) frequencies
    af.plot()                       # -> matplotlib Figure, renders inline
    af.plot(genes=["HLA-DRB1"], top_n=12, mode="composition")
    af.top_alleles("HLA-B", n=5)    # wide allele x group frame
    af.frequency_matrix()           # alleles x populations, for clustering

    # every gene in the table — non-classical HLA and KIR included — with
    # 2-field HLA alleles and 3-field KIR alleles
    af = AlleleFrequencies.from_files(
        "hprc.immuannot.csv", "metadata.csv",
        genes="all", resolution={"HLA": 2, "KIR": 3},
    )
    af.plot(ncols=4)

Or go straight through the functional layer:

    from hla_allele_frequency import load_calls, allele_frequencies, plot_frequencies
    calls = load_calls("hprc.immuannot.csv")           # one row per haplotype call
    freq  = allele_frequencies(calls, "metadata.csv", resolution=2)
    fig   = plot_frequencies(freq, top_n=10)

Everything accepts either paths or already-loaded DataFrames, so you can swap in
your own metadata or a pre-filtered call set without touching this file.

Metadata
--------
Sample ids come from `research_id` and population labels from `ancestry_pred`,
which already holds superpopulation codes (afr, amr, eas, eur, sas, mid, oth);
they are upper-cased and used as they are. Override with `sample_col=` /
`pop_col=` for a table that names those columns differently — finer 1000G/HPRC
population codes are rolled up to superpopulations automatically, and
`group="population"` keeps whatever labels the column holds.

Input layout
------------
An Immuannot per-haplotype allele table (the `hprc.immuannot.csv` layout), which
is tab-separated despite the `.csv` name:

    <index>  Sample  Hap  Gene  Strand  Consensus_Allele  Contig  Spos  Epos  Tool

Each row is one called allele on one assembled haplotype, so allele counts are
haplotype counts and frequencies are haplotype frequencies.

Genes and allele resolution
---------------------------
`genes="classical"` (default) is the eight classical HLA genes; `genes="all"`
takes every gene in the table, including non-classical HLA (DPB2, DOA, MICA, …)
and KIR, ordered classical HLA → other HLA → KIR. `resolution` is a field count
— HLA allows 1-4 fields, KIR 1-3 — given either as one int for everything
(clamped per family) or as `{"HLA": 2, "KIR": 3}`. Novelty flags are stripped
before trimming, so `HLA-DPB2*03:01:01:new` counts as `HLA-DPB2*03:01:01`;
`exclude_novel=True` drops those calls outright instead.

Stacking modes
--------------
  contribution (default)  segment = count(allele, group) / called haplotypes for
                          that gene; the full bar height is the overall allele
                          frequency, and the segments show which groups carry it.
  composition             segment = count(allele, group) / count(allele); every
                          bar reaches 1.0 and shows the population make-up of
                          each allele.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Patch

__all__ = [
    "AlleleFrequencies",
    "CLASSICAL_GENES",
    "MAX_RESOLUTION",
    "SAMPLE_COL",
    "POP_COL",
    "SUPERPOP",
    "allele_frequencies",
    "gene_family",
    "load_calls",
    "load_metadata",
    "order_genes",
    "plot_frequencies",
    "resolution_for",
    "select_genes",
    "to_superpop",
    "trim_allele",
]

# --- defaults -------------------------------------------------------------

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

# Metadata column names. `ancestry_pred` already holds superpopulation labels
# (afr, amr, eas, eur, sas, …), which are upper-cased and used as-is.
SAMPLE_COL = "research_id"
POP_COL = "ancestry_pred"

# Finer 1000G / HPRC population codes -> superpopulation, for metadata tables
# that carry those instead (e.g. HPRC's `population_abbreviation`). Labels that
# are already superpopulations pass straight through — see `to_superpop`.
SUPERPOP = {
    # African
    "ACB": "AFR", "ASW": "AFR", "ESN": "AFR", "GWD": "AFR", "LWK": "AFR",
    "MSL": "AFR", "YRI": "AFR",
    # Admixed American
    "CLM": "AMR", "MXL": "AMR", "PEL": "AMR", "PUR": "AMR",
    # East Asian
    "CDX": "EAS", "CHB": "EAS", "CHS": "EAS", "JPT": "EAS", "KHV": "EAS",
    # European
    "CEU": "EUR", "FIN": "EUR", "GBR": "EUR", "IBS": "EUR", "TSI": "EUR",
    # South Asian
    "BEB": "SAS", "GIH": "SAS", "ITU": "SAS", "PJL": "SAS", "STU": "SAS",
}

# Field counts the nomenclature allows: HLA up to 4, KIR up to 3.
MAX_RESOLUTION = {"HLA": 4, "KIR": 3}

# Novelty flags Immuannot appends (`*03:01:01:new`, `*01:new`) are stripped.
_NEW_RE = re.compile(r"new", re.IGNORECASE)

# MID (Middle Eastern) and OTH (Other) appear in All of Us ancestry calls.
SUPERPOP_ORDER = ["AFR", "AMR", "EAS", "EUR", "SAS", "MID", "OTH"]
UNKNOWN = "Unknown"


def to_superpop(code) -> str | float:
    """Superpopulation for a population label; already-superpop codes pass through."""
    if not isinstance(code, str):
        return pd.NA
    code = code.strip().upper()
    if code in SUPERPOP_ORDER:
        return code
    return SUPERPOP.get(code, pd.NA)

# Validated categorical palette, fixed order — never cycled.
SERIES_COLORS = [
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
    "#e87ba4", "#6250d6", "#008300", "#e34948",
]
OTHER_COLOR = "#898781"
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_MUTED = "#898781"
GRID = "#e1e0d9"

PathOrFrame = "str | Path | pd.DataFrame"


class HLAFrequencyError(ValueError):
    """Raised for malformed inputs, so notebooks get an exception, not a sys.exit."""


# --- loading --------------------------------------------------------------


def _as_frame(obj, reader) -> pd.DataFrame:
    if isinstance(obj, pd.DataFrame):
        return obj.copy()
    return reader(Path(obj).expanduser())


def load_calls(annot: PathOrFrame) -> pd.DataFrame:
    """Read an Immuannot allele table (sniffing tab vs comma) into a DataFrame."""

    def _read(path: Path) -> pd.DataFrame:
        return pd.read_csv(
            path, sep=None, engine="python", dtype=str, encoding_errors="replace"
        )

    df = _as_frame(annot, _read)
    df.columns = [str(c).strip() for c in df.columns]
    required = {"Sample", "Hap", "Gene", "Consensus_Allele"}
    missing = required - set(df.columns)
    if missing:
        raise HLAFrequencyError(f"allele table is missing column(s): {', '.join(sorted(missing))}")
    for col in ("Sample", "Hap", "Gene", "Consensus_Allele"):
        df[col] = df[col].astype(str).str.strip()
    return df.dropna(subset=["Consensus_Allele"])


def load_metadata(
    metadata: PathOrFrame,
    sample_col: str = SAMPLE_COL,
    pop_col: str = POP_COL,
) -> pd.DataFrame:
    """Read sample metadata into a two-column `Sample` / `Population` frame.

    Defaults to the `research_id` / `ancestry_pred` layout. Pass a dict or a
    `{Sample: Population}` Series straight to `allele_frequencies` instead if you
    already have a mapping in the notebook.
    """

    def _read(path: Path) -> pd.DataFrame:
        # Some metadata releases carry stray non-UTF-8 bytes in free-text
        # columns; they never matter for the columns used here.
        return pd.read_csv(path, dtype=str, encoding="utf-8", encoding_errors="replace")

    df = _as_frame(metadata, _read)
    df.columns = [str(c).strip() for c in df.columns]
    for col in (sample_col, pop_col):
        if col not in df.columns:
            raise HLAFrequencyError(
                f"metadata has no column '{col}'. Available: {', '.join(df.columns)}"
            )
    out = df[[sample_col, pop_col]].copy()
    out.columns = ["Sample", "Population"]
    out["Sample"] = out["Sample"].astype(str).str.strip()
    out["Population"] = out["Population"].astype(str).str.strip().str.upper()
    # Some metadata releases leave free-text descriptors or "N/A" in the
    # population column; keep only short codes.
    out.loc[~out["Population"].str.fullmatch(r"[A-Z]{2,4}", na=False), "Population"] = pd.NA
    out.loc[out["Population"].isin({"NA", "NAN", "NONE"}), "Population"] = pd.NA
    return out.drop_duplicates("Sample")


def _coerce_metadata(metadata) -> pd.DataFrame:
    """Accept a path, a DataFrame, a dict, or a Series of Sample -> Population."""
    if isinstance(metadata, pd.Series):
        return metadata.rename("Population").rename_axis("Sample").reset_index()
    if isinstance(metadata, dict):
        return pd.DataFrame(
            {"Sample": list(metadata), "Population": list(metadata.values())}
        )
    return metadata


def gene_family(gene: str) -> str:
    """'KIR' for KIR genes, 'HLA' for everything else (HLA, MIC, TAP, …)."""
    return "KIR" if str(gene).upper().startswith("KIR") else "HLA"


def _split_fields(fields: str, family: str) -> list[str]:
    """Allele fields, normalising the compact IPD-KIR form to colon-separated.

    IPD-KIR writes alleles as `KIR2DL1*0030101` (3-digit family, then 2-digit
    fields, with an optional expression letter); Immuannot may emit either that
    or the colon form. Both end up as ['003', '01', '01'].
    """
    parts = [p for p in fields.split(":") if p != ""]
    if family != "KIR" or len(parts) != 1:
        return parts
    token, suffix = parts[0], ""
    while token and token[-1].isalpha():
        token, suffix = token[:-1], token[-1] + suffix
    if not (token.isdigit() and len(token) >= 3):
        return parts
    chunks = [token[:3]] + [token[i : i + 2] for i in range(3, len(token), 2)]
    if suffix:
        chunks[-1] += suffix
    return chunks


def resolution_for(gene: str, resolution: int | dict[str, int] = 2) -> int:
    """Field count to keep for one gene.

    `resolution` is either an int applied to every gene (clamped to the family
    maximum) or a per-family dict such as `{"HLA": 2, "KIR": 3}`, which is
    validated strictly. HLA alleles have up to 4 fields, KIR alleles up to 3.
    """
    family = gene_family(gene)
    cap = MAX_RESOLUTION[family]
    if isinstance(resolution, dict):
        n = resolution.get(family, resolution.get(family.lower(), 2))
        if not 1 <= int(n) <= cap:
            raise HLAFrequencyError(f"{family} resolution must be 1-{cap}, got {n}")
        return int(n)
    if int(resolution) < 1:
        raise HLAFrequencyError(f"resolution must be >= 1, got {resolution}")
    return min(int(resolution), cap)


def trim_allele(allele: str, resolution: int | dict[str, int] = 2) -> str | None:
    """`HLA-A*01:01:01:01` -> `HLA-A*01:01` at resolution 2.

    Novelty flags are stripped first, so `HLA-DPB2*03:01:01:new` is treated as
    `HLA-DPB2*03:01:01` and `HLA-DPA2*01:new` as `HLA-DPA2*01`.
    """
    if not allele or "*" not in allele:
        return None
    gene, _, fields = allele.partition("*")
    fields = _NEW_RE.sub("", fields)
    family = gene_family(gene)
    parts = [p for p in _split_fields(fields, family) if p != ""]
    if not parts:
        return None
    return f"{gene}*{':'.join(parts[: resolution_for(gene, resolution)])}"


def order_genes(genes: Iterable[str]) -> list[str]:
    """Classical HLA first, then other HLA genes, then KIR, each A-Z."""
    genes = list(dict.fromkeys(genes))
    classical = [g for g in CLASSICAL_GENES if g in genes]
    rest = [g for g in genes if g not in classical]
    other_hla = sorted(g for g in rest if gene_family(g) == "HLA")
    kir = sorted(g for g in rest if gene_family(g) == "KIR")
    return classical + other_hla + kir


def select_genes(available: Iterable[str], genes="classical") -> list[str]:
    """Resolve a `genes=` argument against the genes actually present.

    `"all"` (or None) takes every gene in the table, `"classical"` the eight
    classical HLA genes, or pass an explicit list.
    """
    available = list(dict.fromkeys(available))
    if genes is None or (isinstance(genes, str) and genes.lower() == "all"):
        return order_genes(available)
    if isinstance(genes, str) and genes.lower() == "classical":
        return [g for g in CLASSICAL_GENES if g in available]
    wanted = [genes] if isinstance(genes, str) else list(genes)
    if len(wanted) == 1 and wanted[0].lower() in ("all", "classical"):
        return select_genes(available, wanted[0].lower())
    return [g for g in wanted if g in available]


# --- frequency table ------------------------------------------------------


def allele_frequencies(
    calls: PathOrFrame,
    metadata,
    *,
    genes: Sequence[str] | str | None = "classical",
    resolution: int | dict[str, int] = 2,
    group: str = "superpop",
    exclude_novel: bool = False,
    drop_unknown_pop: bool = False,
    sample_col: str = SAMPLE_COL,
    pop_col: str = POP_COL,
    verbose: bool = True,
) -> pd.DataFrame:
    """Tidy per-(gene, allele, group) counts and frequencies.

    Columns: Gene, Allele, Group, count, gene_total, allele_total, freq,
    within_allele. `freq` is the group's contribution to the gene-wide allele
    frequency; `within_allele` is that group's share of the allele's carriers.

    `genes` takes "classical" (default), "all" for every gene in the table —
    non-classical HLA and KIR included — or an explicit list. `resolution` is a
    field count (HLA 1-4, KIR 1-3) or a per-family dict, e.g. {"HLA": 2, "KIR": 3}.
    """
    if group not in ("superpop", "population"):
        raise HLAFrequencyError("group must be 'superpop' or 'population'")

    df = load_calls(calls)

    wanted = select_genes(df["Gene"].unique(), genes)
    if not wanted:
        raise HLAFrequencyError(
            "none of the requested genes are in the table. Present: "
            + ", ".join(order_genes(df["Gene"].unique())[:40])
        )
    if exclude_novel:
        df = df[~df["Consensus_Allele"].str.contains("new", case=False, na=False)]
    df = df[df["Gene"].isin(wanted)].copy()
    if df.empty:
        raise HLAFrequencyError("no rows left for the requested genes")

    df["Allele"] = [trim_allele(a, resolution) for a in df["Consensus_Allele"]]
    df = df.dropna(subset=["Allele"])

    meta = _coerce_metadata(metadata)
    if not (isinstance(meta, pd.DataFrame) and set(meta.columns) == {"Sample", "Population"}):
        meta = load_metadata(meta, sample_col, pop_col)
    df = df.merge(meta, on="Sample", how="left")

    df["Group"] = df["Population"].map(to_superpop) if group == "superpop" else df["Population"]

    unknown = int(df["Group"].isna().sum())
    if unknown:
        if drop_unknown_pop:
            if verbose:
                print(f"note: dropped {unknown} haplotype calls with no population label", file=sys.stderr)
            df = df.dropna(subset=["Group"])
        else:
            if verbose:
                print(f"note: {unknown} haplotype calls have no population label — grouped as '{UNKNOWN}'", file=sys.stderr)
            df["Group"] = df["Group"].fillna(UNKNOWN)

    counts = (
        df.groupby(["Gene", "Allele", "Group"], dropna=False).size().rename("count").reset_index()
    )
    gene_totals = df.groupby("Gene").size().rename("gene_total")
    allele_totals = df.groupby(["Gene", "Allele"]).size().rename("allele_total")
    counts = counts.merge(gene_totals, on="Gene").merge(allele_totals, on=["Gene", "Allele"])
    counts["freq"] = counts["count"] / counts["gene_total"]
    counts["within_allele"] = counts["count"] / counts["allele_total"]
    # Keep the gene order meaningful (classical HLA, other HLA, then KIR).
    gene_rank = {g: i for i, g in enumerate(order_genes(wanted))}
    counts["_gene_rank"] = counts["Gene"].map(gene_rank)
    counts = counts.sort_values(
        ["_gene_rank", "allele_total", "Group"], ascending=[True, False, True]
    ).drop(columns="_gene_rank")
    counts.attrs["group_by"] = group
    counts.attrs["resolution"] = resolution
    counts.attrs["gene_order"] = order_genes(counts["Gene"])
    return counts


def group_order(groups: Iterable[str], group_by: str = "superpop") -> list[str]:
    groups = list(groups)
    if group_by == "superpop":
        ranked = [g for g in SUPERPOP_ORDER if g in groups]
        ranked += sorted(g for g in groups if g not in SUPERPOP_ORDER and g != UNKNOWN)
    else:
        ranked = sorted(g for g in groups if g != UNKNOWN)
    if UNKNOWN in groups:
        ranked.append(UNKNOWN)
    return ranked


def group_colors(groups: Sequence[str]) -> dict[str, str]:
    """Colors follow the group identity, never its rank in a given panel."""
    return {
        g: OTHER_COLOR if (g == UNKNOWN or i >= len(SERIES_COLORS)) else SERIES_COLORS[i]
        for i, g in enumerate(groups)
    }


def top_alleles(
    freq: pd.DataFrame, gene: str, n: int = 10, mode: str = "contribution"
) -> pd.DataFrame:
    """Wide allele x group frame for one gene, highest frequency first."""
    value_col = "freq" if mode == "contribution" else "within_allele"
    sub = freq[freq["Gene"] == gene]
    if sub.empty:
        raise HLAFrequencyError(f"no calls for {gene}")
    keep = list(sub.groupby("Allele")["allele_total"].first().sort_values(ascending=False).head(n).index)
    groups = group_order(freq["Group"].unique(), freq.attrs.get("group_by", "superpop"))
    return (
        sub[sub["Allele"].isin(keep)]
        .pivot_table(index="Allele", columns="Group", values=value_col, fill_value=0.0)
        .reindex(index=keep, columns=groups, fill_value=0.0)
    )


def frequency_matrix(freq: pd.DataFrame, mode: str = "contribution") -> pd.DataFrame:
    """Every allele x group as one matrix — handy for clustering or heatmaps."""
    value_col = "freq" if mode == "contribution" else "within_allele"
    groups = group_order(freq["Group"].unique(), freq.attrs.get("group_by", "superpop"))
    return freq.pivot_table(
        index=["Gene", "Allele"], columns="Group", values=value_col, fill_value=0.0
    ).reindex(columns=groups, fill_value=0.0)


# --- plotting -------------------------------------------------------------


def plot_gene(
    freq: pd.DataFrame,
    gene: str,
    ax: Axes,
    *,
    top_n: int = 10,
    mode: str = "contribution",
    colors: dict[str, str] | None = None,
    label_bars: bool = True,
) -> Axes:
    """Draw one gene's stacked bars onto an existing Axes."""
    wide = top_alleles(freq, gene, top_n, mode)
    groups = list(wide.columns)
    colors = colors or group_colors(groups)
    keep = list(wide.index)
    y = list(range(len(keep) - 1, -1, -1))

    left = [0.0] * len(keep)
    for g in groups:
        vals = wide[g].tolist()
        ax.barh(
            y,
            vals,
            left=left,
            height=0.62,
            color=colors.get(g, OTHER_COLOR),
            # 2px surface gap between touching segments.
            edgecolor=SURFACE,
            linewidth=1.0,
            zorder=3,
        )
        left = [l + v for l, v in zip(left, vals)]

    n_hap = int(freq.loc[freq["Gene"] == gene, "gene_total"].iloc[0])
    ax.set_yticks(y)
    ax.set_yticklabels([a.split("*", 1)[1] for a in keep], fontsize=9, color=INK)
    # Shared y-extent so bar thickness is identical across panels, even where a
    # gene has fewer than top_n alleles.
    ax.set_ylim(len(keep) - top_n - 0.7, len(keep) - 0.3)
    ax.set_title(f"{gene}  ·  n = {n_hap} haplotypes", fontsize=10, color=INK, loc="left", pad=8)
    ax.tick_params(axis="x", labelsize=9, colors=INK_MUTED, length=0)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 1.0 if mode == "composition" else max(0.05, max(left) * 1.18))
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0%}")
    ax.grid(axis="x", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.set_facecolor(SURFACE)

    if label_bars and mode == "contribution":
        pad = ax.get_xlim()[1] * 0.012
        for yi, total in zip(y, left):
            ax.text(total + pad, yi, f"{total:.1%}", va="center", fontsize=8, color=INK_MUTED)
    return ax


def plot_frequencies(
    freq: pd.DataFrame,
    *,
    genes: Sequence[str] | str | None = None,
    top_n: int = 10,
    mode: str = "contribution",
    ncols: int = 2,
    ax: Axes | None = None,
    figsize: tuple[float, float] | None = None,
) -> Figure:
    """Faceted stacked-bar figure, one panel per gene.

    `genes` defaults to every gene in `freq` (so whatever you asked
    `allele_frequencies` for), and also takes "classical", "all", or a list.
    Pass `ax=` to draw a single gene into an existing Axes instead.
    """
    if mode not in ("contribution", "composition"):
        raise HLAFrequencyError("mode must be 'contribution' or 'composition'")

    # Derived from the rows, not attrs, so a filtered table plots what it holds.
    wanted = select_genes(order_genes(freq["Gene"]), "all" if genes is None else genes)
    if not wanted:
        raise HLAFrequencyError("none of the requested genes are in this frequency table")
    group_by = freq.attrs.get("group_by", "superpop")
    groups = group_order(freq["Group"].unique(), group_by)
    colors = group_colors(groups)

    if ax is not None:
        if len(wanted) != 1:
            raise HLAFrequencyError("ax= draws a single gene; pass genes=['HLA-B']")
        plot_gene(freq, wanted[0], ax, top_n=top_n, mode=mode, colors=colors)
        return ax.figure

    ncols = max(1, min(ncols, len(wanted)))
    nrows = -(-len(wanted) // ncols)
    figsize = figsize or (4.2 * ncols, 1.1 + 0.34 * top_n * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    fig.patch.set_facecolor(SURFACE)

    for axis, gene in zip(axes.ravel(), wanted):
        plot_gene(freq, gene, axis, top_n=top_n, mode=mode, colors=colors)
    for axis in axes.ravel()[len(wanted):]:
        axis.set_visible(False)

    fig.legend(
        handles=[Patch(facecolor=colors[g], label=g) for g in groups],
        loc="upper left",
        bbox_to_anchor=(0.01, 1.0),
        ncol=min(len(groups), 8),
        frameon=False,
        fontsize=9,
        title="Population" if group_by == "population" else "Superpopulation",
        title_fontsize=9,
        handlelength=1.0,
        handleheight=1.0,
    )
    fig.supxlabel(
        "Allele frequency among called haplotypes"
        if mode == "contribution"
        else "Share of carriers by group",
        fontsize=9,
        color=INK_MUTED,
    )
    fig.tight_layout(rect=(0, 0.015, 1, 0.94))
    return fig


# --- the notebook-facing object ------------------------------------------


@dataclass
class AlleleFrequencies:
    """Frequencies plus the plotting helpers, in one object.

    `af.table` is a plain tidy DataFrame — filter or group it however you like
    and hand the result back to `plot_frequencies`.
    """

    table: pd.DataFrame
    group_by: str = "superpop"
    calls: pd.DataFrame | None = field(default=None, repr=False)

    @classmethod
    def from_files(
        cls,
        annot: PathOrFrame,
        metadata,
        **kwargs,
    ) -> "AlleleFrequencies":
        """Load calls + metadata and compute frequencies in one step.

        Accepts every keyword `allele_frequencies` takes (genes, resolution,
        group, exclude_novel, drop_unknown_pop, sample_col, pop_col, verbose).
        """
        calls = load_calls(annot)
        table = allele_frequencies(calls, metadata, **kwargs)
        return cls(table=table, group_by=kwargs.get("group", "superpop"), calls=calls)

    @property
    def genes(self) -> list[str]:
        return order_genes(self.table["Gene"])

    @property
    def groups(self) -> list[str]:
        return group_order(self.table["Group"].unique(), self.group_by)

    @property
    def colors(self) -> dict[str, str]:
        return group_colors(self.groups)

    def top_alleles(self, gene: str, n: int = 10, mode: str = "contribution") -> pd.DataFrame:
        return top_alleles(self.table, gene, n, mode)

    def frequency_matrix(self, mode: str = "contribution") -> pd.DataFrame:
        return frequency_matrix(self.table, mode)

    def plot(self, **kwargs) -> Figure:
        """Faceted figure; takes genes, top_n, mode, ncols, ax, figsize."""
        return plot_frequencies(self.table, **kwargs)

    def save(self, out: str | Path, *, table_out: str | Path | None = None, dpi: int = 300, **kwargs) -> Figure:
        """Plot, write the image, and write the plotted values beside it."""
        out = Path(out).expanduser()
        fig = self.plot(**kwargs)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=dpi, facecolor=SURFACE)
        shown = self.table
        genes = kwargs.get("genes")
        if genes is not None:
            shown = shown[shown["Gene"].isin(select_genes(self.genes, genes))]
        shown.to_csv(
            Path(table_out).expanduser() if table_out else out.with_suffix(".tsv"),
            sep="\t",
            index=False,
            float_format="%.5f",
        )
        return fig


# --- cli ------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    import matplotlib

    matplotlib.use("Agg")

    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--annot", required=True, type=Path, help="Immuannot allele table (TSV/CSV)")
    p.add_argument("--metadata", required=True, type=Path, help="sample metadata table")
    p.add_argument("--sample-col", default=SAMPLE_COL, help=f"sample id column (default {SAMPLE_COL})")
    p.add_argument("--pop-col", default=POP_COL, help=f"population column (default {POP_COL})")
    p.add_argument(
        "--genes",
        nargs="+",
        default=["classical"],
        help="'classical' (default), 'all' for every gene in the table "
        "(non-classical HLA and KIR included), or an explicit list",
    )
    p.add_argument("--top-n", type=int, default=10, help="alleles per gene (default 10)")
    p.add_argument(
        "--resolution", type=int, default=2, help="HLA allele fields to keep, 1-4 (default 2)"
    )
    p.add_argument(
        "--kir-resolution", type=int, default=3, help="KIR allele fields to keep, 1-3 (default 3)"
    )
    p.add_argument("--group", choices=["superpop", "population"], default="superpop")
    p.add_argument("--mode", choices=["contribution", "composition"], default="contribution")
    p.add_argument(
        "--exclude-novel",
        action="store_true",
        help="drop calls flagged 'new' (by default the flag is just stripped)",
    )
    p.add_argument("--drop-unknown-pop", action="store_true", help="drop samples with no population")
    p.add_argument("--ncols", type=int, default=2)
    p.add_argument("--out", type=Path, default=Path("hla_allele_frequency.png"))
    p.add_argument("--table-out", type=Path, help="also write the plotted values as TSV")
    args = p.parse_args(argv)

    try:
        af = AlleleFrequencies.from_files(
            args.annot,
            args.metadata,
            genes=args.genes,
            resolution={"HLA": args.resolution, "KIR": args.kir_resolution},
            group=args.group,
            exclude_novel=args.exclude_novel,
            drop_unknown_pop=args.drop_unknown_pop,
            sample_col=args.sample_col,
            pop_col=args.pop_col,
        )
        af.save(
            args.out,
            table_out=args.table_out,
            genes=args.genes,
            top_n=args.top_n,
            mode=args.mode,
            ncols=args.ncols,
        )
    except HLAFrequencyError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"wrote {args.out}", file=sys.stderr)
    print(f"wrote {args.table_out or args.out.with_suffix('.tsv')}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
