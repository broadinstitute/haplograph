"""
plot_allele_ld — triangular linkage disequilibrium for multi-allelic genes.

Takes the DataFrame returned by `allele_pair_ld.pair_ld` and draws one cell
per allele pair. Genes run left to right in genomic order; every allele of a
gene is its own column, so a locus stays multi-allelic instead of being
collapsed to a single marker. The x-axis is that gene's GRCh38 coordinate.

The upper half is r² and the lower half is |D'|, on a shared scale that
starts at 0.25: weaker pairs are left white, which is what turns the
triangle into LD blocks.

    from allele_pair_ld import pair_ld
    from plot_allele_ld import plot_allele_ld

    pairs = pair_ld(dosage_df, population=pops)
    fig = plot_allele_ld(
        pairs, population="EUR",
        panel="C", title="Linkage Disequilibrium of HLA Region",
    )
    fig.savefig("hla_ld.png", dpi=200, bbox_inches="tight")

matplotlib only, plus the GRCh38 gene table in `allele_dosage_vcf`.

Reading the figure
------------------
Each column is one allele. Alleles of one gene sit together, genes ordered
by chromosome coordinate, and within a gene the more common allele is to
the left. Column width is one allele, not one base, so a gene with more
alleles is wider and the megabase labels are the gene coordinates rather
than a linear ruler.

A cell between columns i and j sits halfway between them, at a height of
half the number of columns that separate them. Above the midline the cell
is r² (teal); below it, |D'| (red). The sign of D' is discarded: a large
negative value is strong linkage too, and the scale runs from 0.25 to 1.
Pairs drawn below `vmin`, and pairs the table does not contain, stay white.
Within-locus pairs are absent unless `pair_ld` was run with
`include_same_locus=True`. HLA and KIR land on separate panels, because
they are not on the same chromosome; a pair that crosses panels is not drawn.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Mapping, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.figure import Figure

_SCRIPTS = str(Path(__file__).resolve().parent)
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)

from allele_dosage_vcf import GENE_COORDINATES

__all__ = ["plot_allele_ld"]


# Pale at the threshold, saturated at 1. Matches a split r² / D' triangle:
# teal above the axis, red below. Values under vmin are masked, not mapped,
# so the floor of the scale stays a tint rather than white.
_CMAP_R2 = LinearSegmentedColormap.from_list(
    "ld_r2", ["#E5F6F6", "#8FD4D4", "#2AABAB", "#0E7C7C"]
)
_CMAP_DPRIME = LinearSegmentedColormap.from_list(
    "ld_dprime", ["#F8E4E4", "#F0A0A0", "#E15B5B", "#C62828"]
)
_CMAP_R2.set_bad("white")
_CMAP_DPRIME.set_bad("white")


def plot_allele_ld(
    pairs: pd.DataFrame,
    *,
    population: Optional[str] = None,
    positions: Union[Mapping, pd.Series, None] = None,
    r2: str = "r2",
    dprime: str = "Dprime",
    vmin: float = 0.25,
    vmax: float = 1.0,
    max_alleles: Optional[int] = None,
    allele_order: str = "frequency",
    title: Optional[str] = None,
    panel: Optional[str] = None,
    figsize: Optional[tuple] = None,
) -> Figure:
    """Draw r² above the axis and |D'| below it, alleles along the genome.

    Parameters
    ----------
    pairs
        Output of `allele_pair_ld.pair_ld` (one row per allele pair).
    population
        Which `population` value to draw. With several populations and no
        argument, ALL is used when that label is present; otherwise this
        raises so a stratum is not chosen silently.
    positions
        Genomic coordinates for loci the builtin GRCh38 table does not know,
        or overrides for ones it does. A locus maps to a base-pair position,
        to ``(start, end)``, or to ``(chrom, start, end)``. A map whose keys
        are allele names (the same strings as in `pairs`) places those
        alleles directly and is ordered by that position.
    r2, dprime
        Columns used for the upper and lower halves. D' is drawn as a
        magnitude.
    vmin, vmax
        Color scale. Pairs below `vmin` are white. The default 0.25 is the
        floor of the reference figure.
    max_alleles
        Keep this many alleles per gene, the most frequent first. ``None``
        keeps every allele that survived `pair_ld`.
    allele_order
        ``"frequency"`` (common allele on the left of its gene) or
        ``"name"``.
    title, panel
        Heading. `panel` is the bold letter at the left of the title.
    figsize
        Width and height in inches. The default widens slightly with the
        number of alleles.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if vmin >= vmax:
        raise ValueError(f"vmin ({vmin}) must be below vmax ({vmax})")
    if allele_order not in {"frequency", "name"}:
        raise ValueError("allele_order must be 'frequency' or 'name'")

    data = _subset_population(pairs, population)
    _require_columns(data, ("allele_a", "allele_b", "locus_a", "locus_b", r2, dprime))

    layout = _allele_layout(
        data,
        positions=positions,
        max_alleles=max_alleles,
        allele_order=allele_order,
    )
    data = data[
        data["allele_a"].isin(layout["allele"]) & data["allele_b"].isin(layout["allele"])
    ]
    if data.empty:
        raise ValueError("no allele pairs left to draw")

    chroms = _chrom_order(layout["chrom"].unique())
    n_alleles = len(layout)
    if figsize is None:
        width = float(np.clip(7.5 + n_alleles / 45.0, 9.5, 14.5))
        height = 4.8 if len(chroms) == 1 else 3.6 * len(chroms)
        figsize = (width, height)

    fig = plt.figure(figsize=figsize, layout="constrained")
    fig.set_facecolor("white")
    gs = fig.add_gridspec(
        len(chroms) * 2,
        2,
        height_ratios=[1, 1] * len(chroms),
        width_ratios=[32, 1],
        wspace=0.025,
        hspace=0.22,
    )

    heading_used = False
    for k, chrom in enumerate(chroms):
        sub = layout[layout["chrom"] == chrom].reset_index(drop=True)
        sub = sub.copy()
        sub["x"] = np.arange(len(sub), dtype=float)
        ax = fig.add_subplot(gs[2 * k: 2 * k + 2, 0])
        cax_r = fig.add_subplot(gs[2 * k, 1])
        cax_d = fig.add_subplot(gs[2 * k + 1, 1])
        _draw_triangle(ax, data, sub, r2, dprime, vmin=vmin, vmax=vmax)
        _colorbars(fig, cax_r, cax_d, vmin=vmin, vmax=vmax)
        ax.set_xlabel(_axis_label(chrom))
        if not heading_used:
            _heading(ax, panel, title)
            heading_used = True

    return fig


# ----------------------------------------------------------------- selection


def _subset_population(pairs: pd.DataFrame, population: Optional[str]) -> pd.DataFrame:
    if not isinstance(pairs, pd.DataFrame):
        raise TypeError("pairs must be the DataFrame returned by pair_ld")
    if "population" not in pairs.columns:
        if population is not None:
            raise ValueError("pairs have no population column; drop population=")
        return pairs.copy()

    present = [str(p) for p in pd.unique(pairs["population"])]
    chosen = population
    if chosen is None:
        if len(present) == 1:
            chosen = present[0]
        elif "ALL" in present:
            chosen = "ALL"
        else:
            shown = ", ".join(present[:8])
            extra = " …" if len(present) > 8 else ""
            raise ValueError(
                f"pairs contain {len(present)} populations ({shown}{extra}); "
                "pass population="
            )
    chosen = str(chosen)
    out = pairs.loc[pairs["population"].astype(str) == chosen].copy()
    if out.empty:
        raise ValueError(
            f"population {chosen!r} is not in pairs; "
            f"present: {', '.join(present[:8])}"
        )
    return out


def _require_columns(pairs: pd.DataFrame, columns: Sequence[str]) -> None:
    missing = [c for c in columns if c not in pairs.columns]
    if missing:
        raise ValueError(
            f"pairs is missing {missing}; expected the columns pair_ld writes "
            f"({', '.join(columns)})"
        )


def _allele_layout(
    pairs: pd.DataFrame,
    *,
    positions: Union[Mapping, pd.Series, None],
    max_alleles: Optional[int],
    allele_order: str,
) -> pd.DataFrame:
    """One row per allele, in the order the triangle's columns are drawn."""
    locus_of = _allele_loci(pairs)
    freq = _allele_frequency(pairs) if allele_order == "frequency" else {}

    rows = []
    for allele, locus in locus_of.items():
        rows.append(
            {
                "allele": allele,
                "locus": locus,
                "freq": float(freq.get(allele, 0.0)),
            }
        )
    layout = pd.DataFrame(rows)
    if layout.empty:
        raise ValueError("no alleles in pairs")

    if max_alleles is not None:
        if max_alleles < 1:
            raise ValueError("max_alleles must be at least 1")
        keep = []
        for _, sub in layout.groupby("locus", sort=False):
            sub = sub.sort_values(["freq", "allele"], ascending=[False, True])
            keep.append(sub.head(int(max_alleles)))
        layout = pd.concat(keep, ignore_index=True)

    coords = _resolve_coordinates(layout, positions)
    located = layout["locus"].map(coords)
    missing = sorted(layout.loc[located.isna(), "locus"].unique())
    if missing:
        shown = ", ".join(missing[:8])
        extra = " …" if len(missing) > 8 else ""
        raise ValueError(
            f"no genomic coordinate for {shown}{extra}. "
            "Pass positions={'GENE': bp} or positions={'GENE': (chrom, start, end)}."
        )
    layout["chrom"] = located.map(lambda c: c[0])
    layout["start"] = located.map(lambda c: c[1])
    layout["end"] = located.map(lambda c: c[2])
    layout["pos"] = (layout["start"] + layout["end"]) / 2.0

    # Per-allele coordinates replace the gene interval for those alleles.
    allele_pos = _allele_position_map(positions, set(layout["allele"]))
    if allele_pos:
        hit = layout["allele"].map(allele_pos)
        layout.loc[hit.notna(), "pos"] = hit.dropna()
        layout.loc[hit.notna(), "start"] = hit.dropna()
        layout.loc[hit.notna(), "end"] = hit.dropna()

    layout["_chrom"] = layout["chrom"].map(_chrom_key)
    sort_cols = ["_chrom", "pos", "locus"]
    ascending = [True, True, True]
    if allele_order == "frequency":
        sort_cols.append("freq")
        ascending.append(False)
    sort_cols.append("allele")
    ascending.append(True)
    layout = layout.sort_values(sort_cols, ascending=ascending, kind="mergesort")
    layout = layout.drop(columns="_chrom").reset_index(drop=True)
    return layout


def _allele_loci(pairs: pd.DataFrame) -> dict:
    """Allele -> locus, refusing an allele that the table assigns two loci."""
    out: dict = {}
    for allele_col, locus_col in (("allele_a", "locus_a"), ("allele_b", "locus_b")):
        for allele, locus in zip(pairs[allele_col].astype(str), pairs[locus_col].astype(str)):
            prev = out.get(allele)
            if prev is not None and prev != locus:
                raise ValueError(f"allele {allele} is listed at both {prev} and {locus}")
            out[allele] = locus
    return out


def _allele_frequency(pairs: pd.DataFrame) -> dict:
    if not {"allele_freq_a", "allele_freq_b"} <= set(pairs.columns):
        return {}
    parts = [
        pairs[["allele_a", "allele_freq_a"]].rename(
            columns={"allele_a": "allele", "allele_freq_a": "freq"}
        ),
        pairs[["allele_b", "allele_freq_b"]].rename(
            columns={"allele_b": "allele", "allele_freq_b": "freq"}
        ),
    ]
    stacked = pd.concat(parts, ignore_index=True)
    stacked["allele"] = stacked["allele"].astype(str)
    return stacked.groupby("allele")["freq"].max().to_dict()


# --------------------------------------------------------------- coordinates


def _bare_locus(name: str) -> str:
    """'HLA-DRB1' / 'HLA_A' -> 'DRB1' / 'A'. KIR names are left intact."""
    s = str(name).strip()
    upper = s.upper()
    for prefix in ("HLA-", "HLA_", "HLA."):
        if upper.startswith(prefix):
            return s[len(prefix):]
    return s


def _builtin_coordinates() -> dict:
    """Bare locus -> (chrom, start, end)."""
    out = {}
    for name, (chrom, start, end) in GENE_COORDINATES.items():
        out[_bare_locus(name)] = (str(chrom), float(start), float(end))
    return out


def _parse_coordinate(value, default_chrom: str):
    """A position, an interval, or (chrom, start, end)."""
    if isinstance(value, str):
        raise TypeError(f"coordinate {value!r} is text; expected a position or interval")
    if isinstance(value, (int, float, np.floating)):
        v = float(value)
        return default_chrom, v, v
    seq = list(value)
    if len(seq) == 2 and not isinstance(seq[0], str):
        return default_chrom, float(seq[0]), float(seq[1])
    if len(seq) >= 3 and isinstance(seq[0], str):
        return str(seq[0]), float(seq[1]), float(seq[2])
    raise TypeError(
        "a coordinate must be bp, (start, end), or (chrom, start, end); "
        f"got {value!r}"
    )


def _resolve_coordinates(layout: pd.DataFrame, positions) -> dict:
    """Locus -> (chrom, start, end), builtin table updated by `positions`."""
    coords = _builtin_coordinates()
    if positions is None:
        return {locus: coords.get(_bare_locus(locus)) for locus in layout["locus"].unique()}

    mapping = positions.to_dict() if isinstance(positions, pd.Series) else dict(positions)
    alleles = set(layout["allele"].astype(str))
    # An allele map is handled separately; don't also treat those keys as loci.
    if mapping and set(map(str, mapping)).issubset(alleles):
        return {locus: coords.get(_bare_locus(locus)) for locus in layout["locus"].unique()}

    for key, value in mapping.items():
        bare = _bare_locus(key)
        default = coords.get(bare, ("", np.nan, np.nan))[0] or ""
        coords[bare] = _parse_coordinate(value, default or "")

    return {locus: coords.get(_bare_locus(locus)) for locus in layout["locus"].unique()}


def _allele_position_map(positions, alleles: set) -> dict:
    if positions is None:
        return {}
    mapping = positions.to_dict() if isinstance(positions, pd.Series) else dict(positions)
    if not mapping or not set(map(str, mapping)).issubset(alleles):
        return {}
    out = {}
    for key, value in mapping.items():
        chrom, start, end = _parse_coordinate(value, "")
        out[str(key)] = (float(start) + float(end)) / 2.0
        del chrom
    return out


def _chrom_key(chrom: str):
    body = str(chrom).lower().removeprefix("chr")
    if body.isdigit():
        return (0, int(body))
    return (1, str(chrom))


def _chrom_order(chroms) -> list:
    return sorted((str(c) for c in chroms), key=_chrom_key)


def _axis_label(chrom: str) -> str:
    if chrom:
        return f"Genomic position ({chrom}, Mb)"
    return "Genomic position (Mb)"


# --------------------------------------------------------------------- diamonds


def _indexed_pairs(pairs: pd.DataFrame, layout: pd.DataFrame, r2: str, dprime: str):
    """Allele columns i < j, with r² and |D'| averaged if a pair is repeated."""
    index = {allele: i for i, allele in enumerate(layout["allele"].astype(str))}
    a = pairs["allele_a"].astype(str).to_numpy()
    b = pairs["allele_b"].astype(str).to_numpy()
    ia = np.fromiter((index.get(x, -1) for x in a), dtype=int, count=len(a))
    ib = np.fromiter((index.get(x, -1) for x in b), dtype=int, count=len(b))
    rv = pairs[r2].to_numpy(dtype=float)
    dv = np.abs(pairs[dprime].to_numpy(dtype=float))

    ok = (ia >= 0) & (ib >= 0) & (ia != ib)
    ia, ib, rv, dv = ia[ok], ib[ok], rv[ok], dv[ok]
    swap = ia > ib
    ia, ib = np.where(swap, ib, ia), np.where(swap, ia, ib)

    acc = pd.DataFrame({"i": ia, "j": ib, "r": rv, "d": dv})
    acc = acc.groupby(["i", "j"], as_index=False).mean()
    if acc.empty:
        empty = np.array([], dtype=int)
        return empty, empty, np.array([]), np.array([])
    return (
        acc["i"].to_numpy(dtype=int),
        acc["j"].to_numpy(dtype=int),
        acc["r"].to_numpy(dtype=float),
        acc["d"].to_numpy(dtype=float),
    )


def _diamonds(ii: np.ndarray, jj: np.ndarray, sign: float) -> np.ndarray:
    """Unit diamonds centered halfway between columns i and j.

    Adjacent diamonds share an edge, so the triangle is solid rather than a
    checkerboard of pixels. `sign` is +1 above the axis and -1 below it.
    """
    x = (ii + jj) / 2.0
    y = sign * (jj - ii) / 2.0
    verts = np.empty((len(ii), 4, 2), dtype=float)
    verts[:, 0, 0] = x
    verts[:, 0, 1] = y - 0.5
    verts[:, 1, 0] = x + 0.5
    verts[:, 1, 1] = y
    verts[:, 2, 0] = x
    verts[:, 2, 1] = y + 0.5
    verts[:, 3, 0] = x - 0.5
    verts[:, 3, 1] = y
    return verts


def _add_diamonds(ax, ii, jj, values, *, sign: float, cmap, vmin: float, vmax: float) -> None:
    keep = np.isfinite(values) & (values >= vmin)
    if not np.any(keep):
        return
    verts = _diamonds(ii[keep], jj[keep], sign)
    colors = cmap(Normalize(vmin=vmin, vmax=vmax, clip=True)(values[keep]))
    ax.add_collection(
        PolyCollection(
            verts,
            facecolors=colors,
            edgecolors="none",
            antialiased=False,
            rasterized=True,
            zorder=1,
        )
    )


def _draw_triangle(ax, pairs, layout, r2, dprime, *, vmin: float, vmax: float) -> None:
    n = len(layout)
    ii, jj, rv, dv = _indexed_pairs(pairs, layout, r2, dprime)
    _add_diamonds(ax, ii, jj, rv, sign=1.0, cmap=_CMAP_R2, vmin=vmin, vmax=vmax)
    _add_diamonds(ax, ii, jj, dv, sign=-1.0, cmap=_CMAP_DPRIME, vmin=vmin, vmax=vmax)

    half = (n - 1) / 2.0 + 0.55
    ax.axhline(0.0, color="#1a1a1a", lw=0.8, zorder=3)
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_ylim(-half, half)
    ax.set_facecolor("white")
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_color("#222222")
    ax.spines["bottom"].set_linewidth(0.6)
    ax.tick_params(axis="y", left=False, labelleft=False)
    _gene_ticks(ax, layout)


def _gene_ticks(ax, layout: pd.DataFrame) -> None:
    """One tick per gene, at the middle of its alleles, labelled with Mb."""
    centers = []
    labels = []
    for locus, sub in layout.groupby("locus", sort=False):
        centers.append(float(sub["x"].mean()))
        pos = float(sub["pos"].mean())
        mb = pos / 1e6 if pos > 1000 else pos
        labels.append(f"{locus}\n{mb:.2f}")

    rotation = 0
    ha = "center"
    fontsize = 8
    if len(labels) > 10:
        labels = [lab.replace("\n", " ") for lab in labels]
        rotation = 60
        ha = "right"
        fontsize = 7
    ax.set_xticks(centers)
    ax.set_xticklabels(labels, fontsize=fontsize, rotation=rotation, ha=ha, color="#1a1a1a")
    ax.tick_params(axis="x", length=3, width=0.6, color="#222222", pad=2)


def _colorbars(fig, cax_r, cax_d, *, vmin: float, vmax: float) -> None:
    rmap = plt.cm.ScalarMappable(cmap=_CMAP_R2)
    dmap = plt.cm.ScalarMappable(cmap=_CMAP_DPRIME)
    rmap.set_clim(vmin, vmax)
    dmap.set_clim(vmin, vmax)
    ticks = np.linspace(vmin, vmax, 4)

    for cax, sm, label in (
        (cax_r, rmap, r"$r^2$"),
        (cax_d, dmap, r"$D'$"),
    ):
        bar = fig.colorbar(sm, cax=cax)
        bar.set_ticks(ticks)
        bar.set_ticklabels([f"{t:.2f}" for t in ticks])
        bar.ax.tick_params(labelsize=8, length=2, width=0.4, color="#444444")
        bar.outline.set_linewidth(0.4)
        bar.outline.set_edgecolor("#444444")
        cax.set_title(label, fontsize=11, pad=4, color="#1a1a1a")


def _heading(ax, panel: Optional[str], title: Optional[str]) -> None:
    if not panel and not title:
        return
    # A blank title reserves the margin; the two pieces are placed in points
    # so the letter cannot sit on top of the first word.
    ax.set_title(" ", loc="left", pad=2, fontsize=6)
    gap = (14 + 9 * len(str(panel))) if panel else 0
    if panel:
        ax.annotate(
            str(panel),
            xy=(0, 1),
            xycoords="axes fraction",
            xytext=(0, 5),
            textcoords="offset points",
            ha="left",
            va="bottom",
            fontsize=15,
            fontweight="bold",
            color="#111111",
            annotation_clip=False,
        )
    if title:
        ax.annotate(
            title,
            xy=(0, 1),
            xycoords="axes fraction",
            xytext=(gap, 6),
            textcoords="offset points",
            ha="left",
            va="bottom",
            fontsize=13,
            color="#161616",
            annotation_clip=False,
        )
