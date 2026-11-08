#!/usr/bin/env python3
"""Stacked bars of significant associations at each haplotype block.

Notebook:
    from plot_nearest_gene_hits import plot_nearest_gene_hits

    fig, df_info = plot_nearest_gene_hits(df_all)
    fig.savefig("nearest_gene_hits.pdf", bbox_inches="tight")

`df_all` needs a marker id, p_bonferroni, and pheno_category. Rows with
p_bonferroni < 0.05 are counted once. The x-axis is the block name. A
MarkerID of the form ``blockname_index`` (``DR_014a83b8``) contributes to
``DR``. A ``blockname`` or ``block`` column is used when it is present.
Each bar is the number of associations at that block, stacked by
pheno_category.
"""

from __future__ import annotations

import math
from typing import Mapping, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator


# MHC block order from hap_block_haplotypes.BLOCKS. Anything else is
# appended after these, in name order.
BLOCK_ORDER = (
    "alpha",
    "beta",
    "RCCX",
    "DR",
    "delta",
    "TAP",
    "epsilon",
    "kappa",
    "KIR",
)

# Legend order and the colours in the reference figure. Any other
# pheno_category is appended after these and given a seaborn colour.
DISEASE_ORDER = (
    "Autoimmune",
    "Cardiometabolic",
    "Infectious",
    "Neoplastic",
    "Neurologic",
    "Organ",
    "Other",
    "Rheumatic",
)
DISEASE_COLORS = {
    "Autoimmune": "#A9D4F0",
    "Cardiometabolic": "#3E73C4",
    "Infectious": "#C5E59A",
    "Neoplastic": "#F6B4B0",
    "Neurologic": "#6FBF55",
    "Organ": "#E23B3B",
    "Other": "#F6C48C",
    "Rheumatic": "#F5A12A",
}
_FALLBACK = (
    "#8DA0CB",
    "#E78AC3",
    "#A6D854",
    "#FFD92F",
    "#E5C494",
    "#66C2A5",
    "#FC8D62",
    "#B3B3B3",
)


def _block_name(marker) -> str:
    """Block from ``blockname_index``. The whole id is kept when it has no ``_``."""
    text = str(marker).strip()
    head, sep, _index = text.rpartition("_")
    return head if sep else text


def _block_series(df: pd.DataFrame, marker_col: str) -> pd.Series:
    for col in ("blockname", "block", "Block"):
        if col in df.columns:
            return df[col].map(lambda value: str(value).strip())
    if marker_col not in df.columns:
        raise ValueError("df_all is missing %s" % marker_col)
    return df[marker_col].map(_block_name)


def block_hit_table(
    df_all: pd.DataFrame,
    *,
    p_col: str = "p_bonferroni",
    alpha: float = 0.05,
    marker_col: str = "MarkerID",
    category_col: str = "pheno_category",
) -> pd.DataFrame:
    """One row per (block, category): Block, Category, Count."""
    missing = [c for c in (p_col, category_col) if c not in df_all.columns]
    if missing:
        raise ValueError("df_all is missing %s" % ", ".join(missing))

    p = pd.to_numeric(df_all[p_col], errors="coerce")
    sig = df_all.loc[p < alpha].copy()
    if sig.empty:
        raise ValueError("no rows with %s < %s" % (p_col, alpha))

    sig["_block"] = _block_series(sig, marker_col)
    sig = sig.dropna(subset=[category_col])
    sig = sig.loc[sig["_block"].ne("") & sig["_block"].ne("nan")]
    if sig.empty:
        raise ValueError("significant rows have no usable block or %s" % category_col)

    counts = (
        sig.groupby(["_block", category_col], observed=True)
        .size()
        .rename("Count")
        .reset_index()
    )
    return counts.rename(columns={"_block": "Block", category_col: "Category"})


def gene_hit_table(*args, **kwargs) -> pd.DataFrame:
    """Alias of ``block_hit_table``. Block is the grouping used on the x-axis."""
    return block_hit_table(*args, **kwargs)


def _ordered(present: Sequence[str], preferred: Sequence[str], requested: Optional[Sequence[str]]) -> list[str]:
    present_list = list(dict.fromkeys(str(x) for x in present))
    by_lower = {name.lower(): name for name in present_list}
    source = list(requested) if requested else list(preferred)
    ordered = []
    for name in source:
        hit = by_lower.get(str(name).lower())
        if hit is not None and hit not in ordered:
            ordered.append(hit)
    rest = [name for name in present_list if name not in ordered]
    if requested is None:
        rest = sorted(rest, key=str.lower)
    ordered.extend(rest)
    return ordered


def _colors(categories: Sequence[str], override: Optional[Mapping[str, str]]) -> dict[str, str]:
    named = {key.lower(): value for key, value in DISEASE_COLORS.items()}
    if override:
        named.update({str(key).lower(): value for key, value in override.items()})
    fallback = iter(_FALLBACK)
    extra = iter(sns.color_palette("colorblind", n_colors=max(len(categories), 1)).as_hex())
    out = {}
    for cat in categories:
        key = str(cat).lower()
        if key in named:
            out[cat] = named[key]
            continue
        try:
            out[cat] = next(fallback)
        except StopIteration:
            out[cat] = next(extra)
    return out


def _tick_step(ymax: float) -> float:
    if ymax <= 8:
        return 2
    if ymax <= 20:
        return 5
    if ymax <= 40:
        return 10
    if ymax <= 100:
        return 20
    if ymax <= 200:
        return 50
    return 10 ** math.floor(math.log10(max(ymax, 1))) / 2


def _font_family() -> str:
    from matplotlib import font_manager

    have = {font.name for font in font_manager.fontManager.ttflist}
    for name in ("Arial", "Helvetica", "DejaVu Sans"):
        if name in have:
            return name
    return "DejaVu Sans"


def plot_nearest_gene_hits(
    df_all: pd.DataFrame,
    *,
    p_col: str = "p_bonferroni",
    alpha: float = 0.05,
    marker_col: str = "MarkerID",
    category_col: str = "pheno_category",
    title: str = "Associations by haplotype block",
    panel: str = "",
    legend_title: str = "Disease Group",
    figsize: tuple[float, float] = (8, 4.2),
    bar_width: float = 0.72,
    category_order: Optional[Sequence[str]] = None,
    block_order: Optional[Sequence[str]] = None,
    colors: Optional[Mapping[str, str]] = None,
):
    """Stacked bar per block. Returns (figure, count table).

    Bar height is the number of rows with ``p_col < alpha``. Segments are
    ``category_col``. Blocks follow the MHC block order unless ``block_order``
    is given.
    """
    info = block_hit_table(
        df_all,
        p_col=p_col,
        alpha=alpha,
        marker_col=marker_col,
        category_col=category_col,
    )
    categories = _ordered(info["Category"].tolist(), DISEASE_ORDER, category_order)
    blocks = _ordered(info["Block"].tolist(), BLOCK_ORDER, block_order)
    palette = _colors(categories, colors)

    wide = (
        info.pivot_table(
            index="Block",
            columns="Category",
            values="Count",
            aggfunc="sum",
            fill_value=0,
            observed=True,
        )
        .reindex(blocks)
        .fillna(0)
    )
    for cat in categories:
        if cat not in wide.columns:
            wide[cat] = 0
    present = [cat for cat in categories if float(wide[cat].sum()) > 0]

    font = _font_family()
    with plt.rc_context({
        "font.family": font,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }):
        fig, ax = plt.subplots(figsize=figsize)
        fig.subplots_adjust(left=0.10, right=0.74, bottom=0.16, top=0.88)
        x = np.arange(len(blocks))
        bottom = np.zeros(len(blocks))
        for cat in present:
            height = wide[cat].to_numpy(dtype=float)
            ax.bar(
                x, height, bottom=bottom, width=bar_width,
                color=palette[cat], edgecolor="none", linewidth=0, zorder=3,
            )
            bottom += height

        ymax = float(bottom.max()) if len(bottom) else 1.0
        step = _tick_step(ymax)
        top = math.ceil(ymax / step) * step
        if top <= ymax:
            top += step
        ax.set_ylim(0, top)
        ax.yaxis.set_major_locator(MultipleLocator(step))
        ax.set_xlim(-0.6, len(blocks) - 0.4)
        ax.set_xticks(x)
        ax.set_xticklabels(blocks, fontsize=10)
        ax.set_xlabel("Block", fontsize=12)
        ax.set_ylabel("Number of associations", fontsize=12)
        ax.tick_params(axis="y", labelsize=10)
        sns.despine(ax=ax)
        if title:
            ax.set_title(title, loc="left", fontsize=13, pad=10, color="black")
        if panel:
            ax.text(
                -0.08, 1.04, panel,
                transform=ax.transAxes,
                fontsize=16, fontweight="bold",
                ha="left", va="bottom", clip_on=False,
            )
        handles = [
            Patch(facecolor=palette[cat], edgecolor="none", label=cat)
            for cat in present
        ]
        legend = ax.legend(
            handles=handles,
            title=legend_title,
            frameon=False,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            borderaxespad=0.0,
            handlelength=0.9,
            handleheight=0.9,
            labelspacing=0.45,
            fontsize=9,
        )
        legend.get_title().set_fontsize(11)

    info_out = info.copy()
    info_out["Block"] = pd.Categorical(info_out["Block"], blocks, ordered=True)
    info_out["Category"] = pd.Categorical(info_out["Category"], categories, ordered=True)
    info_out = info_out.sort_values(["Block", "Category"]).reset_index(drop=True)
    return fig, info_out


if __name__ == "__main__":
    raise SystemExit(
        "import plot_nearest_gene_hits and call it with df_all; "
        "this module does not read a file on its own"
    )
