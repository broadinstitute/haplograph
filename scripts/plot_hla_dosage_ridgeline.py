#!/usr/bin/env python3
"""Publication-ready HLA dosage × ancestry ridgeline plot.

Notebook: ensure df_meta has columns ancestry_pred, hla_dosage, then run this
cell body (or: from plot_hla_dosage_ridgeline import plot_hla_dosage_ridgeline).
"""

from __future__ import annotations

import itertools
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

OKABE_ITO = [
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#999999",
]

# Prefer this order when labels are present; others append alphabetically.
CANONICAL_ANCESTRY_ORDER = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS", "other"]


def _ancestry_order(labels: Sequence[str]) -> list[str]:
    present = list(dict.fromkeys(str(x) for x in labels))
    known = [a for a in CANONICAL_ANCESTRY_ORDER if a in present]
    rest = sorted(set(present) - set(known))
    return known + rest


def _palette_for(order: Sequence[str]) -> dict[str, str]:
    cycle = itertools.cycle(OKABE_ITO)
    return {anc: next(cycle) for anc in order}


def plot_hla_dosage_ridgeline(
    df_meta: pd.DataFrame,
    *,
    ancestry_col: str = "ancestry_pred",
    dosage_col: str = "hla_dosage",
    xlim: tuple[float, float] = (60.0, 80.0),
    out_base: str = "HLA_dosage_ridgeline",
    title: Optional[str] = "HLA dosage by ancestry",
) -> sns.FacetGrid:
    """Draw overlapping KDEs of HLA dosage by ancestry; save PDF + PNG."""
    required = {ancestry_col, dosage_col}
    missing = required - set(df_meta.columns)
    if missing:
        raise ValueError(f"df_meta missing columns: {sorted(missing)}")

    df = df_meta[[ancestry_col, dosage_col]].dropna().copy()
    order = _ancestry_order(df[ancestry_col].astype(str))
    palette = _palette_for(order)

    g = sns.FacetGrid(
        df,
        row=ancestry_col,
        hue=ancestry_col,
        row_order=order,
        hue_order=order,
        aspect=15,
        height=0.5,
        palette=palette,
    )

    g.map(
        sns.kdeplot,
        dosage_col,
        bw_adjust=0.5,
        clip_on=False,
        fill=True,
        alpha=1,
        linewidth=1.5,
    )
    g.map(
        sns.kdeplot,
        dosage_col,
        clip_on=False,
        color="w",
        lw=2,
        bw_adjust=0.5,
    )
    g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

    def label(_x, color, label):
        ax = plt.gca()
        ax.text(
            0.01,
            0.2,
            label,
            fontweight="bold",
            color=color,
            ha="left",
            va="center",
            transform=ax.transAxes,
            clip_on=False,
        )

    g.map(label, dosage_col)

    # Negative hspace needs an explicit left margin so ancestry names are not clipped.
    g.figure.subplots_adjust(left=0.18, hspace=-0.25)

    g.set_titles("")
    g.set(yticks=[], ylabel="", xlim=xlim)
    g.set_axis_labels(x_var="HLA dosage", y_var="")
    g.despine(bottom=True, left=True)

    if title:
        g.figure.suptitle(title, fontweight="bold", y=0.98)

    g.figure.savefig(f"{out_base}.pdf", bbox_inches="tight")
    g.figure.savefig(f"{out_base}.png", dpi=300, bbox_inches="tight")
    return g


def _synthetic_df(n_per: int = 400, seed: int = 0) -> pd.DataFrame:
    rng = __import__("numpy").random.default_rng(seed)
    rows = []
    centers = {
        "AFR": 66.0,
        "AMR": 65.5,
        "EAS": 64.0,
        "EUR": 65.0,
        "MID": 65.2,
        "SAS": 64.5,
        "other": 65.8,
    }
    for anc, mu in centers.items():
        for x in rng.normal(mu, 1.2, size=n_per):
            rows.append({"ancestry_pred": anc, "hla_dosage": float(x)})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    import tempfile
    from pathlib import Path

    out_dir = Path(tempfile.mkdtemp(prefix="hla_ridgeline_"))
    out_base = str(out_dir / "HLA_dosage_ridgeline")
    g = plot_hla_dosage_ridgeline(_synthetic_df(), out_base=out_base)
    pdf, png = Path(f"{out_base}.pdf"), Path(f"{out_base}.png")
    assert pdf.is_file() and png.is_file(), (pdf, png)
    # FacetGrid shares xlim across axes
    xlim = g.axes.flat[0].get_xlim()
    assert abs(xlim[0] - 60) < 1e-6 and abs(xlim[1] - 80) < 1e-6, xlim
    print(f"OK: wrote {pdf} and {png}; xlim={xlim}")
    plt.close(g.figure)
