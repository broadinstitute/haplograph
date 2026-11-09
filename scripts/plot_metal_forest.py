"""Forest plot of a saige_metal.py table.

Notebook:
    from plot_metal_forest import plot_metal_forest

    fig, points = plot_metal_forest(df, top_n=8)
    fig.savefig("forest.pdf", bbox_inches="tight")

    # panel titles are df["phenotype"], not the phecode
    fig, points = plot_metal_forest(df, top_n=6, effect="beta")

One panel per phecode, showing the top alleles by meta-analysis p-value.
Within an allele, META and each ancestry are dodged. Color is the source
of the estimate. Point area grows with -log10(p). `effect="or"` (default)
plots exp(beta) for SAIGE Allele2, with a Wald interval. `effect="beta"`
plots the coefficient. Ancestry columns are whatever `{LABEL}_BETA`,
`{LABEL}_SE`, and `{LABEL}_P` the table has.

The y label is MarkerName. When the table has MarkerID, that string — the
HI haplotype — is wrapped underneath the allele name on the same row.
An estimate with p at or above `alpha` (default 1 − ci, so 0.05 for a 95%
interval) is drawn as an open circle.

`points` is the long table of estimates that were drawn.

One allele across phenotypes, from stacked SAIGE mega-analysis output:

    from plot_metal_forest import plot_allele_phewas_forest

    fig, rows = plot_allele_phewas_forest(df_all, "DR_014a83b8", top_n=15,
                                          title="HLA-DRB1*15:01",
                                          category_col="pheno_category")
"""

from __future__ import annotations

import math
import sys
import textwrap
from statistics import NormalDist
from typing import Mapping, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter


MISSING = {"", "NA", "NaN", "nan", "NAN", ".", "None", "null", "NULL"}
P_CANDIDATES = ("P-value", "P", "p.value", "PVALUE", "p_value", "pval")
# META, then the order used by the ancestry legend in a typical HLA forest.
PREFERRED_STUDIES = ("EUR", "AFR", "AMR", "EAS", "SAS", "MID")
# Ancestry colors match hla_allele_frequency.GROUP_COLORS. META is the
# pooled estimate, which that figure does not draw, so it stays near-black.
GROUP_COLORS = {
    "META": "#222222",
    "AFR": "#FF7F00",
    "AMR": "#ECB01F",
    "EAS": "#778500",
    "EUR": "#377EB8",
    "SAS": "#C44CFD",
    "MID": "#4A90E2",
    "OCE": "#E41A1C",
}
OTHER_COLOR = "#898781"
SERIES_COLORS = (
    "#2a78d6",
    "#eb6834",
    "#1baf7a",
    "#eda100",
    "#e87ba4",
    "#6250d6",
    "#008300",
    "#e34948",
)
# Marker area (pt^2) grows linearly with -log10(p) and saturates at 1e-20
# so a genome-wide hit stays larger than a nominal one without swallowing the row.
P_FLOOR = 1e-300
NLP_SATURATE = 20.0
AREA_MIN = 22.0
AREA_MAX = 210.0
REQUIRED = ("phecode", "MarkerName", "META_BETA", "META_SE")
# MarkerID is the HI haplotype. Wrap it under the allele name so a long
# gene string stays in the margin instead of running off the row.
HAP_WRAP = 36
_HAP_BREAKS = set(",_~|+/;")
PHENO_WRAP = 40


def z_crit(level: float) -> float:
    """Normal quantile for a central Wald interval of this coverage."""
    if not 0.0 < level < 1.0:
        raise ValueError(f"ci must be between 0 and 1, got {level}")
    return NormalDist().inv_cdf((1.0 + level) / 2.0)


def marker_area(p) -> np.ndarray:
    """Point area from a p-value. Larger means smaller p."""
    values = np.clip(np.asarray(p, dtype=float), P_FLOOR, 1.0)
    frac = np.clip(-np.log10(values), 0.0, NLP_SATURATE) / NLP_SATURATE
    return AREA_MIN + (AREA_MAX - AREA_MIN) * frac


def clean_marker_id(value) -> str:
    """HI haplotype text from MarkerID. Empty when the cell is missing."""
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    text = str(value).strip()
    if not text or text in MISSING:
        return ""
    return text


def wrap_marker_id(text: str, width: int = HAP_WRAP) -> str:
    """Break a haplotype after a separator so each line stays near `width`."""
    if len(text) <= width:
        return text
    lines = []
    start = 0
    while start < len(text):
        if len(text) - start <= width:
            lines.append(text[start:])
            break
        end = start + width
        window = text[start:end]
        br = -1
        for i, ch in enumerate(window):
            if ch in _HAP_BREAKS:
                br = i
        if br > 0:
            end = start + br + 1
        lines.append(text[start:end])
        start = end
    return "\n".join(lines)


def wrap_label(text: str, width: int = PHENO_WRAP) -> str:
    """Break a phenotype description on spaces so each line stays near `width`."""
    return textwrap.fill(text, width=width, break_long_words=False, break_on_hyphens=False)


def row_label(name: str, marker_id) -> str:
    """Allele name, with the HI haplotype wrapped on the lines below it."""
    hap = wrap_marker_id(clean_marker_id(marker_id))
    if not hap or hap.replace("\n", "") == name:
        return name
    return name + "\n" + hap


def _label_metrics(df: pd.DataFrame) -> tuple[int, int]:
    """Tallest y label, in lines, and the longest line, in characters."""
    if df.empty:
        return 1, 0
    if "MarkerID" not in df.columns:
        longest = int(df["MarkerName"].astype(str).str.len().max())
        return 1, longest
    n_lines = 1
    longest = 0
    for name, marker_id in zip(df["MarkerName"], df["MarkerID"]):
        parts = row_label(str(name), marker_id).split("\n")
        n_lines = max(n_lines, len(parts))
        longest = max(longest, max(len(part) for part in parts))
    return n_lines, longest


def _panel_span(n_alleles: int, n_lines: int) -> tuple[float, float]:
    """Y pad and axes height (inches).

    The allele name stays on the row. The wrapped HI haplotype hangs under
    that name, so the gap to the next allele has to clear it, and the bottom
    pad has to clear the last one.
    """
    if n_lines <= 1:
        return 0.65, 0.48 * n_alleles
    n_hap = n_lines - 1
    name_pt = 8.0
    hap_pt = n_hap * 6.5 * 1.15
    # Each allele name, the gap under it, the haplotype, and air before the next.
    unit_pt = name_pt + 8.0 + hap_pt + 8.0
    pad_bottom = (8.0 + hap_pt + 6.0) / unit_pt
    pad_top = (name_pt * 0.6 + 10.0) / unit_pt
    pad = max(pad_bottom, pad_top, 0.7)
    span = (n_alleles - 1) + 2 * pad
    return pad, (unit_pt / 72.0) * span


def format_p(p: float) -> str:
    if p >= 0.001:
        return f"{p:g}"
    return f"1e-{int(round(-math.log10(p)))}"


def parse_numeric(series: pd.Series) -> pd.Series:
    text = series.astype(str).str.strip()
    text = text.mask(text.isin(MISSING))
    text = text.mask(text.str.startswith("<", na=False), text.str.slice(1))
    return pd.to_numeric(text, errors="coerce")


def resolve_p_column(columns: Sequence[str]) -> str:
    have = set(columns)
    for name in P_CANDIDATES:
        if name in have:
            return name
    raise ValueError(
        "no p-value column (expected one of %s)" % ", ".join(P_CANDIDATES)
    )


def available_studies(columns: Sequence[str]) -> list[str]:
    """Ancestry labels that have BETA, SE, and P columns."""
    have = set(columns)
    labels = []
    for col in columns:
        if not col.endswith("_BETA"):
            continue
        label = col[: -len("_BETA")]
        if not label or label == "META":
            continue
        if f"{label}_SE" in have and f"{label}_P" in have:
            labels.append(label)
    return labels


def order_studies(labels: Sequence[str], requested: Optional[Sequence[str]]) -> list[str]:
    present = list(dict.fromkeys(labels))
    if requested:
        missing = [label for label in requested if label not in present]
        if missing:
            raise ValueError(
                "requested studies not in the table: %s (have %s)"
                % (", ".join(missing), ", ".join(present) or "none")
            )
        return list(requested)
    known = [label for label in PREFERRED_STUDIES if label in present]
    rest = sorted(set(present) - set(known))
    return known + rest


def group_colors(groups: Sequence[str]) -> dict[str, str]:
    colors = {}
    unknown = 0
    for group in groups:
        if group in GROUP_COLORS:
            colors[group] = GROUP_COLORS[group]
        elif group in ("OTH", "Other", "Unknown"):
            colors[group] = OTHER_COLOR
        else:
            colors[group] = SERIES_COLORS[unknown % len(SERIES_COLORS)]
            unknown += 1
    return colors


def dodge_offsets(groups: Sequence[str], span: float = 0.56) -> dict[str, float]:
    """Vertical offsets around an allele row. First group sits highest."""
    n = len(groups)
    if n <= 1:
        return {groups[0]: 0.0} if n == 1 else {}
    step = span / (n - 1)
    return {group: span / 2.0 - i * step for i, group in enumerate(groups)}


def _attach_numbers(df: pd.DataFrame, studies: Sequence[str]) -> pd.DataFrame:
    work = df.copy()
    work["phecode"] = work["phecode"].astype(str).str.strip()
    work["MarkerName"] = work["MarkerName"].astype(str).str.strip()
    p_col = resolve_p_column(work.columns)
    work["META_BETA"] = parse_numeric(work["META_BETA"])
    work["META_SE"] = parse_numeric(work["META_SE"])
    work["_meta_p"] = parse_numeric(work[p_col])
    if "META_AF" in work.columns:
        work["META_AF"] = parse_numeric(work["META_AF"])
    if "N" in work.columns:
        work["N"] = parse_numeric(work["N"])
    if "HetISq" in work.columns:
        work["HetISq"] = parse_numeric(work["HetISq"])
    for label in studies:
        for suffix in ("BETA", "SE", "P", "AF", "N"):
            col = f"{label}_{suffix}"
            if col in work.columns:
                work[col] = parse_numeric(work[col])
    return work


def select_top(
    df: pd.DataFrame,
    studies: Sequence[str],
    phecode_order: Sequence[str],
    top_n: int,
    min_af: Optional[float],
    min_n: Optional[float],
) -> pd.DataFrame:
    """Wide rows: the top_n markers per phecode, ranked by meta p-value."""
    if top_n < 1:
        raise ValueError(f"top_n must be >= 1, got {top_n}")
    work = _attach_numbers(df, studies)
    keep = np.isfinite(work["META_BETA"]) & np.isfinite(work["META_SE"]) & np.isfinite(work["_meta_p"])
    keep &= work["META_SE"] > 0
    keep &= work["_meta_p"].between(0.0, 1.0)
    keep &= work["MarkerName"] != ""
    if min_af is not None:
        if "META_AF" not in work.columns:
            raise ValueError("min_af was set but META_AF is not in the table")
        keep &= work["META_AF"] >= min_af
    if min_n is not None:
        if "N" not in work.columns:
            raise ValueError("min_n was set but N is not in the table")
        keep &= work["N"] >= min_n
    work = work.loc[keep].copy()
    if work.empty:
        raise ValueError("no rows with a finite META_BETA, META_SE, and p-value")
    work["_abs_beta"] = work["META_BETA"].abs()
    # Lowest p wins a duplicate marker. Ties break toward a larger |beta|, then name.
    work = work.sort_values(
        ["_meta_p", "_abs_beta", "MarkerName"],
        ascending=[True, False, True],
        kind="mergesort",
    )
    work = work.drop_duplicates(["phecode", "MarkerName"], keep="first")
    parts = []
    for phecode in phecode_order:
        sub = work.loc[work["phecode"] == phecode]
        if sub.empty:
            print(f"warning: {phecode} has no plottable alleles", file=sys.stderr)
            continue
        sub = sub.head(top_n).copy()
        sub["rank"] = np.arange(1, len(sub) + 1)
        parts.append(sub)
    if not parts:
        raise ValueError("no alleles to plot")
    return pd.concat(parts, ignore_index=True)


def _estimate(beta: float, se: float, effect: str, z: float) -> tuple[float, float, float]:
    low = beta - z * se
    high = beta + z * se
    if effect == "beta":
        return beta, low, high
    return math.exp(beta), math.exp(low), math.exp(high)


def forest_frame(
    selected: pd.DataFrame,
    studies: Sequence[str],
    effect: str,
    ci: float,
) -> pd.DataFrame:
    """One row per drawn point: estimate and Wald interval on the plot scale."""
    if effect not in ("or", "beta"):
        raise ValueError(f"effect must be 'or' or 'beta', got {effect}")
    z = z_crit(ci)
    groups = ["META", *studies]
    n_by_phecode = selected.groupby("phecode", sort=False)["MarkerName"].nunique().to_dict()
    records = []
    for data in selected.to_dict("records"):
        # rank 1 sits on the top row of its panel.
        y = int(n_by_phecode[data["phecode"]]) - int(data["rank"])
        for group in groups:
            if group == "META":
                beta, se, p = data["META_BETA"], data["META_SE"], data["_meta_p"]
                af = data.get("META_AF", np.nan)
                n = data.get("N", np.nan)
            else:
                beta = data.get(f"{group}_BETA", np.nan)
                se = data.get(f"{group}_SE", np.nan)
                p = data.get(f"{group}_P", np.nan)
                af = data.get(f"{group}_AF", np.nan)
                n = data.get(f"{group}_N", np.nan)
            if not (np.isfinite(beta) and np.isfinite(se) and se > 0 and np.isfinite(p)):
                continue
            if p < 0 or p > 1:
                continue
            estimate, ci_low, ci_high = _estimate(float(beta), float(se), effect, z)
            records.append(
                {
                    "phecode": data["phecode"],
                    "MarkerName": data["MarkerName"],
                    "MarkerID": clean_marker_id(data.get("MarkerID", "")),
                    "SAIGE_Allele1": data.get("SAIGE_Allele1", ""),
                    "SAIGE_Allele2": data.get("SAIGE_Allele2", ""),
                    "rank": int(data["rank"]),
                    "y": y,
                    "group": group,
                    "beta": float(beta),
                    "se": float(se),
                    "p": float(p),
                    "af": af if np.isfinite(af) else np.nan,
                    "n": n if np.isfinite(n) else np.nan,
                    "estimate": estimate,
                    "ci_low": ci_low,
                    "ci_high": ci_high,
                    "HetISq": data.get("HetISq", np.nan),
                }
            )
    if not records:
        raise ValueError("no finite estimates to draw")
    return pd.DataFrame.from_records(records)


def _data_xlim(plotted: pd.DataFrame, effect: str) -> tuple[float, float]:
    """Shared axis window. Every point is inside; extreme CI tails may clip."""
    est = plotted["estimate"].to_numpy(dtype=float)
    ends = np.concatenate(
        [plotted["ci_low"].to_numpy(dtype=float), plotted["ci_high"].to_numpy(dtype=float)]
    )
    null = 1.0 if effect == "or" else 0.0
    if effect == "or":
        est = est[np.isfinite(est) & (est > 0)]
        ends = ends[np.isfinite(ends) & (ends > 0)]
        if est.size == 0:
            return 0.25, 4.0
        lo = min(float(np.percentile(ends, 5)) if ends.size else est.min(), float(est.min()), null)
        hi = max(float(np.percentile(ends, 95)) if ends.size else est.max(), float(est.max()), null)
        lo = max(lo, 1.0 / 64.0)
        hi = min(max(hi, lo * 1.5), 64.0)
        hi = max(hi, float(est.max()), null)
        log_lo, log_hi = math.log(lo), math.log(hi)
        span = max(log_hi - log_lo, math.log(1.5))
        pad = 0.07 * span
        return math.exp(log_lo - pad), math.exp(log_hi + pad)
    finite_est = est[np.isfinite(est)]
    finite_ends = ends[np.isfinite(ends)]
    if finite_est.size == 0:
        return -1.0, 1.0
    lo = min(float(np.percentile(finite_ends, 5)) if finite_ends.size else finite_est.min(), float(finite_est.min()), null)
    hi = max(float(np.percentile(finite_ends, 95)) if finite_ends.size else finite_est.max(), float(finite_est.max()), null)
    span = max(hi - lo, 0.5)
    pad = 0.07 * span
    return lo - pad, hi + pad


def _full_xlim(plotted: pd.DataFrame, effect: str) -> tuple[float, float]:
    """Axis window that holds every interval and the null line."""
    ends = np.concatenate(
        [plotted["ci_low"].to_numpy(dtype=float), plotted["ci_high"].to_numpy(dtype=float)]
    )
    ends = ends[np.isfinite(ends)]
    if effect == "or":
        ends = ends[ends > 0]
        if ends.size == 0:
            return 0.25, 4.0
        log_lo = math.log(max(min(float(ends.min()), 1.0), 1.0 / 64.0))
        log_hi = math.log(min(max(float(ends.max()), 1.0), 64.0))
        pad = 0.06 * max(log_hi - log_lo, math.log(1.5))
        return math.exp(log_lo - pad), math.exp(log_hi + pad)
    if ends.size == 0:
        return -1.0, 1.0
    lo, hi = min(float(ends.min()), 0.0), max(float(ends.max()), 0.0)
    pad = 0.06 * max(hi - lo, 0.5)
    return lo - pad, hi + pad


def _format_or(value: float, _pos: int) -> str:
    if value <= 0:
        return ""
    if value >= 1 and abs(value - round(value)) < 1e-6:
        return str(int(round(value)))
    return f"{value:.3g}"


def _style_x(ax, effect: str, xlim: tuple[float, float]) -> None:
    null = 1.0 if effect == "or" else 0.0
    if effect == "or":
        ax.set_xscale("log", base=2)
        decades = [tick for tick in (0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64) if xlim[0] <= tick <= xlim[1]]
        if len(decades) < 2:
            fine = (0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2)
            decades = [tick for tick in fine if xlim[0] <= tick <= xlim[1]]
        if len(decades) >= 2:
            ax.xaxis.set_major_locator(FixedLocator(decades))
        ax.xaxis.set_major_formatter(FuncFormatter(_format_or))
        ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(xlim)
    ax.axvline(null, color="#666666", linestyle=(0, (4, 2.5)), linewidth=0.8, zorder=1)
    ax.set_xlabel("Odds ratio" if effect == "or" else "β")


def _scatter_group(ax, x, y, p, color, alpha: Optional[float], filled=None) -> None:
    """Filled circles are significant. Open circles are at or above `alpha`.

    `filled` overrides the p < alpha test, e.g. with a Bonferroni-adjusted call.
    """
    area = marker_area(p)
    if alpha is None and filled is None:
        ax.scatter(x, y, s=area, c=color, linewidths=0, zorder=3)
        return
    significant = np.asarray(filled, dtype=bool) if filled is not None else p < alpha
    if np.any(significant):
        ax.scatter(x[significant], y[significant], s=area[significant], c=color, linewidths=0, zorder=3)
    if np.any(~significant):
        # White fill keeps the interval from showing through the ring.
        ax.scatter(
            x[~significant],
            y[~significant],
            s=area[~significant],
            facecolors="white",
            edgecolors=color,
            linewidths=1.05,
            zorder=4,
        )


def _draw_panel(ax, panel: pd.DataFrame, offsets, colors, effect: str, xlim, alpha: Optional[float]) -> None:
    allele_order = panel[["MarkerName", "rank", "y"]].drop_duplicates("MarkerName").sort_values("rank")
    y_of = dict(zip(allele_order["MarkerName"], allele_order["y"]))
    # Guides sit on the allele, behind the dodged points.
    for y in y_of.values():
        ax.axhline(y, color="#d5d5d5", linestyle=(0, (1, 1.4)), linewidth=0.8, zorder=0)
    # Ancestry first, META last, so the pooled estimate is drawn on top.
    draw_groups = [group for group in offsets if group != "META"]
    if "META" in offsets:
        draw_groups.append("META")
    for group in draw_groups:
        offset = offsets[group]
        sub = panel.loc[panel["group"] == group]
        if sub.empty:
            continue
        x = sub["estimate"].to_numpy(dtype=float)
        lo = sub["ci_low"].to_numpy(dtype=float)
        hi = sub["ci_high"].to_numpy(dtype=float)
        y = sub["y"].to_numpy(dtype=float) + offset
        p = sub["p"].to_numpy(dtype=float)
        if effect == "or":
            positive = (x > 0) & (lo > 0) & (hi > 0) & np.isfinite(x)
            x, lo, hi, y, p = x[positive], lo[positive], hi[positive], y[positive], p[positive]
        if x.size == 0:
            continue
        xerr = np.vstack([np.maximum(x - lo, 0.0), np.maximum(hi - x, 0.0)])
        color = colors[group]
        ax.errorbar(
            x,
            y,
            xerr=xerr,
            fmt="none",
            ecolor=color,
            elinewidth=0.9,
            capsize=2.0,
            capthick=0.8,
            zorder=2,
        )
        _scatter_group(ax, x, y, p, color, alpha)
    markers = allele_order["MarkerName"].tolist()
    id_of = {}
    if "MarkerID" in panel.columns:
        for name, marker_id in zip(panel["MarkerName"], panel["MarkerID"]):
            id_of.setdefault(name, marker_id)
    labels = [row_label(name, id_of.get(name, "")) for name in markers]
    ax.set_yticks([y_of[marker] for marker in markers])
    ax.set_yticklabels(labels)
    has_hap = any("\n" in label for label in labels)
    for tick in ax.get_yticklabels():
        tick.set_fontsize(7 if has_hap else 8)
        tick.set_linespacing(1.15)
        # The tick text only reserves the margin. The visible name and the
        # haplotype under it are drawn separately so the name stays on the row.
        if has_hap:
            tick.set_alpha(0)
        else:
            tick.set_fontstyle("italic")
    n_lines = max(label.count("\n") + 1 for label in labels)
    pad, _height = _panel_span(len(markers), n_lines)
    ax.set_ylim(min(y_of.values()) - pad, max(y_of.values()) + pad)
    if has_hap:
        ax._haplograph_size_labels = True
        for marker, label in zip(markers, labels):
            name, _, hap = label.partition("\n")
            y = y_of[marker]
            ax.annotate(
                name,
                xy=(0, y),
                xycoords=("axes fraction", "data"),
                xytext=(-3, 0),
                textcoords="offset points",
                ha="right",
                va="center",
                fontsize=8,
                fontstyle="italic",
                color="#222222",
                annotation_clip=False,
            )
            if hap:
                ax.annotate(
                    hap,
                    xy=(0, y),
                    xycoords=("axes fraction", "data"),
                    xytext=(-3, -7),
                    textcoords="offset points",
                    ha="right",
                    va="top",
                    fontsize=6.5,
                    color="#444444",
                    linespacing=1.15,
                    annotation_clip=False,
                )
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", labelsize=8, which="minor", length=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    _style_x(ax, effect, xlim)
    ax.set_axisbelow(True)


def _legend_handles(groups, colors, p_values: Sequence[float], alpha: Optional[float]):
    color_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markerfacecolor=colors[group],
            markeredgecolor=colors[group],
            markersize=7,
            label=group,
        )
        for group in groups
    ]
    # Fixed anchors so the legend is comparable across figures. Drop the
    # saturated end when nothing in the data is that small. Largest first.
    anchors = [1e-20, 1e-10, 1e-5, 0.05]
    smallest = float(np.min(p_values)) if len(p_values) else 0.05
    if smallest > 1e-12:
        anchors = anchors[1:]
    size_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markerfacecolor="#4a4a4a",
            markeredgecolor="#4a4a4a",
            markersize=math.sqrt(float(marker_area(p))),
            label=format_p(p),
        )
        for p in anchors
    ]
    if alpha is not None:
        size_handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="None",
                markerfacecolor="white",
                markeredgecolor="#4a4a4a",
                markeredgewidth=1.05,
                markersize=7,
                label="p ≥ %s" % format_p(alpha),
            )
        )
    return color_handles, size_handles


def _style_legend(legend) -> None:
    legend.get_title().set_fontsize(10)
    for text in legend.get_texts():
        text.set_fontsize(8.5)


def _add_legends(host, color_handles, size_handles) -> None:
    # One block in the spare cell: ancestry, then the p-value size key.
    # Anchors are axes fractions so a tall panel does not pin them apart.
    ancestry = host.legend(
        handles=color_handles,
        title="Ancestry",
        loc="lower left",
        bbox_to_anchor=(0.06, 0.54),
        frameon=False,
        borderaxespad=0.0,
        labelspacing=0.35,
        handletextpad=0.45,
    )
    host.add_artist(ancestry)
    p_value = host.legend(
        handles=size_handles,
        title="P-value",
        loc="upper left",
        bbox_to_anchor=(0.06, 0.48),
        frameon=False,
        borderaxespad=0.0,
        labelspacing=0.55,
        handletextpad=0.65,
    )
    _style_legend(ancestry)
    _style_legend(p_value)


def _panel_title(selected: pd.DataFrame, phecode: str, labels: Optional[Mapping[str, str]]) -> str:
    """Phenotype name for a panel. An explicit label wins, then df['phenotype']."""
    if labels:
        titled = {str(key).strip(): value for key, value in labels.items()}
        if phecode in titled and str(titled[phecode]).strip():
            return str(titled[phecode]).strip()
    if "phenotype" in selected.columns:
        values = selected.loc[selected["phecode"] == phecode, "phenotype"]
        for value in values:
            if value is None or (isinstance(value, float) and math.isnan(value)):
                continue
            text = str(value).strip()
            if text and text.lower() not in MISSING:
                return text
    return phecode


def plot_metal_forest(
    df: pd.DataFrame,
    *,
    top_n: int = 8,
    effect: str = "or",
    ci: float = 0.95,
    studies: Optional[Sequence[str]] = None,
    phecodes: Optional[Sequence[str]] = None,
    labels: Optional[Mapping[str, str]] = None,
    ncols: int = 3,
    min_af: Optional[float] = None,
    min_n: Optional[float] = None,
    xlim: Optional[tuple[float, float]] = None,
    alpha: Optional[float] = None,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Draw one forest panel per phecode. Returns the figure and the points drawn.

    The panel title is the phenotype. A `phenotype` column on `df` is used
    when present. `labels` maps phecode to a title and overrides that column.
    `studies` is the ancestry order; the default is every ancestry column in
    the table. `phecodes` selects and orders panels. A `MarkerID` column is
    drawn under each allele name as the HI haplotype.

    `alpha` is the p-value threshold. An estimate with p ≥ alpha is an open
    circle; a smaller p stays filled. The default is 1 − `ci`, so a 95%
    interval marks p ≥ 0.05.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("plot_metal_forest expects a DataFrame")
    missing = [col for col in REQUIRED if col not in df.columns]
    if missing:
        raise ValueError("table is missing %s" % ", ".join(missing))
    if alpha is None:
        alpha = 1.0 - ci
    if not 0.0 < alpha <= 1.0:
        raise ValueError(f"alpha must be between 0 and 1, got {alpha}")
    table = df
    found = order_studies(available_studies(table.columns), studies)
    groups = ["META", *found]
    phecode_values = table["phecode"].astype(str).str.strip()
    if phecodes is None:
        phecode_order = list(dict.fromkeys(phecode_values.tolist()))
    else:
        phecode_order = [str(phecode).strip() for phecode in phecodes]
        missing = [phecode for phecode in phecode_order if phecode not in set(phecode_values)]
        if missing:
            raise ValueError("phecodes not in the table: %s" % ", ".join(missing))
    selected = select_top(table, found, phecode_order, top_n, min_af, min_n)
    plotted = forest_frame(selected, found, effect, ci)
    plotted["significant"] = plotted["p"] < alpha
    window = xlim if xlim is not None else _data_xlim(plotted, effect)
    colors = group_colors(groups)
    offsets = dodge_offsets(groups)
    panels = [phecode for phecode in phecode_order if (selected["phecode"] == phecode).any()]
    n_panels = len(panels)
    ncols = max(1, min(int(ncols), n_panels))
    nrows = math.ceil(n_panels / ncols)
    # A spare cell holds the legend (a 5-phenotype figure on 3 columns). A
    # full grid reserves a right margin instead.
    legend_in_grid = n_panels < nrows * ncols
    max_n = int(selected.groupby("phecode", sort=False).size().max())
    n_label_lines, longest_line = _label_metrics(selected)
    _pad, row_block = _panel_span(max_n, n_label_lines)
    # The left margin grows with the longest line of the wrapped haplotype.
    char_in = 7.0 / 72.0 * 0.56
    extra_label = max(0.0, longest_line * char_in + 0.35 - 1.15)
    width = ncols * 4.2 + extra_label + (0.0 if legend_in_grid else 2.4)
    height = nrows * (row_block + 0.95)
    with plt.rc_context(
        {
            "font.size": 9,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    ):
        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(width, height),
            squeeze=False,
        )
        fig.patch.set_facecolor("white")
        flat = list(axes.ravel())
        for ax, phecode in zip(flat, panels):
            panel = plotted.loc[plotted["phecode"] == phecode]
            _draw_panel(ax, panel, offsets, colors, effect, window, alpha)
            ax.set_title(_panel_title(selected, phecode, labels), fontsize=11, pad=6)
        spare = flat[n_panels:]
        color_handles, size_handles = _legend_handles(groups, colors, plotted["p"].tolist(), alpha)
        if spare:
            host = spare[0]
            host.set_axis_off()
            _add_legends(host, color_handles, size_handles)
            for ax in spare[1:]:
                ax.set_axis_off()
        fig.tight_layout()
        # The sizing labels have done their job. Drop them so the figure carries
        # only the allele name and the haplotype drawn under it.
        for ax in flat:
            if getattr(ax, "_haplograph_size_labels", False):
                for tick in ax.get_yticklabels():
                    tick.set_visible(False)
        if not spare:
            # Keep the legend outside the top-right panel. subplots_adjust
            # opens the margin; the anchor is in that panel's axes fraction.
            fig.subplots_adjust(right=0.76)
            anchor = flat[ncols - 1]
            ancestry = fig.legend(
                handles=color_handles,
                title="Ancestry",
                loc="upper left",
                bbox_to_anchor=(1.12, 1.02),
                bbox_transform=anchor.transAxes,
                frameon=False,
                borderaxespad=0.0,
                labelspacing=0.35,
                handletextpad=0.45,
            )
            p_value = fig.legend(
                handles=size_handles,
                title="P-value",
                loc="upper left",
                bbox_to_anchor=(1.12, 0.38),
                bbox_transform=anchor.transAxes,
                frameon=False,
                borderaxespad=0.0,
                labelspacing=0.55,
                handletextpad=0.65,
            )
            _style_legend(ancestry)
            _style_legend(p_value)
    plotted.attrs["groups"] = groups
    plotted.attrs["effect"] = effect
    plotted.attrs["ci"] = ci
    plotted.attrs["alpha"] = alpha
    return fig, plotted


def _first_column(columns: Sequence[str], candidates: Sequence[str], exclude: Sequence[str] = ()) -> Optional[str]:
    have = set(columns)
    for name in candidates:
        if name in have and name not in exclude:
            return name
    return None


def _clean_text(value) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    text = str(value).strip()
    return "" if text in MISSING else text


def _format_p_cell(p: float) -> str:
    """P-value for the text column: 0.012, or 3.2×10⁻¹² in mathtext."""
    if p >= 0.001:
        return f"{p:.3f}"
    if p <= 0:
        return "<1×10$^{-300}$"
    mantissa, exponent = f"{p:.1e}".split("e")
    return f"{mantissa}×10$^{{{int(exponent)}}}$"


def _format_estimate_cell(estimate: float, low: float, high: float) -> str:
    return f"{estimate:.2f} ({low:.2f}–{high:.2f})"


def plot_allele_phewas_forest(
    df: pd.DataFrame,
    allele: str,
    *,
    top_n: int = 15,
    effect: str = "or",
    ci: float = 0.95,
    alpha: float = 0.05,
    marker_col: str = "MarkerID",
    phecode_col: Optional[str] = None,
    label_col: Optional[str] = None,
    p_bonf_col: str = "p_bonferroni",
    labels: Optional[Mapping[str, str]] = None,
    title: Optional[str] = None,
    show_phecode: bool = True,
    category_col: Optional[str] = None,
    colors: Optional[Mapping[str, str]] = None,
    significant_only: bool = False,
    annotate: bool = True,
    xlim: Optional[tuple[float, float]] = None,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Forest plot of one allele across phenotypes. Returns the figure and the rows drawn.

    `df` is stacked SAIGE single-variant output for the mega-analysis: one
    row per marker and phenotype, with BETA, SE, and p.value (BETA is the
    effect of Allele2). The rows whose `marker_col` equals `allele` are
    ranked by p and the `top_n` phenotypes are drawn, the smallest p on top.
    The title is `title`, else `allele`.

    A point is filled when p_bonferroni < `alpha` and open otherwise. When
    the table has no `p_bonf_col`, it is p times the number of phenotypes
    tested for this allele, capped at 1, which is only right when `df` holds
    every test for the allele. `significant_only` keeps the filled rows.
    Point area grows with -log10(p).

    The row label is `labels[phecode]`, else `label_col`, else the phecode.
    `phecode_col` defaults to the first of phecode / phenotype / phenocode /
    trait; `label_col` to the first description-like column that is not the
    phecode. `show_phecode` appends the phecode to a descriptive label.
    `category_col` (e.g. pheno_category) colors points; `colors` maps a
    category to a color. `annotate` writes the estimate, CI, and p on the right.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("plot_allele_phewas_forest expects a DataFrame")
    if effect not in ("or", "beta"):
        raise ValueError(f"effect must be 'or' or 'beta', got {effect}")
    if not 0.0 < alpha <= 1.0:
        raise ValueError(f"alpha must be between 0 and 1, got {alpha}")
    if top_n < 1:
        raise ValueError(f"top_n must be >= 1, got {top_n}")
    z = z_crit(ci)
    columns = list(df.columns)
    missing = [col for col in (marker_col, "BETA", "SE") if col not in df.columns]
    if missing:
        raise ValueError("table is missing %s" % ", ".join(missing))
    id_col = phecode_col or _first_column(columns, ("phecode", "phenotype", "phenocode", "trait"))
    if id_col is None or id_col not in df.columns:
        raise ValueError("no phenotype column (set phecode_col)")
    if label_col is None:
        label_col = _first_column(
            columns,
            ("description", "phecode_string", "phenotype_description", "phenotype_name", "phenotype"),
            exclude=(id_col,),
        )
    elif label_col not in df.columns:
        raise ValueError(f"label_col {label_col} is not in the table")
    if category_col is not None and category_col not in df.columns:
        raise ValueError(f"category_col {category_col} is not in the table")
    p_col = resolve_p_column(columns)

    target = str(allele).strip()
    work = df.loc[df[marker_col].astype(str).str.strip() == target].copy()
    if work.empty:
        raise ValueError(f"{target} is not in {marker_col}")
    work["phecode"] = work[id_col].map(_clean_text)
    work["beta"] = parse_numeric(work["BETA"])
    work["se"] = parse_numeric(work["SE"])
    work["p"] = parse_numeric(work[p_col])
    keep = np.isfinite(work["beta"]) & np.isfinite(work["se"]) & (work["se"] > 0)
    keep &= work["p"].between(0.0, 1.0)
    keep &= work["phecode"] != ""
    work = work.loc[keep].copy()
    if work.empty:
        raise ValueError(f"{target} has no row with a finite BETA, SE, and p-value")
    work["_abs_beta"] = work["beta"].abs()
    work = work.sort_values(
        ["p", "_abs_beta", "phecode"],
        ascending=[True, False, True],
        kind="mergesort",
    ).drop_duplicates("phecode", keep="first")
    if p_bonf_col in work.columns:
        work["p_bonferroni"] = parse_numeric(work[p_bonf_col])
    else:
        work["p_bonferroni"] = np.minimum(work["p"] * len(work), 1.0)
    work["significant"] = work["p_bonferroni"] < alpha
    if significant_only:
        work = work.loc[work["significant"]]
        if work.empty:
            raise ValueError(f"{target} has no phenotype with {p_bonf_col} < {alpha}")
    shown = work.head(top_n).copy()
    n_rows = len(shown)
    shown["rank"] = np.arange(1, n_rows + 1)
    shown["y"] = n_rows - shown["rank"]
    low = shown["beta"] - z * shown["se"]
    high = shown["beta"] + z * shown["se"]
    if effect == "or":
        shown["estimate"] = np.exp(shown["beta"])
        shown["ci_low"] = np.exp(low)
        shown["ci_high"] = np.exp(high)
    else:
        shown["estimate"] = shown["beta"]
        shown["ci_low"] = low
        shown["ci_high"] = high

    titled = {str(key).strip(): str(value).strip() for key, value in (labels or {}).items()}
    row_labels = []
    for record in shown.to_dict("records"):
        phecode = record["phecode"]
        text = titled.get(phecode) or (_clean_text(record.get(label_col)) if label_col else "")
        if not text:
            text = phecode
        elif show_phecode and text != phecode:
            text = f"{text} ({phecode})"
        row_labels.append(text)
    shown["label"] = row_labels
    wrapped = [wrap_label(text) for text in row_labels]

    if category_col is not None:
        shown["category"] = shown[category_col].map(_clean_text).replace("", "Other")
        categories = list(dict.fromkeys(shown["category"].tolist()))
        palette = dict(colors or {})
        unknown = 0
        for category in categories:
            if category not in palette:
                palette[category] = SERIES_COLORS[unknown % len(SERIES_COLORS)]
                unknown += 1
    else:
        shown["category"] = ""
        categories = []
        palette = {"": GROUP_COLORS["META"]}

    window = xlim if xlim is not None else _full_xlim(shown, effect)
    color_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markerfacecolor=palette[category],
            markeredgecolor=palette[category],
            markersize=7,
            label=category,
        )
        for category in categories
    ]
    _unused, size_handles = _legend_handles([], {}, shown["p"].tolist(), alpha)
    size_handles[-1].set_label("P$_{Bonf}$ ≥ %s" % format_p(alpha))

    # Fixed inches so the row pitch and text column do not depend on n_rows.
    max_lines = max(text.count("\n") + 1 for text in wrapped)
    longest = max(len(line) for text in wrapped for line in text.split("\n"))
    row_in = 0.30 + 0.13 * (max_lines - 1)
    legend_in = 0.35 + 0.27 * len(size_handles)
    if color_handles:
        legend_in += 0.45 + 0.2 * len(color_handles)
    ax_h = max(row_in * n_rows, legend_in, 1.2)
    ax_w = 3.4
    left = 0.3 + longest * 8.0 / 72.0 * 0.55
    text_w = 2.0 if annotate else 0.0
    legend_w = 1.8
    top, bottom = 0.75 if annotate else 0.55, 0.55
    width = left + ax_w + text_w + legend_w
    height = ax_h + top + bottom
    with plt.rc_context(
        {
            "font.size": 9,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    ):
        fig = plt.figure(figsize=(width, height))
        fig.patch.set_facecolor("white")
        ax = fig.add_axes([left / width, bottom / height, ax_w / width, ax_h / height])
        for y in shown["y"]:
            ax.axhline(y, color="#d5d5d5", linestyle=(0, (1, 1.4)), linewidth=0.8, zorder=0)
        for category in dict.fromkeys(shown["category"]):
            sub = shown.loc[shown["category"] == category]
            x = sub["estimate"].to_numpy(dtype=float)
            lo = sub["ci_low"].to_numpy(dtype=float)
            hi = sub["ci_high"].to_numpy(dtype=float)
            y = sub["y"].to_numpy(dtype=float)
            color = palette[category]
            ax.errorbar(
                x,
                y,
                xerr=np.vstack([np.maximum(x - lo, 0.0), np.maximum(hi - x, 0.0)]),
                fmt="none",
                ecolor=color,
                elinewidth=0.9,
                capsize=2.0,
                capthick=0.8,
                zorder=2,
            )
            _scatter_group(
                ax,
                x,
                y,
                sub["p"].to_numpy(dtype=float),
                color,
                alpha,
                filled=sub["significant"].to_numpy(dtype=bool),
            )
        ax.set_yticks(shown["y"].tolist())
        ax.set_yticklabels(wrapped, fontsize=8, linespacing=1.1)
        ax.set_ylim(-0.6, n_rows - 0.4)
        ax.tick_params(axis="y", length=0, pad=4)
        ax.tick_params(axis="x", labelsize=8, which="minor", length=0)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        _style_x(ax, effect, window)
        ax.set_axisbelow(True)
        ax.set_title(title or target, fontsize=11, fontstyle="italic", pad=22 if annotate else 6)

        if annotate:
            est_x = 1.0 + 0.15 / ax_w
            p_x = 1.0 + 1.45 / ax_w
            header = "OR (%d%% CI)" % round(ci * 100) if effect == "or" else "β (%d%% CI)" % round(ci * 100)
            for x_frac, text in ((est_x, header), (p_x, "P")):
                ax.text(x_frac, 1.0, text, transform=ax.transAxes, ha="left", va="bottom", fontsize=8, fontweight="bold")
            row_transform = ax.get_yaxis_transform()
            for record in shown.to_dict("records"):
                weight = "bold" if record["significant"] else "normal"
                ax.text(
                    est_x,
                    record["y"],
                    _format_estimate_cell(record["estimate"], record["ci_low"], record["ci_high"]),
                    transform=row_transform,
                    ha="left",
                    va="center",
                    fontsize=7.5,
                    clip_on=False,
                )
                ax.text(
                    p_x,
                    record["y"],
                    _format_p_cell(record["p"]),
                    transform=row_transform,
                    ha="left",
                    va="center",
                    fontsize=7.5,
                    fontweight=weight,
                    clip_on=False,
                )

        legend_x = (left + ax_w + text_w + 0.2) / width
        legend_top = (bottom + ax_h) / height
        if color_handles:
            category_legend = fig.legend(
                handles=color_handles,
                title="Category",
                loc="upper left",
                bbox_to_anchor=(legend_x, legend_top),
                frameon=False,
                borderaxespad=0.0,
                labelspacing=0.35,
                handletextpad=0.45,
            )
            _style_legend(category_legend)
            legend_top -= (0.45 + 0.2 * len(color_handles)) / height
        p_legend = fig.legend(
            handles=size_handles,
            title="P-value",
            loc="upper left",
            bbox_to_anchor=(legend_x, legend_top),
            frameon=False,
            borderaxespad=0.0,
            labelspacing=0.55,
            handletextpad=0.65,
        )
        _style_legend(p_legend)

    keep_cols = [
        "rank",
        "y",
        "phecode",
        "label",
        "category",
        "beta",
        "se",
        "p",
        "p_bonferroni",
        "significant",
        "estimate",
        "ci_low",
        "ci_high",
    ]
    extra = [col for col in ("Allele1", "Allele2", "AF_Allele2", "N_case", "N_ctrl", "N") if col in shown.columns]
    plotted = shown[keep_cols + extra].reset_index(drop=True)
    plotted.insert(0, "MarkerID", target)
    plotted.attrs["effect"] = effect
    plotted.attrs["ci"] = ci
    plotted.attrs["alpha"] = alpha
    return fig, plotted
