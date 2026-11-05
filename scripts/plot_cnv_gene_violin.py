#!/usr/bin/env python3
"""Copy-number × ancestry violin plots for copy-number-variable genes.

Notebook:
    from plot_cnv_gene_violin import plot_cnv_gene_violin, rank_cnv_genes, tidy_copy_number

    # wide matrix (samples x genes) + per-sample ancestry
    fig, stats = plot_cnv_gene_violin(cn_wide, ancestry=df_meta["ancestry_pred"])

    # or an already-long frame with sample / gene / copy_number / ancestry columns
    fig, stats = plot_cnv_gene_violin(long_df)

    # the 12 most copy-number-variable genes only, as boxes instead of violins
    fig, stats = plot_cnv_gene_violin(cn_wide, anc, top_n=12, rank_by="var", kind="box")

    # or just the ranking, to pick genes yourself
    rank_cnv_genes(cn_wide, anc, metric="cv").head(20)

One subplot per gene, one violin (or box) per ancestry within a panel, with the
per-sample points jittered on top of it. Panels run most variable first when
`top_n` is used, by median copy number otherwise. `stats` is the
per (gene, ancestry) summary that was plotted (n, median, quartiles, mean, range,
fraction non-diploid), plus a Kruskal-Wallis p-value per gene comparing ancestry
groups — returned, never annotated onto the figure, so the caller decides what
multiple-testing correction to show.
"""

from __future__ import annotations

import itertools
import math
from typing import Mapping, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.container import ErrorbarContainer
from matplotlib.lines import Line2D

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

AncestryArg = Union[pd.Series, pd.DataFrame, Mapping, None]


def _ancestry_order(labels: Sequence[str]) -> list[str]:
    present = list(dict.fromkeys(str(x) for x in labels))
    known = [a for a in CANONICAL_ANCESTRY_ORDER if a in present]
    rest = sorted(set(present) - set(known))
    return known + rest


def _palette_for(order: Sequence[str]) -> dict[str, str]:
    cycle = itertools.cycle(OKABE_ITO)
    return {anc: next(cycle) for anc in order}


def _seaborn_version() -> tuple[int, ...]:
    parts = []
    for piece in str(sns.__version__).split("."):
        digits = "".join(itertools.takewhile(str.isdigit, piece))
        if not digits:
            break
        parts.append(int(digits))
    return tuple(parts)


# seaborn >= 0.13 colours categories through `hue` and needs legend=False to stay
# quiet; older versions have neither keyword and colour the x levels from `palette`
# directly (passing `legend` there lands in Axes.boxplot and raises TypeError).
_SNS_HUE_API = _seaborn_version() >= (0, 13)


def _color_kwargs(col: str, order: Sequence[str], palette: Mapping[str, str]) -> dict:
    if _SNS_HUE_API:
        return dict(hue=col, hue_order=list(order), palette=dict(palette), legend=False)
    return dict(palette=dict(palette))


def _density_kwargs(bw_adjust: float) -> dict:
    """Equal-width violins, with the kernel width named as this seaborn spells it.

    `gridsize` is the number of points the outline is evaluated on; the default
    100 shows visible facets on a narrow panel, so raise it.
    """
    if _SNS_HUE_API:
        return dict(density_norm="width", bw_adjust=bw_adjust, gridsize=300)
    # Pre-0.13 has no bw_adjust: `bw` is scipy's bw_method, where a scalar is the
    # factor itself rather than a multiple of Scott's. Scott's rule for these
    # group sizes lands near 0.3, so scale that to keep the two branches close.
    return dict(scale="width", bw=max(0.05, 0.3 * bw_adjust), gridsize=300)


def _mean_ci(values: np.ndarray, ci_level: float, n_boot: int, rng: np.random.Generator):
    """Mean and a percentile bootstrap CI. A one-value group has a zero-width interval."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return math.nan, math.nan, math.nan
    mean = float(v.mean())
    if v.size == 1 or n_boot <= 0:
        return mean, mean, mean
    boots = v[rng.integers(0, v.size, size=(n_boot, v.size))].mean(axis=1)
    tail = (100.0 - ci_level) / 2.0
    lo, hi = np.percentile(boots, [tail, 100.0 - tail])
    return mean, float(lo), float(hi)


def _as_ancestry_series(ancestry: AncestryArg, id_col: str, anc_col: str) -> pd.Series:
    """Sample -> ancestry label, as a Series indexed by sample ID (string)."""
    if isinstance(ancestry, pd.DataFrame):
        cols = {str(c).lower(): c for c in ancestry.columns}
        key = cols.get(id_col.lower())
        val = cols.get(anc_col.lower())
        if val is None:
            val = next(
                (cols[k] for k in ("ancestry", "ancestry_pred", "population", "pop", "superpop")
                 if k in cols),
                None,
            )
        if val is None:
            raise ValueError(
                f"ancestry DataFrame has no {anc_col!r} column; columns are {list(ancestry.columns)}"
            )
        if key is not None:
            s = pd.Series(ancestry[val].values, index=ancestry[key].astype(str).values)
        else:
            s = pd.Series(ancestry[val].values, index=ancestry.index.astype(str))
    elif isinstance(ancestry, pd.Series):
        s = pd.Series(ancestry.values, index=ancestry.index.astype(str))
    elif isinstance(ancestry, Mapping):
        s = pd.Series({str(k): v for k, v in ancestry.items()})
    else:
        raise TypeError(f"cannot read ancestry from {type(ancestry).__name__}")
    return s.dropna().astype(str)


def tidy_copy_number(
    copy_number: pd.DataFrame,
    ancestry: AncestryArg = None,
    *,
    genes: Optional[Sequence[str]] = None,
    sample_col: str = "sample",
    gene_col: str = "gene",
    cn_col: str = "copy_number",
    ancestry_col: str = "ancestry",
    genes_as_rows: bool = False,
) -> pd.DataFrame:
    """Return a long sample/gene/copy_number/ancestry frame.

    `copy_number` is either already long (it has `gene_col` and `cn_col`) or a
    wide matrix of samples x genes — pass `genes_as_rows=True` for genes x samples.
    `ancestry` is only needed when the frame does not already carry `ancestry_col`.
    """
    df = copy_number
    cols = {str(c).lower(): c for c in df.columns}
    is_long = gene_col.lower() in cols and cn_col.lower() in cols

    if is_long:
        gcol, vcol = cols[gene_col.lower()], cols[cn_col.lower()]
        scol = cols.get(sample_col.lower())
        out = pd.DataFrame(
            {
                sample_col: (df[scol] if scol is not None else df.index).astype(str),
                gene_col: df[gcol].astype(str),
                cn_col: pd.to_numeric(df[vcol], errors="coerce"),
            }
        )
        acol = cols.get(ancestry_col.lower()) or cols.get("ancestry_pred")
        if acol is not None:
            out[ancestry_col] = df[acol].astype(str).values
    else:
        wide = df.T if genes_as_rows else df
        out = (
            wide.rename_axis(index=sample_col, columns=gene_col)
            .apply(pd.to_numeric, errors="coerce")
            .stack()
            .rename(cn_col)
            .reset_index()
        )
        out[sample_col] = out[sample_col].astype(str)
        out[gene_col] = out[gene_col].astype(str)

    if ancestry_col not in out.columns:
        if ancestry is None:
            raise ValueError(
                f"no {ancestry_col!r} column in the copy-number frame; pass ancestry=..."
            )
        labels = _as_ancestry_series(ancestry, sample_col, ancestry_col)
        out[ancestry_col] = out[sample_col].map(labels)

    if genes is not None:
        wanted = [str(g) for g in genes]
        unknown = sorted(set(wanted) - set(out[gene_col]))
        if unknown:
            raise ValueError(f"genes not in the copy-number data: {unknown}")
        out = out[out[gene_col].isin(wanted)]
        out[gene_col] = pd.Categorical(out[gene_col], categories=wanted, ordered=True)

    return out.dropna(subset=[cn_col, ancestry_col])


def _summarise(
    df: pd.DataFrame, gene_col: str, ancestry_col: str, cn_col: str, diploid_cn: float
) -> pd.DataFrame:
    rows = []
    for (gene, anc), g in df.groupby([gene_col, ancestry_col], observed=True, sort=False):
        v = g[cn_col].to_numpy(dtype=float)
        rows.append(
            {
                gene_col: gene,
                ancestry_col: anc,
                "n": v.size,
                "mean": v.mean(),
                "median": np.median(v),
                "q1": np.percentile(v, 25),
                "q3": np.percentile(v, 75),
                "min": v.min(),
                "max": v.max(),
                "frac_non_diploid": float(np.mean(~np.isclose(v, diploid_cn))),
            }
        )
    stats = pd.DataFrame(rows)

    pvals = {}
    for gene, g in df.groupby(gene_col, observed=True, sort=False):
        groups = [
            v[cn_col].to_numpy(dtype=float)
            for _, v in g.groupby(ancestry_col, observed=True)
            if len(v) >= 2
        ]
        p = math.nan
        if len(groups) >= 2:
            try:
                from scipy.stats import kruskal

                p = float(kruskal(*groups).pvalue)
            except Exception:  # scipy missing, or all values identical
                p = math.nan
        pvals[gene] = p
    stats["kruskal_p"] = stats[gene_col].map(pvals)
    return stats.sort_values([gene_col, ancestry_col]).reset_index(drop=True)


def _adjust_p(pvals: Sequence[float], method: str) -> np.ndarray:
    """Multiple-testing correction that keeps NaNs (untested comparisons) as NaN."""
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    m = int(ok.sum())
    if m == 0 or method == "none":
        out[ok] = p[ok]
        return out
    if method == "bonferroni":
        out[ok] = np.minimum(p[ok] * m, 1.0)
        return out
    if method != "fdr_bh":
        raise ValueError(f"p_adjust must be none, bonferroni or fdr_bh; got {method!r}")
    idx = np.argsort(p[ok], kind="stable")
    ranked = p[ok][idx] * m / np.arange(1, m + 1)
    # Step-up: an adjusted value is never larger than a less significant one.
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(m)
    q[idx] = np.minimum(ranked, 1.0)
    out[ok] = q
    return out


def _format_p(p: float, style: str) -> str:
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return "p = n/a"
    if style == "star":
        for cut, mark in ((1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (0.05, "*")):
            if p < cut:
                return mark
        return "ns"
    text = f"p = {p:.2g}" if p >= 1e-3 else f"p = {p:.1e}"
    if style == "value":
        return text
    if style == "both":
        return f"{text} {_format_p(p, 'star')}"
    raise ValueError(f"p_format must be value, star or both; got {style!r}")


def _mann_whitney(a: np.ndarray, b: np.ndarray) -> float:
    if min(a.size, b.size) < 3:
        return math.nan
    try:
        from scipy.stats import mannwhitneyu

        return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except Exception:  # scipy missing, or both groups constant and identical
        return math.nan


def _pairwise_rows(
    sub: pd.DataFrame,
    ancestry_col: str,
    cn_col: str,
    order: Sequence[str],
    ref: Optional[str],
) -> list[dict]:
    """Mann-Whitney U per ancestry pair (or each group against `ref`)."""
    values = {
        a: sub.loc[sub[ancestry_col] == a, cn_col].to_numpy(dtype=float) for a in order
    }
    if ref is not None:
        if ref not in values:
            raise ValueError(f"pairwise_ref {ref!r} is not among the plotted groups {list(order)}")
        pairs = [(ref, b) for b in order if b != ref]
    else:
        pairs = list(itertools.combinations(order, 2))
    return [
        {
            "group1": a,
            "group2": b,
            "n1": values[a].size,
            "n2": values[b].size,
            "p": _mann_whitney(values[a], values[b]),
        }
        for a, b in pairs
    ]


def _draw_bracket(ax, i: int, j: int, y: float, h: float, text: str) -> None:
    ax.plot(
        [i, i, j, j], [y, y + h, y + h, y], lw=0.9, color="0.25", clip_on=False, zorder=6
    )
    ax.text(
        (i + j) / 2,
        y + h,
        text,
        ha="center",
        va="bottom",
        fontsize=7.5,
        color="0.15",
        clip_on=False,
        zorder=6,
    )


def _annotate_panel(
    ax,
    gene: str,
    order: Sequence[str],
    values: np.ndarray,
    omnibus_q: float,
    pairs: pd.DataFrame,
    *,
    p_format: str,
    p_threshold: float,
    max_brackets: int,
    fixed_ylim: bool,
) -> None:
    """Kruskal-Wallis p in the corner, significant pairwise brackets above the data."""
    drawn = []
    if len(pairs):
        sig = pairs[pairs["q"].notna() & (pairs["q"] <= p_threshold)]
        drawn = list(sig.nsmallest(max_brackets, "q").itertuples(index=False))

    lo, hi = ax.get_ylim()
    top = float(np.nanmax(values)) if values.size else hi
    span = max(hi - lo, 1e-9)
    step = 0.09 * span
    y = top + 0.5 * step
    pos = {a: i for i, a in enumerate(order)}

    for row in drawn:
        i, j = sorted((pos[row.group1], pos[row.group2]))
        _draw_bracket(ax, i, j, y, 0.25 * step, _format_p(row.q, p_format))
        y += step

    show_kw = not (isinstance(omnibus_q, float) and math.isnan(omnibus_q))
    if not fixed_ylim:
        # Headroom for the brackets, plus a clear band above them for the KW label.
        needed = (y if drawn else top) + 0.4 * step
        if show_kw:
            needed += 0.16 * span
        ax.set_ylim(lo, max(hi, needed))

    if show_kw:
        ax.text(
            0.02,
            0.98,
            f"KW {_format_p(omnibus_q, p_format)}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.5,
            color="0.25",
            zorder=6,
        )


RANK_METRICS = ("var", "sd", "cv", "iqr", "frac_non_diploid", "ancestry_range")


def rank_cnv_genes(
    copy_number: pd.DataFrame,
    ancestry: AncestryArg = None,
    *,
    metric: str = "var",
    top_n: Optional[int] = None,
    sample_col: str = "sample",
    gene_col: str = "gene",
    cn_col: str = "copy_number",
    ancestry_col: str = "ancestry",
    genes_as_rows: bool = False,
    diploid_cn: float = 2.0,
) -> pd.DataFrame:
    """Rank genes by how copy-number-variable they are, most variable first.

    One row per gene with every metric computed, so you can re-sort without
    recomputing; `metric` only decides the sort (and what `top_n` cuts on):

      var / sd              variance / SD of copy number across samples
      cv                    SD / mean — compares genes with different baselines
      iqr                   q3 - q1, i.e. robust to a few extreme samples
      frac_non_diploid      fraction of samples away from `diploid_cn`
      ancestry_range        spread of the per-ancestry means (differentiation
                            between groups, not variability within them)

    `ancestry` is only needed for `ancestry_range`; without it that column is NaN.
    """
    if metric not in RANK_METRICS:
        raise ValueError(f"metric must be one of {RANK_METRICS}; got {metric!r}")

    try:
        df = tidy_copy_number(
            copy_number,
            ancestry,
            sample_col=sample_col,
            gene_col=gene_col,
            cn_col=cn_col,
            ancestry_col=ancestry_col,
            genes_as_rows=genes_as_rows,
        )
    except ValueError:
        if metric == "ancestry_range" or ancestry is not None:
            raise
        df = tidy_copy_number(
            copy_number,
            pd.Series("ALL", index=_sample_index(copy_number, genes_as_rows)),
            sample_col=sample_col,
            gene_col=gene_col,
            cn_col=cn_col,
            ancestry_col=ancestry_col,
            genes_as_rows=genes_as_rows,
        )

    by_gene = df.groupby(df[gene_col].astype(str), observed=True)[cn_col]
    out = pd.DataFrame(
        {
            "n": by_gene.size(),
            "mean": by_gene.mean(),
            "median": by_gene.median(),
            "min": by_gene.min(),
            "max": by_gene.max(),
            "sd": by_gene.std(ddof=1),
            "var": by_gene.var(ddof=1),
            "iqr": by_gene.quantile(0.75) - by_gene.quantile(0.25),
        }
    )
    out["cv"] = (out["sd"] / out["mean"]).replace([np.inf, -np.inf], np.nan)
    out["frac_non_diploid"] = (
        df.assign(_off=~np.isclose(df[cn_col].to_numpy(dtype=float), diploid_cn))
        .groupby(df[gene_col].astype(str), observed=True)["_off"]
        .mean()
    )
    # Means, not medians: copy numbers are integers, so per-ancestry medians tie
    # constantly and would rank genes by nothing but their arbitrary order.
    anc_mean = df.groupby(
        [df[gene_col].astype(str), df[ancestry_col]], observed=True
    )[cn_col].mean()
    out["ancestry_range"] = anc_mean.groupby(level=0).agg(lambda s: s.max() - s.min())

    out = out.rename_axis(gene_col).sort_values(
        metric, ascending=False, kind="stable", na_position="last"
    )
    return out.head(top_n) if top_n else out


def _sample_index(copy_number: pd.DataFrame, genes_as_rows: bool) -> pd.Index:
    wide = copy_number.T if genes_as_rows else copy_number
    return wide.index.astype(str)


def _gene_order(
    df: pd.DataFrame,
    gene_col: str,
    ancestry_col: str,
    cn_col: str,
    genes: Optional[Sequence[str]],
    order_by: str,
    ranking: Optional[pd.DataFrame] = None,
) -> list[str]:
    present = df[gene_col].astype(str)
    if genes is not None:
        return [g for g in (str(x) for x in genes) if (present == g).any()]
    if order_by == "rank":
        if ranking is None:
            raise ValueError("gene_order_by='rank' needs a ranking; pass top_n or rank_by")
        return [g for g in ranking.index.astype(str) if (present == g).any()]
    if order_by == "name":
        return sorted(present.unique())
    if order_by == "median":
        key = df.groupby(present, observed=True)[cn_col].median()
    elif order_by == "range":
        per_group = df.groupby([present, df[ancestry_col]], observed=True)[cn_col].mean()
        key = per_group.groupby(level=0).agg(lambda s: s.max() - s.min())
    else:
        raise ValueError(
            f"gene_order_by must be rank, median, range or name; got {order_by!r}"
        )
    return list(key.sort_values(kind="stable").index.astype(str))


def plot_cnv_gene_violin(
    copy_number: pd.DataFrame,
    ancestry: AncestryArg = None,
    *,
    genes: Optional[Sequence[str]] = None,
    top_n: Optional[int] = None,
    rank_by: str = "var",
    kind: str = "violin",
    sample_col: str = "sample",
    gene_col: str = "gene",
    cn_col: str = "copy_number",
    ancestry_col: str = "ancestry",
    genes_as_rows: bool = False,
    min_group_n: int = 3,
    show_points: bool = True,
    gene_order_by: Optional[str] = None,
    ncols: Optional[int] = None,
    panel_size: tuple[float, float] = (2.9, 2.9),
    figsize: Optional[tuple[float, float]] = None,
    point_size: float = 3.0,
    violin_alpha: float = 0.55,
    bw_adjust: float = 1.1,
    ci_level: float = 95,
    n_boot: int = 1000,
    ci_seed: Optional[int] = 0,
    annotate_p: str = "omnibus",
    p_adjust: str = "fdr_bh",
    p_format: str = "value",
    pairwise_ref: Optional[str] = None,
    p_threshold: float = 0.05,
    max_brackets: int = 6,
    sharey: bool = False,
    ylim: Optional[tuple[float, float]] = None,
    diploid_cn: float = 2.0,
    title: Optional[str] = "Gene copy number by ancestry",
    out_base: Optional[str] = "CNV_gene_copy_number_by_ancestry",
) -> tuple[plt.Figure, pd.DataFrame]:
    """Violin + jitter of per-sample copy number, one subplot per gene.

    Parameters
    ----------
    copy_number : DataFrame
        Long (sample/gene/copy_number) or wide (samples x genes) copy numbers.
    ancestry : Series / DataFrame / dict, optional
        Sample -> ancestry label. Not needed if `copy_number` already has the column.
    genes : sequence, optional
        Genes to plot, in this order (overrides `top_n` and `gene_order_by`).
    top_n : int, optional
        Plot only the `top_n` most copy-number-variable genes, ranked by `rank_by`
        (see `rank_cnv_genes`); panels then run most variable first. Default: all
        genes. Use this when the matrix has more genes than fit on one figure.
    rank_by : {"var", "sd", "cv", "iqr", "frac_non_diploid", "ancestry_range"}
        Which measure of "variable" `top_n` cuts on. "var"/"sd" favour genes with
        large absolute swings, "cv" normalises by the gene's own mean, "iqr" ignores
        a handful of outlying samples, "ancestry_range" picks genes that differ
        between ancestries rather than within them.
    kind : {"violin", "box", "point"}
        Violin (KDE), box, or a point plot of the group mean with a bootstrap CI.
        All three keep the jittered points. A box also draws the mean and its
        bootstrap CI on the box's centre line and joins the means across
        ancestries, because integer copy numbers often collapse the box onto one
        value and hide the ancestry shift — the slope of that line is the shift.
        "point" keeps every sample as a tab10 point and draws the mean as a
        white-filled dot with a capped bootstrap CI. The y-axis is the CI
        itself, padded by 0.5 on each side, so a short interval is readable;
        pass `ylim` to override that window. The default title for that kind
        names the interval.
    min_group_n : int
        Ancestry groups with fewer samples (for a given gene) are dropped, since a
        density over one or two points is not readable.
    show_points : bool
        Overlay jittered per-sample points on each violin.
    ci_level, n_boot, ci_seed
        Confidence level (95 by default), bootstrap resamples, and the seed that
        makes the interval reproducible (pass None for a fresh bootstrap each
        call). Used by `kind="point"` and by the mean marker on `kind="box"`.
    violin_alpha, bw_adjust : float
        Fill opacity, and the KDE bandwidth scale. Copy numbers are integers, so a
        narrow kernel turns each integer into its own bump and the outline into a
        comb; the default 1.1 smooths over that. Lower it (0.6-0.7) when you want
        the individual copy-number modes resolved instead.
    gene_order_by : {"rank", "median", "range", "name"}, optional
        Panel order when `genes` is not given: by the `rank_by` metric (most variable
        first), by median copy number, by the spread of ancestry means, or
        alphabetically. Defaults to "rank" with `top_n`, "median" without it.
    ncols, panel_size, figsize, sharey
        Grid layout. `figsize` overrides `panel_size x grid`; `sharey=True` is worth
        it when the genes are on a comparable copy-number scale.
    diploid_cn : float
        Copy number treated as the diploid baseline: drawn as a dashed reference
        line and used for the `frac_non_diploid` column of the summary.
    annotate_p : {"omnibus", "pairwise", "none"}
        "omnibus" prints one Kruskal-Wallis p per panel (does copy number differ
        across ancestries at all?). "pairwise" adds Mann-Whitney U brackets between
        ancestry pairs on top of it. Needs scipy; without it the p-values are NaN
        and nothing is drawn.
    p_adjust : {"fdr_bh", "bonferroni", "none"}
        Multiple-testing correction, applied across all panels of the figure at
        once — omnibus p-values over genes, pairwise p-values over every pair in
        every panel. Annotated values are the adjusted ones; both raw and adjusted
        are in the returned tables.
    p_format : {"value", "star", "both"}
        How the annotation reads: "p = 0.0012", "**", or both.
    pairwise_ref : str, optional
        With `annotate_p="pairwise"`, compare every ancestry against this one group
        instead of testing all pairs — far fewer brackets and a cleaner figure.
    p_threshold, max_brackets
        Only pairs reaching `p_threshold` (after adjustment) get a bracket, and at
        most `max_brackets` per panel, most significant first. Raise `p_threshold`
        to 1.0 to bracket everything tested.
    out_base : str, optional
        Writes `<out_base>.pdf`, `<out_base>.png` and `<out_base>.tsv` (the summary).
        Pass None to skip writing.

    Returns
    -------
    (figure, summary) — summary has one row per (gene, ancestry) plus the per-gene
    Kruskal-Wallis p (`kruskal_p`) and its adjusted value (`kruskal_q`). When
    pairwise tests run, the full pair table (including pairs too weak to draw) is
    `summary.attrs["pairwise"]`, and is written to `<out_base>.pairwise.tsv`.
    """
    df = tidy_copy_number(
        copy_number,
        ancestry,
        genes=genes,
        sample_col=sample_col,
        gene_col=gene_col,
        cn_col=cn_col,
        ancestry_col=ancestry_col,
        genes_as_rows=genes_as_rows,
    )
    if df.empty:
        raise ValueError("no (copy number, ancestry) pairs left after dropping missing values")

    if min_group_n > 1:
        sizes = df.groupby([gene_col, ancestry_col], observed=True)[cn_col].transform("size")
        df = df[sizes >= min_group_n]
        if df.empty:
            raise ValueError(f"no gene x ancestry group has at least {min_group_n} samples")

    if kind not in ("violin", "box", "point"):
        raise ValueError(f"kind must be 'violin', 'box' or 'point'; got {kind!r}")
    if kind == "point" and title == "Gene copy number by ancestry":
        title = f"Mean copy number by ancestry ({ci_level:g}% CI)"
    if annotate_p not in ("omnibus", "pairwise", "none"):
        raise ValueError(
            f"annotate_p must be omnibus, pairwise or none; got {annotate_p!r}"
        )

    ranking = None
    if genes is None and (top_n or gene_order_by == "rank"):
        ranking = rank_cnv_genes(
            df,
            metric=rank_by,
            top_n=top_n,
            sample_col=sample_col,
            gene_col=gene_col,
            cn_col=cn_col,
            ancestry_col=ancestry_col,
            diploid_cn=diploid_cn,
        )
        df = df[df[gene_col].astype(str).isin(set(ranking.index.astype(str)))]
        if df.empty:
            raise ValueError(f"no gene survived top_n={top_n} on {rank_by!r}")

    if gene_order_by is None:
        gene_order_by = "rank" if ranking is not None else "median"
    gene_order = _gene_order(
        df, gene_col, ancestry_col, cn_col, genes, gene_order_by, ranking
    )
    anc_order = _ancestry_order(df[ancestry_col])
    # Point panels color the samples with tab10 and keep the mean black.
    palette = "tab10" if kind == "point" else _palette_for(anc_order)
    n = len(gene_order)

    ncols = ncols or min(4, n)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=figsize or (panel_size[0] * ncols, panel_size[1] * nrows),
        sharey=sharey or bool(ylim),
        squeeze=False,
    )

    panels: list[tuple] = []
    pairwise_rows: list[dict] = []
    point_windows: list[tuple[float, float]] = []

    for ax, gene in itertools.zip_longest(axes.flat, gene_order):
        if gene is None:
            ax.set_visible(False)
            continue
        sub = df[df[gene_col].astype(str) == gene]
        order = [a for a in anc_order if (sub[ancestry_col] == a).any()]
        panels.append((ax, gene, order, sub))

        if annotate_p == "pairwise":
            for row in _pairwise_rows(sub, ancestry_col, cn_col, order, pairwise_ref):
                pairwise_rows.append({gene_col: gene, **row})

        # The violin/box carries the shape, the points the sample-level spread, and
        # the IQR bar with its white median marker goes on top of both.
        common = dict(
            data=sub,
            x=ancestry_col,
            y=cn_col,
            order=order,
            linewidth=0.8,
            ax=ax,
            zorder=2,
        )
        if kind != "point":
            common.update(_color_kwargs(ancestry_col, order, palette))
        if kind == "violin":
            first_coll = len(ax.collections)
            sns.violinplot(
                **common, **_density_kwargs(bw_adjust), cut=0, inner=None, width=0.9
            )
            for coll in ax.collections[first_coll:]:
                coll.set_alpha(violin_alpha)
                coll.set_edgecolor("0.35")
                coll.set_joinstyle("round")  # no mitre spikes at the integer modes
        elif kind == "box":
            # Narrow and pale: the mean ± CI is drawn on the box's own centre
            # line, so the fill has to stay quiet enough to read through.
            sns.boxplot(
                **common,
                width=0.6,
                fliersize=0 if show_points else 2,
                boxprops=dict(alpha=0.45),
                medianprops=dict(color="0.15", linewidth=1.6),
                whiskerprops=dict(color="0.55"),
                capprops=dict(color="0.55"),
            )
        # The point estimate is drawn by hand below, for every kind that needs
        # one: seaborn's pointplot keyword names for the marker, the CI and the
        # error bar all differ between 0.12 and 0.13, and none of them can be
        # sized in data space.

        if show_points:
            point_color = (
                dict(hue=ancestry_col, hue_order=order, palette=palette, legend=False)
                if kind == "point" and _SNS_HUE_API
                else dict(palette=palette)
                if kind == "point"
                else _color_kwargs(ancestry_col, order, palette)
            )
            sns.stripplot(
                data=sub,
                x=ancestry_col,
                y=cn_col,
                order=order,
                **point_color,
                jitter=0.22,
                size=3.5 if kind == "point" else point_size,
                alpha=0.8 if kind == "point" else 0.85,
                linewidth=0.0 if kind == "point" else 0.3,
                edgecolor="white",
                ax=ax,
                zorder=2 if kind == "point" else 3,
            )

        if ax.get_legend() is not None:  # older seaborn draws one per axes
            ax.get_legend().remove()

        # Violin: an IQR bar with a white median dot, since inner=None left the
        # shape without any summary. Box: the box keeps the median, and the mean
        # with its CI sits on the same centre line, joined across ancestries so a
        # shift the (often flat, integer-valued) box hides reads as a slope.
        # Point: the mean and its CI are the whole plot. All three sit at the tick.
        mean_x, mean_y, ci_lo, ci_hi = [], [], [], []
        rng = np.random.default_rng(ci_seed)
        for i, anc in enumerate(order):
            v = sub.loc[sub[ancestry_col] == anc, cn_col].to_numpy(dtype=float)
            if kind == "violin":
                q1, med, q3 = np.percentile(v, [25, 50, 75])
                ax.vlines(i, q1, q3, color="0.15", linewidth=3.0, zorder=4)
                ax.plot(
                    i, med, marker="o", markersize=4.5, color="white",
                    markeredgecolor="0.15", markeredgewidth=0.8, zorder=5,
                )
            if kind == "violin":
                continue
            mean, lo, hi = _mean_ci(v, ci_level, n_boot, rng)
            if not np.isfinite(mean):
                continue
            ci_lo.append(lo)
            ci_hi.append(hi)
            x = float(i)
            mean_x.append(x)
            mean_y.append(mean)
            ax.errorbar(
                x,
                mean,
                yerr=None if hi <= lo else [[mean - lo], [hi - mean]],
                fmt="o",
                color="0.1",
                ecolor="0.1",
                elinewidth=1.6,
                capsize=3.5,
                capthick=1.4,
                markersize=5.5,
                markerfacecolor="white",
                markeredgecolor="0.1",
                markeredgewidth=1.3,
                zorder=7,
            )
        # The line is the point of the box variant: it turns six separate means
        # into one readable trend. On a point plot the groups are unordered
        # categories, so joining them would invent a trend that isn't there.
        if kind == "box" and len(mean_x) > 1:
            ax.plot(mean_x, mean_y, color="0.1", linewidth=1.15, zorder=6)

        ax.axhline(diploid_cn, color="0.45", linestyle="--", linewidth=0.9, zorder=1)
        ax.set_title(gene, fontweight="bold", fontsize=11)
        ax.set_xlabel("")
        ax.set_ylabel("Copy number", fontweight="bold")
        counts = sub.groupby(ancestry_col, observed=True)[sample_col].nunique()
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(
            [f"{a}\n{counts[a]}" for a in order], fontsize=8, linespacing=1.4
        )
        ax.tick_params(axis="x", length=0)
        ax.set_xlim(-0.6, len(order) - 0.4)
        if ylim:
            ax.set_ylim(*ylim)
        elif kind == "point" and ci_lo:
            # Zoom to the estimates. The sample cloud spans whole copy numbers
            # and would otherwise swallow an interval of a few hundredths.
            # A shared y-axis uses the union of these windows, applied once
            # after every panel has reported its interval.
            window = (min(ci_lo) - 0.5, max(ci_hi) + 0.5)
            point_windows.append(window)
            if not sharey:
                ax.set_ylim(*window)
        ax.yaxis.grid(True, which="major", color="0.85", linewidth=0.6)
        ax.xaxis.grid(True, color="0.88", linewidth=0.8)
        ax.set_axisbelow(True)
        sns.despine(ax=ax)

    if sharey and point_windows and ylim is None:
        axes.flat[0].set_ylim(
            min(lo for lo, hi in point_windows),
            max(hi for lo, hi in point_windows),
        )

    stats = _summarise(df, gene_col, ancestry_col, cn_col, diploid_cn)
    per_gene_p = stats.drop_duplicates(gene_col).set_index(gene_col)["kruskal_p"]
    per_gene_q = pd.Series(
        _adjust_p(per_gene_p.to_numpy(), p_adjust), index=per_gene_p.index
    )
    stats["kruskal_q"] = stats[gene_col].map(per_gene_q)

    pairs = pd.DataFrame(pairwise_rows)
    if not pairs.empty:
        # Corrected over every pair in the figure, not panel by panel.
        pairs["q"] = _adjust_p(pairs["p"].to_numpy(), p_adjust)
        pairs = pairs.sort_values([gene_col, "q"], kind="stable").reset_index(drop=True)
    stats.attrs["pairwise"] = pairs

    if annotate_p != "none":
        for ax, gene, order, sub in panels:
            values = sub[cn_col].to_numpy(dtype=float)
            _annotate_panel(
                ax,
                gene,
                order,
                values,
                float(per_gene_q.get(gene, math.nan)),
                pairs[pairs[gene_col] == gene] if not pairs.empty else pairs,
                p_format=p_format,
                p_threshold=p_threshold,
                max_brackets=max_brackets,
                fixed_ylim=ylim is not None or kind == "point",
            )

    # Only the leftmost panel of each row needs the y label.
    for row in axes:
        for ax in row[1:]:
            if ax.get_visible():
                ax.set_ylabel("")

    fig.supxlabel(f"Ancestry (sample count under each {kind})", fontsize=9, y=0.005)
    if title:
        fig.suptitle(title, fontweight="bold", y=0.995)
    if kind == "box":
        # Hang the key just under the title; the top of the grid is left clear
        # of both so the two don't land on the same line.
        fig.legend(
            handles=[
                Line2D(
                    [0], [0], color="0.25", linewidth=1.6, label="median (box)",
                ),
                Line2D(
                    [0], [0], marker="o", color="0.1", markerfacecolor="white",
                    markeredgecolor="0.1", markeredgewidth=1.3, markersize=6,
                    linewidth=1.6, label=f"mean ± {ci_level:g}% CI",
                ),
            ],
            loc="upper center",
            bbox_to_anchor=(0.5, 0.955),
            ncol=2,
            frameon=False,
            fontsize=8,
        )
        fig.tight_layout(rect=(0, 0.012, 1, 0.90))
    else:
        fig.tight_layout()

    if out_base:
        fig.savefig(f"{out_base}.pdf", bbox_inches="tight")
        fig.savefig(f"{out_base}.png", dpi=300, bbox_inches="tight")
        stats.to_csv(f"{out_base}.tsv", sep="\t", index=False, float_format="%.4g")
        if not pairs.empty:
            pairs.to_csv(
                f"{out_base}.pairwise.tsv", sep="\t", index=False, float_format="%.4g"
            )
    return fig, stats


def _synthetic(n_per: int = 120, seed: int = 0) -> tuple[pd.DataFrame, pd.Series]:
    """Wide CN matrix + ancestry, with a real AFR excess at C4A and none at HLA-DRB1."""
    rng = np.random.default_rng(seed)
    centers = {"AFR": 0.9, "AMR": 0.2, "EAS": 0.0, "EUR": 0.1, "MID": 0.3, "SAS": 0.2}
    rows, labels = {}, {}
    for anc, shift in centers.items():
        for i in range(n_per):
            sid = f"{anc}_{i:03d}"
            labels[sid] = anc
            rows[sid] = {
                "AMY1": float(max(0, rng.poisson(7))),  # most variable by var/sd/iqr
                "C4A": float(max(0, rng.poisson(2 + shift))),
                "C4B": float(max(0, rng.poisson(1.8))),
                "HLA-DRB1": 2.0,  # invariant: last on every metric
                "KIR2DL5A": float(rng.integers(0, 3)),
                "SMN1": float(2 + (rng.random() < 0.05)),  # nearly invariant
            }
    wide = pd.DataFrame.from_dict(rows, orient="index")
    return wide, pd.Series(labels, name="ancestry")


if __name__ == "__main__":
    import tempfile
    from pathlib import Path

    wide, anc = _synthetic()
    out_dir = Path(tempfile.mkdtemp(prefix="cnv_violin_"))
    out_base = str(out_dir / "CNV_gene_copy_number_by_ancestry")
    fig, stats = plot_cnv_gene_violin(wide, anc, out_base=out_base)

    for ext in ("pdf", "png", "tsv"):
        assert Path(f"{out_base}.{ext}").is_file(), ext
    assert set(stats["gene"]) == set(wide.columns), sorted(set(stats["gene"]))
    assert (stats["n"] == 120).all(), stats["n"].unique()

    c4a = stats.set_index(["gene", "ancestry"]).loc["C4A"]
    assert c4a.loc["AFR", "mean"] > c4a.loc["EAS", "mean"], c4a[["mean"]]

    # The long-form path must give the same summary as the wide-form one.
    long_df = tidy_copy_number(wide, anc)
    _, stats_long = plot_cnv_gene_violin(long_df, out_base=None)
    pd.testing.assert_frame_equal(stats, stats_long)

    # Ranking: the engineered extremes must land at the ends on every metric.
    ranked = rank_cnv_genes(wide, anc)
    assert ranked.index[0] == "AMY1", ranked.index.tolist()
    assert ranked.index[-1] == "HLA-DRB1", ranked.index.tolist()
    assert list(ranked.columns[-2:]) == ["frac_non_diploid", "ancestry_range"]
    assert math.isclose(ranked.loc["HLA-DRB1", "var"], 0.0), ranked.loc["HLA-DRB1"]
    assert rank_cnv_genes(wide, anc, metric="iqr").index[0] == "AMY1"
    by_range = rank_cnv_genes(wide, anc, metric="ancestry_range")
    assert by_range.index[0] == "C4A", by_range["ancestry_range"].to_dict()
    assert len(rank_cnv_genes(wide, anc, top_n=3)) == 3

    # top_n keeps the most variable genes, most variable first, boxes included.
    top_base = str(out_dir / "top3_box")
    fig_top, stats_top = plot_cnv_gene_violin(
        wide, anc, top_n=3, kind="box", out_base=top_base
    )
    plotted = list(dict.fromkeys(stats_top["gene"]))
    assert set(plotted) == set(ranked.index[:3]), plotted
    assert [t.get_text() for t in fig_top.axes[0].get_xticklabels()][0].startswith("AFR")
    assert fig_top.axes[0].get_title() == "AMY1", fig_top.axes[0].get_title()
    assert Path(f"{top_base}.png").is_file()
    plt.close(fig_top)

    # Benjamini-Hochberg: monotone, never below the raw p, NaNs preserved.
    raw = [0.001, 0.008, 0.039, 0.041, 0.9, math.nan]
    q = _adjust_p(raw, "fdr_bh")
    assert math.isnan(q[-1]) and np.all(np.diff(q[:-1]) >= -1e-12), q
    assert np.all(q[:-1] >= np.asarray(raw[:-1]) - 1e-12), q
    assert math.isclose(_adjust_p([0.01, 0.02], "bonferroni")[0], 0.02), "bonferroni"
    assert _format_p(1e-9, "star") == "****" and _format_p(0.4, "star") == "ns"

    # Omnibus p-values: the engineered C4A shift is significant, C4B is not.
    kw = stats.drop_duplicates("gene").set_index("gene")
    assert kw.loc["C4A", "kruskal_q"] < 1e-3, kw.loc["C4A", "kruskal_q"]
    assert kw.loc["C4B", "kruskal_q"] > 0.05, kw.loc["C4B", "kruskal_q"]
    tested = kw.dropna(subset=["kruskal_p"])
    assert (tested["kruskal_q"] >= tested["kruskal_p"] - 1e-12).all(), tested
    assert math.isnan(kw.loc["HLA-DRB1", "kruskal_q"]), "invariant gene stays untested"

    # Pairwise: AFR vs EAS at C4A is the comparison that must survive correction.
    pair_base = str(out_dir / "pairwise_vs_EAS")
    fig_p, stats_p = plot_cnv_gene_violin(
        wide, anc, genes=["C4A", "C4B"], annotate_p="pairwise",
        pairwise_ref="EAS", p_format="both", out_base=pair_base,
    )
    pw = stats_p.attrs["pairwise"]
    assert set(pw["group1"]) == {"EAS"}, sorted(set(pw["group1"]))
    hit = pw[(pw["gene"] == "C4A") & (pw["group2"] == "AFR")].iloc[0]
    assert hit["q"] < 0.01 and hit["q"] >= hit["p"], hit.to_dict()
    assert Path(f"{pair_base}.pairwise.tsv").is_file()
    # Brackets drawn only where the correction kept the pair significant.
    c4a_ax = next(a for a in fig_p.axes if a.get_title() == "C4A")
    c4b_ax = next(a for a in fig_p.axes if a.get_title() == "C4B")
    n_sig = int((pw[pw["gene"] == "C4A"]["q"] <= 0.05).sum())
    assert sum("p = " in t.get_text() for t in c4a_ax.texts) == min(n_sig, 6) + 1
    assert sum("p = " in t.get_text() for t in c4b_ax.texts) == 1, "omnibus only"
    plt.close(fig_p)

    _, stats_none = plot_cnv_gene_violin(wide, anc, annotate_p="none", out_base=None)
    assert stats_none.attrs["pairwise"].empty

    # Point plot: a mean marker per group, with a bootstrap CI that brackets it
    # and is far narrower than the sample spread.
    fig_pt, stats_pt = plot_cnv_gene_violin(
        wide, anc, genes=["C4A"], kind="point", show_points=False,
        annotate_p="none", out_base=None,
    )
    pt_ax = fig_pt.axes[0]
    bars = {}  # group position -> (lo, hi) of its error bar
    for container in pt_ax.containers:
        if not isinstance(container, ErrorbarContainer) or not container.lines[2]:
            continue
        (seg,) = container.lines[2][0].get_segments()
        bars[round(float(seg[0][0]))] = (float(seg[:, 1].min()), float(seg[:, 1].max()))
    assert sorted(bars) == [0, 1, 2, 3, 4, 5], sorted(bars)

    afr = stats_pt[stats_pt["ancestry"] == "AFR"].iloc[0]
    lo_ci, hi_ci = bars[0]
    assert lo_ci < afr["mean"] < hi_ci, (lo_ci, afr["mean"], hi_ci)
    assert hi_ci - lo_ci < afr["q3"] - afr["q1"], "CI must be narrower than the IQR"
    plt.close(fig_pt)

    # Point style: tab10 sample colors, a black mean, and a y-axis zoomed to
    # the confidence intervals (padded by 0.5), not the sample cloud.
    fig_lo, _ = plot_cnv_gene_violin(
        wide, anc, genes=["SMN1", "C4A"], kind="point", show_points=True,
        annotate_p="none", out_base=None,
    )
    assert fig_lo._suptitle.get_text() == "Mean copy number by ancestry (95% CI)"
    smn_ax = next(a for a in fig_lo.axes if a.get_title() == "SMN1")
    ylo, yhi = smn_ax.get_ylim()
    ci_span = [
        (float(seg[:, 1].min()), float(seg[:, 1].max()))
        for c in smn_ax.containers
        if isinstance(c, ErrorbarContainer) and c.lines[2]
        for (seg,) in [c.lines[2][0].get_segments()]
    ]
    ci_min = min(lo for lo, hi in ci_span)
    ci_max = max(hi for lo, hi in ci_span)
    assert math.isclose(ylo, ci_min - 0.5), (ylo, ci_min)
    assert math.isclose(yhi, ci_max + 0.5), (yhi, ci_max)
    assert yhi - ylo < 2, (ylo, yhi)
    means = [
        c.lines[0] for c in smn_ax.containers if isinstance(c, ErrorbarContainer)
    ]
    assert len(means) == 6, len(means)
    assert means[0].get_markerfacecolor() == "white", means[0].get_markerfacecolor()
    assert means[0].get_markersize() < 8, means[0].get_markersize()
    tab_blue = plt.get_cmap("tab10")(0)
    sample_colors = smn_ax.collections[0].get_facecolors()
    assert any(np.allclose(c[:3], tab_blue[:3], atol=0.02) for c in sample_colors)
    plt.close(fig_lo)

    try:
        plot_cnv_gene_violin(wide, anc, kind="bar", out_base=None)
    except ValueError as exc:
        assert "point" in str(exc), exc
    else:
        raise AssertionError("an unknown kind must be rejected")

    print(f"OK: wrote {out_base}.{{pdf,png,tsv}} and {top_base}.{{pdf,png,tsv}}")
    print(ranked.round(3).to_string())
    plt.close(fig)
