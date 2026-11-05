"""
hla_linkage_plot — Sankey/alluvial plot of classical HLA linkage.

Takes the DataFrame returned by `allele_pair_ld.pair_ld` and draws the classical
HLA loci as columns in genomic order, each allele as a node, and each allele pair
as a ribbon whose width is its r2. Reading left to right across the columns traces
the common haplotypes.

    from allele_pair_ld import pair_ld
    from hla_linkage_plot import plot_hla_linkage

    pairs = pair_ld(dosage_df, population=pops)
    fig = plot_hla_linkage(pairs, population="EUR")
    fig.savefig("hla_linkage_EUR.png", dpi=200, bbox_inches="tight")

matplotlib only — no plotly, no other dependency.

Reading the figure
------------------
Node HEIGHT is the allele's FREQUENCY in that population, so a tall block is a
common allele. Frequencies come from the `allele_freq_a` / `allele_freq_b` columns
`pair_ld` already emits; pass `node_size="carrier_freq"` for carrier frequency
instead, or `node_size="linkage"` for the old behaviour (height = total linked
weight).

Ribbon WIDTH is the pair's LD strength — whichever column you name in `weight`
(r2 by default, Dprime for strength of association). The key at the bottom
gives both scales. Ribbons run only between neighbouring columns by default,
because that is what a Sankey can place honestly — a link from HLA-A straight to
DQB1 would have to cross columns it has no relation to. Pass `adjacent_only=False`
to draw the long-range ones too, as arcs riding above the columns.

Height and width are now two independent quantities, so they carry two independent
scales. Frequencies set the node heights; the ribbon scale is then the largest one
that still lets every node's ribbons stack inside it. A ribbon of a given value is
therefore the same thickness everywhere in the figure, and a node with room to
spare simply shows slack — that slack is real, it is linkage the threshold or the
neighbouring column does not account for.

Columns are coloured by locus, in a palette validated for colour-vision deficiency
on the adjacent-pair gates. Every node and column is directly labelled, so identity
never rests on colour alone.
"""

from __future__ import annotations

import re
from typing import Iterable, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, PathPatch
from matplotlib.path import Path
from matplotlib.transforms import blended_transform_factory

__all__ = [
    "plot_hla_linkage",
    "CLASSICAL_HLA_ORDER",
    "HLA_ALL_ORDER",
    "KIR_ORDER",
    "LOCUS_PRESETS",
    "make_hla_demo_data",
]

# Classical class I then class II, in chromosome 6p21 order (A ~29.9 Mb -> DPB1 ~33.1 Mb).
CLASSICAL_HLA_ORDER: tuple = (
    "A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1",
)

# Every HLA and HLA-region gene worth drawing, in GRCh38 chr6 coordinate order
# (HLA-F ~29.72 Mb -> HLA-DPB2 ~33.11 Mb). Note C precedes B, and E sits between
# A and C rather than next to the other non-classical class I genes — the column
# order here is the chromosome, not the class.
HLA_ALL_ORDER: tuple = (
    # class I: non-classical and pseudogenes interleaved with the classical three
    "F", "G", "H", "A", "E", "C", "B", "MICA", "MICB",
    # class II: DR first, then DQ, DO/DM, DP
    "DRA", "DRB5", "DRB4", "DRB3", "DRB1",
    "DQA1", "DQB1", "DQA2", "DQB2",
    "DOB", "DMB", "DMA", "DOA",
    "DPA1", "DPB1", "DPA2", "DPB2",
)

# KIR at 19q13.4, centromeric -> telomeric along the reference haplotype. The
# gene content itself is haplotype-variable, so this is an order to lay columns
# out in, not a claim that any sample carries all of them. 2DS3 and 2DS5 can each
# sit at the centromeric or the telomeric site; they are placed one at each here,
# which is the common arrangement.
KIR_ORDER: tuple = (
    "KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DL3", "KIR2DL5B", "KIR2DS3",
    "KIR2DP1", "KIR2DL1", "KIR3DP1", "KIR2DL4", "KIR3DL1", "KIR3DS1",
    "KIR2DL5A", "KIR2DS5", "KIR2DS1", "KIR2DS4", "KIR3DL2",
)

LOCUS_PRESETS: dict = {
    "hla-classic": CLASSICAL_HLA_ORDER,
    "hla-all": HLA_ALL_ORDER,
    "kir": KIR_ORDER,
    # "auto" is resolved against the data, not a fixed list
    "auto": None,
}

# Categorical palette, validated for the adjacent pairlist in both modes.
_PALETTE = {
    "light": ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
              "#e87ba4", "#008300", "#6250d6", "#e34948"],
    "dark": ["#3987e5", "#d95926", "#199e70", "#c98500",
             "#d55181", "#008300", "#9085e9", "#e66767"],
}
_THEME = {
    "light": dict(surface="#fcfcfb", primary="#0b0b0b", secondary="#52514e",
                  muted="#8a8880", ribbon="#4a4a46"),
    "dark": dict(surface="#1a1a19", primary="#f0efec", secondary="#c3c2b7",
                 muted="#8a8880", ribbon="#cfcec6"),
}


# --------------------------------------------------------------------- helpers


_KIR_BARE = re.compile(r"^\d+D[LSP]\d+[AB]?$")


def _norm_locus(name: str) -> str:
    """'HLA-DRB1' / 'HLA_DRB1' / 'drb1' -> 'DRB1'; '2DL1' -> 'KIR2DL1'.

    HLA loses its gene-family prefix and KIR gains one, so both end up written
    the way the presets write them whichever convention the caller's data uses.
    """
    s = str(name).strip().upper()
    for prefix in ("HLA-", "HLA_", "HLA."):
        if s.startswith(prefix):
            s = s[len(prefix):]
            break
    for prefix in ("KIR-", "KIR_", "KIR."):
        if s.startswith(prefix):
            s = "KIR" + s[len(prefix):]
            break
    if _KIR_BARE.match(s):
        s = "KIR" + s
    return s


def _resolve_loci(loci, pairs: pd.DataFrame, verbose: bool) -> list:
    """Turn a preset name, or an explicit sequence, into the column order.

    "auto" reads the loci out of `pairs` and orders them by the presets, so a
    panel that mixes classical HLA with DRB3/4/5 or MICA lands in chromosome
    order without anyone having to type the list out. Loci no preset knows go
    after the ones that are known, alphabetically — arbitrary, but stable.
    """
    if loci is None:
        loci = "hla-classic"
    if not isinstance(loci, str):
        return list(loci)

    key = loci.strip().lower().replace("_", "-")
    key = {"classic": "hla-classic", "classical": "hla-classic", "hla": "hla-classic",
           "all": "auto"}.get(key, key)
    if key not in LOCUS_PRESETS:
        raise ValueError(
            f"unknown locus preset {loci!r}; choose from "
            f"{sorted(LOCUS_PRESETS)} or pass an explicit sequence of loci"
        )
    if key != "auto":
        return list(LOCUS_PRESETS[key])

    present = {
        _norm_locus(x)
        for col in ("locus_a", "locus_b")
        for x in pd.unique(pairs[col])
    }
    known = [
        loc for loc in (*HLA_ALL_ORDER, *KIR_ORDER) if loc in present
    ]
    rest = sorted(present - set(known))
    if verbose and rest:
        print(
            f"[info] loci='auto': {len(rest)} locus/loci no preset places "
            f"({', '.join(rest[:6])}{' …' if len(rest) > 6 else ''}) "
            "appended alphabetically; pass an explicit `loci` to order them"
        )
    return known + rest


def _wsym(weight: str) -> str:
    """Pretty symbol for the weight column, for the key."""
    return {"r2": "r²", "r2_dosage": "r²", "r2_freq": "r²", "Dprime": "|D'|", "Wn": "|Wn|"}.get(
        weight, weight
    )


# Node-size modes -> the pair of columns `pair_ld` writes them into.
_FREQ_COLS = {
    "allele_freq": ("allele_freq_a", "allele_freq_b"),
    "carrier_freq": ("carrier_freq_a", "carrier_freq_b"),
}

# Columns that are signed: width is their magnitude. D' and Wn run -1..1, and a
# strong negative D' is strong linkage (repulsion), not a weak link — filtering on
# the raw value would silently drop exactly those pairs.
_SIGNED = {"Dprime", "Wn", "r", "D", "r_dosage", "r_freq", "delta"}


def _trim_alleles(links: pd.DataFrame, max_alleles: int, trim_by: str,
                  verbose: bool) -> pd.DataFrame:
    """Keep only the most important alleles at each locus.

    A long tail of rare alleles is what turns this figure into a comb: each one
    contributes a hairline node, and with D' weighting especially they all clear
    a low threshold, because D' is inflated when one allele is rare.

    Which alleles to keep is the whole question, and ranking them by summed link
    weight — the obvious choice, and what this did originally — is actively wrong
    for a saturating statistic. |D'| hits its ceiling of 1 whenever one of the
    four haplotype combinations is unobserved, which a rare allele all but
    guarantees. Ranking on that sum therefore keeps the RAREST alleles and throws
    the common ones away: an observed D' panel of eight loci retained A*66:01 and
    B*57:01 as its tallest European blocks and dropped A*02:01, B*07:02,
    DRB1*15:01 and DPB1*04:01 entirely, so not one classical haplotype could be
    drawn. Ranking by frequency is immune to that, and matches what the node
    heights already mean.
    """
    by = trim_by
    have_freq = {"src_freq", "tgt_freq"} <= set(links.columns)
    if by == "auto":
        by = "frequency" if have_freq else "weight"
    if by == "frequency" and not have_freq:
        if verbose:
            print(
                "[warn] trim_by='frequency' needs the frequency columns, which these "
                "pairs lack; trimming by summed weight instead"
            )
        by = "weight"

    col_src, col_tgt = (
        ("src_freq", "tgt_freq") if by == "frequency" else ("weight", "weight")
    )
    # frequency is a property of the allele, so take the max over its rows;
    # summed weight is a property of its links, so take the sum
    agg = "max" if by == "frequency" else "sum"

    keep: set = set()
    for loc in pd.unique(
        np.concatenate([links["src_locus"].to_numpy(), links["tgt_locus"].to_numpy()])
    ):
        w = pd.concat(
            [
                links.loc[links["src_locus"] == loc].groupby("src_allele")[col_src].agg(agg),
                links.loc[links["tgt_locus"] == loc].groupby("tgt_allele")[col_tgt].agg(agg),
            ]
        )
        w = w.groupby(level=0).max().sort_values(ascending=False)
        keep.update((loc, a) for a in w.head(int(max_alleles)).index)

    mask = [
        (s, sa) in keep and (t, ta) in keep
        for s, sa, t, ta in zip(
            links["src_locus"], links["src_allele"],
            links["tgt_locus"], links["tgt_allele"],
        )
    ]
    out = links[mask]
    if out.empty:
        raise ValueError(f"no links survive max_alleles={max_alleles}")
    return out


def _short_allele(allele: str) -> str:
    """'B*07:02' -> '07:02'; leaves names without a '*' alone."""
    s = str(allele)
    return s.split("*", 1)[1] if "*" in s else s


def _bezier_ribbon(x0, y0a, y0b, x1, y1a, y1b, curve: float = 0.42) -> Path:
    """Filled ribbon from a slice of one node to a slice of the next."""
    xm0 = x0 + (x1 - x0) * curve
    xm1 = x1 - (x1 - x0) * curve
    verts = [
        (x0, y0a),
        (xm0, y0a), (xm1, y1a), (x1, y1a),      # top edge
        (x1, y1b),                               # down the target face
        (xm1, y1b), (xm0, y0b), (x0, y0b),      # bottom edge back
        (x0, y0a),
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.LINETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.CLOSEPOLY,
    ]
    return Path(verts, codes)


def _place_headers(fig, ax, jobs, fontsize: float, color: str):
    """Locus names above the columns, shrunk then rotated if they would collide.

    The long presets put 17+ columns in the same width and KIR names run to nine
    characters, so a fixed type size either overlaps or has to be set small enough
    to be useless at eight columns. Measure the column pitch and fit to it.
    """
    texts = [
        ax.text(x, y, name, ha="center", va="bottom", fontsize=fontsize,
                fontweight="bold", color=color)
        for x, y, name in jobs
    ]
    if not texts:
        return
    fig.canvas.draw()
    pitch = abs(
        ax.transData.transform((1, 0))[0] - ax.transData.transform((0, 0))[0]
    )
    widest = max(t.get_window_extent().width for t in texts)
    if widest <= pitch * 0.92 or pitch <= 0:
        return

    # Shrink, but only so far — past ~70% the names stop being readable and
    # rotating is the better trade.
    shrink = max(0.70, pitch * 0.92 / widest)
    for t in texts:
        t.set_fontsize(fontsize * shrink)
    fig.canvas.draw()
    widest = max(t.get_window_extent().width for t in texts)
    if widest <= pitch * 0.92:
        return

    for t in texts:
        t.set_rotation(35)
        t.set_rotation_mode("anchor")
        t.set_ha("left")
    fig.canvas.draw()
    # rotated names reach up and right; give the tallest one room
    px_per_y = abs(
        ax.transData.transform((0, 1))[1] - ax.transData.transform((0, 0))[1]
    )
    if px_per_y > 0:
        rise = max(t.get_window_extent().height for t in texts) / px_per_y
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, max(hi, jobs[0][1] + rise * 1.05))


def _place_labels(fig, ax, jobs, fontsize: float, node_width: float, color: str):
    """Label as many nodes as fit, biggest allele first.

    Must run AFTER the limits and tight_layout are settled: how much vertical room
    one line of text takes in DATA units depends on the axes' final size on the
    page, and that is the whole point — ask for a bigger font and the labels thin
    themselves out instead of overprinting each other.
    """
    if not jobs:
        return
    fig.canvas.draw()
    bbox = ax.get_window_extent()
    span = ax.get_ylim()[1] - ax.get_ylim()[0]
    if bbox.height <= 0:
        return
    line_h = fontsize * fig.dpi / 72.0 * 1.22      # px per line, with leading
    min_gap = line_h / bbox.height * span

    taken: dict = {}
    for loc, x, y, h, text in sorted(jobs, key=lambda j: -j[3]):
        ys = taken.setdefault(loc, [])
        if any(abs(y - other) < min_gap for other in ys):
            continue
        ys.append(y)
        ax.text(
            x + node_width / 2 + 0.018, y, text,
            ha="left", va="center", fontsize=fontsize, color=color, zorder=4,
        )


def _arc(x0, y0, x1, y1, height: float) -> Path:
    """Open arc riding above the columns, for a non-adjacent link."""
    top = max(y0, y1) + height
    verts = [(x0, y0), ((x0 + x1) / 2, top), (x1, y1)]
    return Path(verts, [Path.MOVETO, Path.CURVE3, Path.CURVE3])


def _prepare_links(
    pairs: pd.DataFrame,
    population: Optional[str],
    weight: str,
    loci: Sequence[str],
    min_weight: float,
    top_links: Optional[int],
    adjacent_only: bool,
    node_size: str,
    verbose: bool,
) -> pd.DataFrame:
    need = {"allele_a", "allele_b", "locus_a", "locus_b", weight}
    missing = need - set(pairs.columns)
    if missing:
        avail = [
            c for c in pairs.columns
            if c.startswith(("r2", "r_")) or c in ("Dprime", "Wn", "delta")
        ]
        raise ValueError(
            f"pairs is missing column(s): {sorted(missing)}"
            + (f"; weight-like columns present: {avail}" if avail and weight in missing
               else "")
        )

    # Frequency columns are only needed when they drive the node heights, and we
    # degrade to linkage-sized nodes rather than refusing to draw anything.
    # `freq_kind` says what the carried src_freq/tgt_freq actually are, so nothing
    # downstream has to re-derive it from `node_size` and get it wrong.
    freq_kind = node_size if node_size in _FREQ_COLS else None
    fcols = _FREQ_COLS.get(freq_kind)
    if fcols is not None and not set(fcols) <= set(pairs.columns):
        if verbose:
            print(
                f"[warn] node_size={node_size!r} needs columns {list(fcols)}, which "
                "pairs does not have; sizing nodes by linked weight instead"
            )
        fcols, freq_kind = None, None
    if fcols is None:
        # Even when the heights are not frequencies, carry them if they are there:
        # `trim_by="frequency"` needs them, and that is the one defence against a
        # saturating weight picking the rare alleles. They are NOT used for the
        # heights in that case — freq_kind stays None, so `_node_sizes` still
        # falls through to linked weight and the key still says so.
        cand = _FREQ_COLS["allele_freq"]
        if set(cand) <= set(pairs.columns):
            fcols = cand

    df = pairs.copy()

    # ---- one population -----------------------------------------------------
    if "population" in df.columns:
        available = list(pd.unique(df["population"]))
        if population is None:
            if len(available) > 1:
                raise ValueError(
                    f"pairs holds several populations ({available}); name one with "
                    "population=, or pass population='ALL'"
                )
            population = available[0]
        elif population not in available:
            raise ValueError(f"population {population!r} not in pairs; have {available}")
        df = df[df["population"] == population]
    if df.empty:
        raise ValueError("no rows for that population")

    if "same_locus" in df.columns:
        df = df[~df["same_locus"].astype(bool)]

    # ---- restrict to the named loci, in the given order ----------------------
    order = [_norm_locus(x) for x in loci]
    rank = {loc: i for i, loc in enumerate(order)}
    df = df.assign(
        _la=df["locus_a"].map(_norm_locus), _lb=df["locus_b"].map(_norm_locus)
    )
    df_present = set(df["_la"]) | set(df["_lb"])
    df = df[df["_la"].isin(rank) & df["_lb"].isin(rank)]
    if df.empty:
        have = sorted(df_present)
        raise ValueError(
            "no pairs fall between the requested loci. Asked for "
            f"{order[:6]}{' …' if len(order) > 6 else ''}; pairs holds "
            f"{have[:12]}{' …' if len(have) > 12 else ''}. Pick another preset, "
            "pass loci='auto', or check the locus naming."
        )

    # ---- orient every link so the source is the earlier locus ----------------
    flip = df["_la"].map(rank) > df["_lb"].map(rank)
    src_allele = np.where(flip, df["allele_b"], df["allele_a"])
    tgt_allele = np.where(flip, df["allele_a"], df["allele_b"])
    src_locus = np.where(flip, df["_lb"], df["_la"])
    tgt_locus = np.where(flip, df["_la"], df["_lb"])

    w = pd.to_numeric(df[weight], errors="coerce").to_numpy(dtype=float)
    if weight in _SIGNED:
        w = np.abs(w)

    links = pd.DataFrame(
        {
            "src_locus": src_locus,
            "tgt_locus": tgt_locus,
            "src_allele": src_allele,
            "tgt_allele": tgt_allele,
            "weight": w,
        }
    )
    if fcols is not None:
        fa, fb = fcols
        links["src_freq"] = pd.to_numeric(
            np.where(flip, df[fb], df[fa]), errors="coerce"
        )
        links["tgt_freq"] = pd.to_numeric(
            np.where(flip, df[fa], df[fb]), errors="coerce"
        )
    links = links.dropna(subset=["weight"])
    links = links[links["weight"] > 0]
    links["src_rank"] = links["src_locus"].map(rank)
    links["tgt_rank"] = links["tgt_locus"].map(rank)

    if adjacent_only:
        links = links[links["tgt_rank"] - links["src_rank"] == 1]
        if links.empty:
            raise ValueError(
                "no links between neighbouring loci; the loci present may not be "
                "consecutive in `loci`, or try adjacent_only=False"
            )

    links = links[links["weight"] >= float(min_weight)]
    if links.empty:
        raise ValueError(f"no links left at min_weight={min_weight}")

    if top_links:
        links = (
            links.sort_values("weight", ascending=False)
            .groupby(["src_locus", "tgt_locus"], sort=False)
            .head(int(top_links))
        )

    links.attrs["population"] = population
    links.attrs["freq_kind"] = freq_kind
    return links.sort_values("weight", ascending=False).reset_index(drop=True)


def _node_sizes(links: pd.DataFrame, cols: Sequence[str], nodes: dict,
                node_size: str, min_node: float) -> tuple:
    """Height quantity per node, and what that quantity is called.

    With frequencies available a node's height is the allele's frequency in this
    population — the same allele has the same height in every figure drawn from
    the same data, whatever the LD threshold. Falling back to linkage, it is
    max(incoming, outgoing) total weight: the conventional Sankey choice, which
    keeps both faces wide enough for the ribbons without double-counting a middle
    column against an end one.
    """
    have_freq = {"src_freq", "tgt_freq"} <= set(links.columns)
    freq_kind = links.attrs.get("freq_kind")
    if freq_kind is not None and have_freq:
        src = links.groupby(["src_locus", "src_allele"])["src_freq"].max()
        tgt = links.groupby(["tgt_locus", "tgt_allele"])["tgt_freq"].max()
        size = {}
        for loc in cols:
            for a in nodes[loc]:
                v = max(
                    float(src.get((loc, a), 0.0) or 0.0),
                    float(tgt.get((loc, a), 0.0) or 0.0),
                )
                # A frequency of 0 would collapse the block to a hairline and
                # take its ribbons with it; keep a floor so the link stays legible.
                size[(loc, a)] = max(v, float(min_node))
        return size, freq_kind

    out_w = links.groupby(["src_locus", "src_allele"])["weight"].sum()
    in_w = links.groupby(["tgt_locus", "tgt_allele"])["weight"].sum()
    size = {
        (loc, a): float(max(out_w.get((loc, a), 0.0), in_w.get((loc, a), 0.0)))
        for loc in cols
        for a in nodes[loc]
    }
    return size, "linkage"


def _layout(links: pd.DataFrame, loci: Sequence[str], gap_unit: float,
            gap_budget: float, sweeps: int, node_size: str, min_node: float,
            balance: str):
    """Node order and y extents per column.

    Every column shares ONE size -> height scale, so a node of a given frequency
    is the same height anywhere in the figure and the node key means something.
    The ribbon scale is a separate matter, settled by the caller.
    """
    nodes: dict = {}
    for loc in loci:
        present = pd.unique(
            np.concatenate(
                [
                    links.loc[links["src_locus"] == loc, "src_allele"].to_numpy(),
                    links.loc[links["tgt_locus"] == loc, "tgt_allele"].to_numpy(),
                ]
            )
            if ((links["src_locus"] == loc).any() or (links["tgt_locus"] == loc).any())
            else np.array([], dtype=object)
        )
        if len(present):
            nodes[loc] = list(present)
    if not nodes:
        raise ValueError("no nodes to draw")

    cols = [loc for loc in loci if loc in nodes]

    size, size_kind = _node_sizes(links, cols, nodes, node_size, min_node)

    # start each column biggest-first, then sweep
    for loc in cols:
        nodes[loc] = sorted(nodes[loc], key=lambda a: -size[(loc, a)])

    def positions(loc):
        return {a: i for i, a in enumerate(nodes[loc])}

    for _ in range(max(0, int(sweeps))):
        for direction in ("forward", "back"):
            seq = cols[1:] if direction == "forward" else cols[-2::-1]
            for loc in seq:
                if direction == "forward":
                    prev = cols[cols.index(loc) - 1]
                    sub = links[(links["src_locus"] == prev) & (links["tgt_locus"] == loc)]
                    key_col, self_col = "src_allele", "tgt_allele"
                else:
                    prev = cols[cols.index(loc) + 1]
                    sub = links[(links["src_locus"] == loc) & (links["tgt_locus"] == prev)]
                    key_col, self_col = "tgt_allele", "src_allele"
                if sub.empty:
                    continue
                pos = positions(prev)
                bary = {}
                for a in nodes[loc]:
                    s = sub[sub[self_col] == a]
                    if s.empty:
                        bary[a] = np.inf        # unlinked nodes sink to the bottom
                    else:
                        w = s["weight"].to_numpy()
                        p = s[key_col].map(pos).to_numpy(dtype=float)
                        bary[a] = float((p * w).sum() / w.sum())
                nodes[loc] = sorted(nodes[loc], key=lambda a: (bary[a], -size[(loc, a)]))

    # Per-column gap, shrunk so the gaps can never eat more than `gap_budget` of
    # the height. Without this cap a column of 50 alleles spends 60%+ of the axes
    # on whitespace, every other column is squeezed to a fraction of the height,
    # and the figure drifts away from its header.
    gaps = {
        loc: min(gap_unit, gap_budget / max(len(nodes[loc]) - 1, 1))
        for loc in cols
    }

    # Per-column scale that would make that column exactly fill the axes.
    fills: dict = {}
    for loc in cols:
        tot = sum(size[(loc, a)] for a in nodes[loc])
        avail = 1.0 - gaps[loc] * max(len(nodes[loc]) - 1, 0)
        fills[loc] = avail / tot if (tot > 0 and avail > 0) else np.inf

    # One scale for the whole figure: the tallest column fills the axes and every
    # other column is honestly shorter. A locus whose alleles are mostly rare, or
    # mostly filtered out, is then a short column with white space around it —
    # that is the price of heights being comparable between loci.
    scale = min(fills.values())
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError(
            "cannot scale the layout; every node has zero size — check the "
            f"{'frequency' if size_kind != 'linkage' else 'weight'} column"
        )
    # balance="column" spends that comparability to buy the white space back:
    # each column fills the axes on its own scale, so a block is the allele's
    # share WITHIN its locus among the alleles drawn, not its frequency.
    scales = (
        {loc: scale for loc in cols}
        if balance == "global"
        else {loc: (fills[loc] if np.isfinite(fills[loc]) else scale) for loc in cols}
    )

    extents: dict = {}
    for loc in cols:
        heights = [size[(loc, a)] * scales[loc] for a in nodes[loc]]
        drawn = sum(heights) + gaps[loc] * max(len(heights) - 1, 0)
        y = 0.5 + drawn / 2.0                       # centre each column
        for a, h in zip(nodes[loc], heights):
            extents[(loc, a)] = (y - h, y)
            y -= h + gaps[loc]
    return cols, nodes, size, extents, scale, size_kind


# ------------------------------------------------------------------- main plot


def plot_hla_linkage(
    pairs: pd.DataFrame,
    population: Optional[str] = None,
    *,
    weight: str = "r2",
    loci: Union[str, Sequence[str], None] = "hla-classic",
    min_weight: float = 0.05,
    top_links: Optional[int] = None,
    adjacent_only: bool = True,
    node_size: str = "allele_freq",
    min_node: float = 0.004,
    ribbon_fill: float = 0.88,
    max_arcs: Optional[int] = 12,
    link_color: str = "locus",
    label_alleles: bool = True,
    short_labels: bool = True,
    min_label_frac: float = 0.0,
    label_fontsize: float = 11.0,
    locus_fontsize: float = 14.0,
    key_fontsize: float = 10.0,
    title_fontsize: float = 16.0,
    balance: str = "global",
    show_key: bool = True,
    theme: str = "light",
    node_width: float = 0.055,
    gap_unit: float = 0.010,
    gap_budget: float = 0.18,
    max_alleles: Optional[int] = None,
    trim_by: str = "auto",
    sweeps: int = 4,
    figsize: tuple = (13.0, 7.5),
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    verbose: bool = True,
):
    """Draw classical HLA linkage as a Sankey, ribbons weighted by r2.

    Parameters
    ----------
    pairs : DataFrame
        Output of `allele_pair_ld.pair_ld`.
    population : str, optional
        Which population's rows to draw. Required when `pairs` holds more than
        one; pass "ALL" for the all-samples rows.
    weight : str
        Column giving the ribbon width. "r2" by default, the squared correlation
        of presence and absence. "Dprime" gives strength of association rather
        than predictability. Signed columns ("Dprime", "Wn", "r", "D") are taken
        as magnitudes, so a strong negative D' draws as a thick ribbon rather
        than being filtered away.
    node_size : {"allele_freq", "carrier_freq", "linkage"}
        What an allele block's height means. "allele_freq" (default) is the
        allele's frequency in this population, so block height is population
        composition and ribbon width is linkage — two independent readings.
        "carrier_freq" is the share of samples carrying it. "linkage" restores
        the old behaviour, height = total linked weight. Falls back to "linkage"
        with a warning if `pairs` carries no frequency columns.
    min_node : float
        Floor on a node's size quantity, so a very rare allele still draws as a
        visible block with its ribbons attached rather than a hairline.
    ribbon_fill : float
        Share of the tightest node's height the ribbons may fill. Below 1 the
        stacked ribbons keep a little air between them; 1.0 packs them to the
        node edges. Thickness stays proportional to `weight` either way — this
        only sets the one global factor.
    loci : str or sequence
        Which loci to draw, left to right. A preset name:

          "hla-classic"  the eight classical genes (default)
          "hla-all"      every HLA and HLA-region gene, chr6 coordinate order:
                         the non-classical class I (F, G, E, MICA/MICB), the DRB
                         paralogues, DQA2/DQB2, the DO/DM genes and DPA2/DPB2
          "kir"          the KIR cluster at 19q13.4, centromeric to telomeric
          "auto"         whatever loci are in `pairs`, ordered by the presets

        or an explicit sequence, which is drawn in the order given. Names are
        matched after normalising the gene-family prefix, so "HLA-DRB1", "DRB1",
        "KIR2DL1" and "2DL1" all land where you would expect.

        A note on "kir": KIR is on chromosome 19 and HLA on chromosome 6, so a
        KIR-HLA pair has no linkage to draw — their relationship is receptor and
        ligand, not LD. Draw them as separate figures; mixing them into one
        Sankey makes the near-zero cross-chromosome ribbons look meaningful.
    min_weight : float
        Drop links below this r2. The main lever on how busy the figure is.
    top_links : int, optional
        Keep only the N strongest links between each pair of columns.
    adjacent_only : bool
        Ribbons only between neighbouring columns. False also draws long-range
        links as arcs above the columns — readable only for a handful of them,
        so raise `min_weight` or lean on `max_arcs` in that mode.
    max_arcs : int, optional
        With adjacent_only=False, draw only this many of the strongest
        long-range links. None draws them all, which is usually a hairball.
    link_color : {"locus", "gray"}
        Ribbons tinted by source locus, or neutral grey.
    min_label_frac : float
        Hard floor: never label a node shorter than this share of the column. 0 by
        default, because labels are now placed by collision instead — within a
        column they are laid down largest allele first and one is dropped only
        when it would overlap a label already placed, so raising the font size
        thins the labels out rather than turning the column into a smear.
    label_fontsize, locus_fontsize, key_fontsize, title_fontsize : float
        Type sizes in points. `label_fontsize` drives the label-collision spacing,
        so a larger value directly means fewer, more legible allele labels.
    balance : {"global", "column"}
        "global" (default) puts every column on ONE height scale, so a block of a
        given frequency is the same height at every locus and short columns are
        honestly short. "column" rescales each column to fill the axes: it removes
        nearly all the white space, at the cost that heights only compare within a
        locus — a block is then the allele's share among the alleles drawn there.
    show_key : bool
        Draw the scale key under the figure. False reclaims that strip entirely.
    gap_unit : float
        Vertical gap between stacked nodes, in axes units. Capped per column by
        `gap_budget`, so a locus with many alleles gets a smaller gap rather than
        spending the whole column on whitespace.
    gap_budget : float
        Largest share of a column's height the gaps may consume. Lower it to make
        the allele blocks taller; raise it to separate them more.
    max_alleles : int, optional
        Keep only this many alleles per locus. The lever to reach for when a long
        tail of rare alleles combs the figure out.
    trim_by : {"auto", "frequency", "weight"}
        How `max_alleles` picks which alleles to keep. "auto" (default) uses
        frequency when `pairs` carries the frequency columns and falls back to
        summed link weight otherwise. **Do not use "weight" with a saturating
        statistic such as D'** — |D'| = 1 is forced by a rare allele, so ranking
        on it keeps the rare alleles and discards the common ones. "weight" is
        kept for r2 figures, where the sum is a reasonable importance measure,
        and for reproducing older output.
    theme : {"light", "dark"}
        Both are validated against their own surface.

    Returns
    -------
    matplotlib Figure.
    """
    if theme not in _THEME:
        raise ValueError("theme must be 'light' or 'dark'")
    if link_color not in {"locus", "gray"}:
        raise ValueError("link_color must be 'locus' or 'gray'")
    if node_size not in set(_FREQ_COLS) | {"linkage"}:
        raise ValueError(
            "node_size must be 'allele_freq', 'carrier_freq' or 'linkage'"
        )
    if balance not in {"global", "column"}:
        raise ValueError("balance must be 'global' or 'column'")
    if trim_by not in {"auto", "frequency", "weight"}:
        raise ValueError("trim_by must be 'auto', 'frequency' or 'weight'")
    pal = _PALETTE[theme]
    th = _THEME[theme]

    loci = _resolve_loci(loci, pairs, verbose)
    links = _prepare_links(
        pairs, population, weight, loci, min_weight, top_links, adjacent_only,
        node_size, verbose,
    )
    pop_label = links.attrs.get("population")
    order = [_norm_locus(x) for x in loci]
    if max_alleles is not None:
        if int(max_alleles) < 1:
            raise ValueError("max_alleles must be at least 1, or None for no limit")
        before = len(links)
        links = _trim_alleles(links, max_alleles, trim_by, verbose)
        if verbose and len(links) < before:
            print(
                f"[info] max_alleles={max_alleles}: kept {len(links)} of {before} links"
            )

    # A statistic that has piled up on its own ceiling cannot rank anything, and
    # the figure then shows whichever tied links the sort happened to surface.
    # |D'| does this routinely once rare alleles are in play.
    if verbose and len(links):
        at_max = float((links["weight"] >= links["weight"].max() - 1e-9).mean())
        if at_max > 0.25:
            print(
                f"[warn] {at_max:.0%} of links sit at the maximum "
                f"{_wsym(weight)}={links['weight'].max():.3g}: ribbon thickness "
                "cannot rank them, and which ones are drawn is decided by "
                "tie-break order. Raise `min_freq` in pair_ld to drop the rare "
                "alleles that saturate D', or weight the figure by r2."
            )
    cols, nodes, size, extents, scale, size_kind = _layout(
        links, order, gap_unit, gap_budget, sweeps, node_size, min_node, balance
    )

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, facecolor=th["surface"])
    else:
        fig = ax.figure
    ax.set_facecolor(th["surface"])

    xpos = {loc: i for i, loc in enumerate(cols)}
    # keyed to the locus's place in the requested order, not its index among the
    # columns that survived filtering — so raising min_weight never repaints the
    # loci that remain
    color_of = {loc: pal[order.index(loc) % len(pal)] for loc in cols}
    # With the longer presets the palette wraps, and two loci 8 apart in the
    # order can end up neighbours once the absent ones drop out. Nudge the
    # repeat along so no two touching columns share a colour.
    for prev, loc in zip(cols, cols[1:]):
        if color_of[loc] == color_of[prev]:
            i = (pal.index(color_of[loc]) + 1) % len(pal)
            color_of[loc] = pal[i]
    wmax = float(links["weight"].max())

    # ---- ribbon scale, independent of the node scale -------------------------
    # Node heights are now frequencies, so the ribbons meeting a node no longer
    # add up to its height. Pick the largest single factor that still lets every
    # node's ribbons stack inside it: one factor for the whole figure, so a given
    # value is the same thickness everywhere and the key below is honest.
    adjacent = links[links["tgt_rank"] - links["src_rank"] == 1]
    drawable = adjacent[
        [
            (s, sa) in extents and (t, ta) in extents
            for s, sa, t, ta in zip(
                adjacent["src_locus"], adjacent["src_allele"],
                adjacent["tgt_locus"], adjacent["tgt_allele"],
            )
        ]
    ]
    link_scale = scale
    if not drawable.empty:
        face = pd.concat(
            [
                drawable.groupby(["src_locus", "src_allele"])["weight"].sum(),
                drawable.groupby(["tgt_locus", "tgt_allele"])["weight"].sum(),
            ]
        ).groupby(level=[0, 1]).max()
        link_scale = np.inf
        for key, need in face.items():
            lo, hi = extents[key]
            if need > 0:
                link_scale = min(link_scale, (hi - lo) / float(need))
        if not np.isfinite(link_scale) or link_scale <= 0:
            link_scale = scale
        # stacked ribbons otherwise meet the node edge-to-edge and read as one block
        link_scale *= float(np.clip(ribbon_fill, 0.05, 1.0))

    # start each face's stack centred in the node, so the slack splits evenly
    start_out: dict = {}
    start_in: dict = {}
    for grp, start in (
        (["src_locus", "src_allele"], start_out),
        (["tgt_locus", "tgt_allele"], start_in),
    ):
        if drawable.empty:
            break
        for key, tot in drawable.groupby(grp)["weight"].sum().items():
            lo, hi = extents[key]
            start[key] = hi - max((hi - lo) - tot * link_scale, 0.0) / 2.0

    # ---- ribbons, heaviest first so thin ones stay visible on top ------------
    cursor_out: dict = {}
    cursor_in: dict = {}
    for row in adjacent.sort_values("weight", ascending=False).itertuples():
        sloc, tloc = row.src_locus, row.tgt_locus
        if sloc not in xpos or tloc not in xpos:
            continue
        s_key, t_key = (sloc, row.src_allele), (tloc, row.tgt_allele)
        if s_key not in extents or t_key not in extents:
            continue

        # one global scale, so thickness is the LD value itself wherever it is drawn
        s_share = t_share = row.weight * link_scale

        s_top = cursor_out.get(s_key, start_out.get(s_key, extents[s_key][1]))
        t_top = cursor_in.get(t_key, start_in.get(t_key, extents[t_key][1]))
        cursor_out[s_key] = s_top - s_share
        cursor_in[t_key] = t_top - t_share

        frac = row.weight / wmax
        path = _bezier_ribbon(
            xpos[sloc] + node_width / 2, s_top, s_top - s_share,
            xpos[tloc] - node_width / 2, t_top, t_top - t_share,
        )
        ax.add_patch(
            PathPatch(
                path,
                facecolor=color_of[sloc] if link_color == "locus" else th["ribbon"],
                edgecolor="none",
                alpha=float(np.clip(0.22 + 0.5 * frac, 0, 1)),
                zorder=1 + frac,
            )
        )

    # ---- long-range links as arcs -------------------------------------------
    arc_top = 1.0
    if not adjacent_only:
        far = links[links["tgt_rank"] - links["src_rank"] > 1]
        far = far.sort_values("weight", ascending=False)
        if max_arcs is not None and len(far) > int(max_arcs):
            if verbose:
                print(
                    f"[info] {len(far)} long-range links; drawing the {int(max_arcs)} "
                    "strongest (raise max_arcs, or min_weight, to change that)"
                )
            far = far.head(int(max_arcs))
        for row in far.itertuples():
            if row.src_locus not in xpos or row.tgt_locus not in xpos:
                continue
            s_key, t_key = (row.src_locus, row.src_allele), (row.tgt_locus, row.tgt_allele)
            if s_key not in extents or t_key not in extents:
                continue
            frac = row.weight / wmax
            height = 0.10 + 0.05 * (row.tgt_rank - row.src_rank)
            arc_top = max(
                arc_top,
                max(np.mean(extents[s_key]), np.mean(extents[t_key])) + height,
            )
            ax.add_patch(
                PathPatch(
                    _arc(
                        xpos[row.src_locus], np.mean(extents[s_key]),
                        xpos[row.tgt_locus], np.mean(extents[t_key]),
                        height=height,
                    ),
                    facecolor="none",
                    edgecolor=color_of[row.src_locus] if link_color == "locus" else th["ribbon"],
                    linewidth=0.6 + 3.2 * frac,
                    alpha=float(np.clip(0.30 + 0.45 * frac, 0, 1)),
                    zorder=0.5,
                )
            )

    # ---- nodes and labels ----------------------------------------------------
    top_drawn = max(hi for (_, _), (lo, hi) in extents.items())
    header_y = max(top_drawn, arc_top) + 0.028
    label_jobs: list = []
    header_jobs: list = []
    for loc in cols:
        x = xpos[loc]
        for a in nodes[loc]:
            lo, hi = extents[(loc, a)]
            ax.add_patch(
                plt.Rectangle(
                    (x - node_width / 2, lo), node_width, max(hi - lo, 1e-4),
                    facecolor=color_of[loc], edgecolor=th["surface"],
                    linewidth=0.8, zorder=3,
                )
            )
            if label_alleles and (hi - lo) >= min_label_frac:
                # Biggest first, so when labels have to be thinned it is the rare
                # alleles that lose theirs.
                label_jobs.append(
                    (loc, x, (lo + hi) / 2, hi - lo,
                     _short_allele(a) if short_labels else str(a))
                )

        header_jobs.append((x, header_y, loc))

    # ---- weight key: reference ribbons at true thickness ---------------------
    # nice round values spanning what is actually drawn
    # The key strip is dead space in a figure that is mostly ribbons, so the
    # reference values are the small ones: the scales are linear, and a true-height
    # swatch for the LARGEST value would cost a fifth of the figure to say what a
    # 5% swatch says. KEY_CAP is the ceiling on a swatch, in axes units.
    KEY_CAP = 0.06
    # The key is laid out in AXES FRACTION across, so it reads the same whether
    # there are 8 columns or 25, and in DATA units down, so the swatches keep
    # their true height. x = 0 is the left spine.
    keyt = blended_transform_factory(ax.transAxes, ax.transData)
    left = 0.0
    swatch_x = 0.22
    swatch_step = 0.17
    key_bottom = -0.04

    def _key_row(y_from, label, vals, scl, draw):
        """One key row, swatches drawn at true scale. Returns the new bottom."""
        if not vals:
            return y_from
        h_max = max(v * scl for v in vals)
        y = y_from - h_max / 2 - 0.018
        ax.text(
            left, y, label, ha="left", va="center", transform=keyt,
            fontsize=key_fontsize, color=th["secondary"],
        )
        x = swatch_x
        for v in vals:
            draw(x, y, v * scl, v)
            x += swatch_step
        return y - h_max / 2 - 0.022

    if show_key:
        key_vals = [
            v for v in (0.05, 0.1, 0.2, 0.25, 0.4, 0.5, 0.6, 0.75, 0.8, 1.0)
            if min_weight <= v <= wmax and v * link_scale <= KEY_CAP
        ]
        if len(key_vals) > 3:
            key_vals = [key_vals[0], key_vals[len(key_vals) // 2], key_vals[-1]]
        if not key_vals:
            key_vals = [round(min(wmax, KEY_CAP / link_scale), 2)]

        def _ribbon_swatch(x, y, h, v):
            ax.add_patch(
                plt.Rectangle(
                    (x, y - h / 2), 0.035, max(h, 1e-4), transform=keyt,
                    facecolor=th["ribbon"], edgecolor="none",
                    alpha=float(np.clip(0.22 + 0.5 * (v / wmax if wmax else 1.0), 0, 1)),
                    clip_on=False, zorder=3,
                )
            )
            ax.text(
                x + 0.045, y, f"{_wsym(weight)} = {v:.2f}", ha="left", va="center",
                transform=keyt, fontsize=key_fontsize, color=th["secondary"],
            )

        key_bottom = _key_row(
            key_bottom, f"ribbon thickness · {_wsym(weight)} ({weight})",
            key_vals, link_scale, _ribbon_swatch,
        )

        # ---- node-height key -------------------------------------------------
        if size_kind == "linkage":
            node_label = f"block height · total linked {_wsym(weight)}"
            node_vals = []
        elif balance == "column":
            # Heights are column-relative here, so a percentage swatch would be a
            # lie: say what the height means and show no scale.
            node_label = (
                "block height · share of the locus, each column scaled to fill"
            )
            node_vals = []
        else:
            node_label = (
                "block height · allele frequency"
                if size_kind == "allele_freq"
                else "block height · carrier frequency"
            )
            smax = max(size.values())
            node_vals = [
                v for v in (0.01, 0.02, 0.05, 0.1, 0.2)
                if v <= smax and v * scale <= KEY_CAP
            ]
            if len(node_vals) > 3:
                node_vals = [node_vals[0], node_vals[len(node_vals) // 2], node_vals[-1]]
            if not node_vals:
                node_vals = [round(min(smax, KEY_CAP / scale), 3)]

        def _node_swatch(x, y, h, v):
            ax.add_patch(
                plt.Rectangle(
                    (x, y - h / 2), 0.008, max(h, 1e-4), transform=keyt,
                    facecolor=th["muted"], edgecolor="none",
                    clip_on=False, zorder=3,
                )
            )
            ax.text(
                x + 0.018, y, f"{v * 100:g}%", ha="left", va="center",
                transform=keyt, fontsize=key_fontsize, color=th["secondary"],
            )

        if node_vals:
            key_bottom = _key_row(key_bottom, node_label, node_vals, scale, _node_swatch)
        else:
            key_bottom -= 0.030
            ax.text(
                left, key_bottom, node_label, ha="left", va="center",
                transform=keyt, fontsize=key_fontsize, color=th["secondary"],
            )
            key_bottom -= 0.022

    if title is None:
        # from what was actually drawn, not from the preset asked for: filtering
        # can leave a "hla-all" call showing only the classical eight
        if all(c.startswith("KIR") for c in cols):
            title = "KIR linkage"
        elif set(cols) <= set(CLASSICAL_HLA_ORDER):
            title = "Classical HLA linkage"
        elif any(c.startswith("KIR") for c in cols):
            title = "KIR and HLA — note: different chromosomes, no linkage"
        else:
            title = "HLA region linkage"
        if pop_label:
            title += f" — {pop_label}"
    ax.set_title(title, fontsize=title_fontsize, color=th["primary"], loc="left", pad=10)

    # Margins pared to what the ink needs: a little air left of the first column,
    # and room on the right for the last column's labels.
    label_room = 0.12 + 0.016 * max(
        (len(t) for *_, t in label_jobs), default=6
    ) * (label_fontsize / 11.0)
    ax.set_xlim(-0.10, len(cols) - 1 + label_room)
    ax.set_ylim(min(key_bottom, -0.02), header_y + 0.055)
    ax.set_axis_off()
    fig.patch.set_facecolor(th["surface"])
    fig.tight_layout(pad=0.4)

    _place_labels(fig, ax, label_jobs, label_fontsize, node_width, th["primary"])
    _place_headers(fig, ax, header_jobs, locus_fontsize, th["primary"])
    return fig


# ----------------------------------------------------------------------- demo


def make_hla_demo_data(n: int = 500, seed: int = 0):
    """Synthetic classical-HLA dosages built from a handful of real-ish haplotypes.

    Returns (dosage, population) shaped for `pair_ld`. Because the samples are drawn
    from whole haplotypes, neighbouring loci carry strong LD that the plot should
    show as a few dominant ribbons — in particular the long A1-B8-DR3 style
    haplotype running the width of the figure.
    """
    rng = np.random.default_rng(seed)
    loci = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]

    haplos = [
        (["01:01", "07:01", "08:01", "03:01", "05:01", "02:01", "01:03", "01:01"], 0.16),
        (["03:01", "07:02", "07:02", "15:01", "01:02", "06:02", "01:03", "04:01"], 0.15),
        (["02:01", "05:01", "44:02", "04:01", "03:01", "03:02", "01:03", "04:01"], 0.13),
        (["29:02", "16:01", "44:03", "07:01", "02:01", "02:02", "01:03", "04:02"], 0.10),
        (["02:01", "03:04", "40:01", "13:01", "01:03", "06:03", "02:01", "04:01"], 0.10),
        (["24:02", "12:03", "52:01", "15:02", "01:03", "06:01", "02:02", "09:01"], 0.09),
        (["11:01", "04:01", "35:01", "01:01", "01:01", "05:01", "01:03", "02:01"], 0.09),
        (["33:03", "14:03", "58:01", "03:01", "05:01", "02:01", "02:01", "13:01"], 0.09),
        (["68:01", "06:02", "57:01", "07:01", "02:01", "03:03", "01:03", "02:01"], 0.09),
    ]
    probs = np.array([h[1] for h in haplos], dtype=float)
    probs /= probs.sum()

    pops = np.array(["EUR"] * (n // 2) + ["EAS"] * (n - n // 2))
    # the two populations draw the same haplotypes at different frequencies
    p_eur = probs.copy()
    p_eas = probs[::-1].copy()
    p_eas /= p_eas.sum()

    counts: dict = {}
    for i in range(n):
        p = p_eur if pops[i] == "EUR" else p_eas
        for _ in range(2):          # two haplotypes per diploid sample
            h = haplos[rng.choice(len(haplos), p=p)][0]
            for loc, allele in zip(loci, h):
                # occasional recombinant, so the LD is strong but not absolute
                if rng.random() < 0.06:
                    allele = haplos[rng.integers(len(haplos))][0][loci.index(loc)]
                counts.setdefault(f"{loc}*{allele}", np.zeros(n))[i] += 1

    samples = [f"S{i:04d}" for i in range(n)]
    dosage = pd.DataFrame(counts, index=samples).T.sort_index()
    return dosage, pd.Series(pops, index=samples, name="population")
