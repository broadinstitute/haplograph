"""
hla_kir_plots — the two KIR/HLA views that are biologically meaningful.

HLA (chr6) and KIR (chr19) are unlinked, so within-individual r2 between them is
~0 by construction and an LD plot across the two is not informative. What IS
informative is the receptor-ligand relationship: which HLA alleles carry which
KIR ligand epitope, and whether an individual carrying a given KIR gene also
carries its ligand.

    from hla_kir_ligand import build_ligand_table
    from hla_kir_plots import plot_kir_ligand_sankey, plot_ligand_dosage_heatmap

    tbl = build_ligand_table(cache_dir="ipd")
    fig1 = plot_kir_ligand_sankey(dosage_df, tbl, population=pops, group="EUR")
    fig2 = plot_ligand_dosage_heatmap(dosage_df, tbl)

`dosage_df` is the same allele x sample matrix pair_ld takes, holding BOTH the
HLA allele dosages and the KIR gene copy numbers. KIR rows are recognised by a
name starting with "KIR"; HLA rows by an A*/B*/C* prefix.
"""

from __future__ import annotations

import re
from typing import Dict, Iterable, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from hla_kir_ligand import KIR_LIGAND_RECEPTORS, assign_ligands
from hla_linkage_plot import _PALETTE, _THEME, _bezier_ribbon, _short_allele

__all__ = [
    "ligand_dosage",
    "ligand_called",
    "locus_dosage",
    "call_status",
    "call_summary",
    "excess_calls",
    "kir_dosage",
    "kir_ligand_association",
    "plot_kir_ligand_sankey",
    "plot_ligand_dosage_heatmap",
]

LIGAND_ORDER = ["C1", "C2", "Bw4-80I", "Bw4-80T", "A3/A11"]

# single-hue sequential ramp (blue), light -> dark, for the heatmap
_SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#2a78d6", "#256abf",
        "#1c5cab", "#184f95", "#104281", "#0d366b"]


# ------------------------------------------------------------------- reshaping


def _subset(dosage: pd.DataFrame, samples: Optional[Sequence]) -> pd.DataFrame:
    return dosage if samples is None else dosage.loc[:, list(samples)]


# Columns that mark a `pair_ld` output rather than a dosage matrix. Passing the
# pair table here is an easy slip, because both are DataFrames and the failure
# that follows ("0 samples", "no HLA rows") points nowhere near the real cause.
_PAIR_TABLE_COLS = {"src_allele", "tgt_allele", "src_locus", "tgt_locus", "weight"}


def _check_dosage(dosage) -> pd.DataFrame:
    """Fail early, and by name, when this is not an allele x sample matrix."""
    if not isinstance(dosage, pd.DataFrame):
        raise TypeError(
            f"dosage must be an allele x sample DataFrame, got {type(dosage).__name__}"
        )
    hit = _PAIR_TABLE_COLS & set(map(str, dosage.columns))
    if hit:
        raise ValueError(
            f"this looks like a pair_ld() pair table, not a dosage matrix "
            f"(it has columns {sorted(hit)}). Pass the allele x sample dosage "
            "matrix -- the same object you gave pair_ld as its first argument."
        )
    if dosage.empty:
        raise ValueError("dosage matrix is empty")
    return dosage


def _as_sample_id(x) -> str:
    """Normalise a sample identifier for matching.

    Person IDs survive a round trip through a CSV as floats often enough that
    '1000001' and '1000001.0' both turn up, and one side being int and the other
    str is the usual reason a population map silently matches nothing.
    """
    s = str(x).strip()
    return s[:-2] if s.endswith(".0") and s[:-2].isdigit() else s


def _pop_series(population) -> pd.Series:
    """Coerce a population map (Series / dict / one-column DataFrame) to a Series."""
    if isinstance(population, pd.DataFrame):
        if population.shape[1] == 1:
            pop = population.iloc[:, 0]
        else:
            cand = [c for c in population.columns
                    if re.search(r"anc|pop|group|ethn", str(c), re.I)]
            if len(cand) != 1:
                raise ValueError(
                    "population is a DataFrame with columns "
                    f"{list(population.columns)}; pass a single Series indexed by "
                    "sample id, e.g. df.set_index('person_id')['ancestry_pred']"
                )
            pop = population[cand[0]]
    else:
        pop = pd.Series(population)
    pop = pop.copy()
    pop.index = [_as_sample_id(i) for i in pop.index]
    # a duplicated id would make .loc return a Series and break the comparison
    return pop[~pop.index.duplicated(keep="first")]


def _resolve_samples(
    dosage: pd.DataFrame,
    population,
    group: Optional[str],
    *,
    min_samples: int = 3,
) -> Optional[list]:
    """Columns of `dosage` belonging to `group`, or None for "use everything".

    Matching is case- and whitespace-insensitive, because ancestry labels arrive
    as 'afr' from one table and 'AFR' from the next. When nothing matches, the
    error says which of the two possible reasons it was: the ids do not line up,
    or the label does not exist.
    """
    if population is None:
        return None
    if group is None:
        raise ValueError("pass group= to say which population to use")

    pop = _pop_series(population)
    col_ids = {_as_sample_id(c): c for c in dosage.columns}

    shared = [s for s in col_ids if s in pop.index]
    if not shared:
        raise ValueError(
            "the population map and the dosage matrix share no sample ids, so "
            "every population comes out empty.\n"
            f"  dosage columns  (n={dosage.shape[1]}): {list(dosage.columns)[:5]}\n"
            f"  population index (n={len(pop)}): {list(pop.index)[:5]}\n"
            "The usual cause is a population Series that kept its default "
            "RangeIndex -- it must be indexed by sample id "
            "(df.set_index('person_id')['ancestry_pred'])."
        )

    want = str(group).strip().casefold()
    samples = [
        col_ids[s] for s in shared
        if str(pop.loc[s]).strip().casefold() == want
    ]
    if not samples:
        counts = (
            pop.loc[shared].astype(str).str.strip().value_counts().to_dict()
        )
        raise ValueError(
            f"population {group!r} is not a label in the population map.\n"
            f"  labels present among the {len(shared)} shared samples: {counts}"
        )
    if len(samples) < min_samples:
        raise ValueError(
            f"population {group!r} has {len(samples)} samples "
            f"(need at least {min_samples})"
        )
    return samples


def _hla_rows(dosage: pd.DataFrame) -> list:
    return [a for a in dosage.index if re.match(r"^(HLA-)?[ABC]\*", str(a))]


def _locus_of(allele) -> str:
    m = re.match(r"^(?:HLA-)?([ABC])\*", str(allele))
    return m.group(1) if m else ""


# Which class I loci have to be fully typed before "this person has no copies of
# that ligand" is a statement about them rather than about the assay. C1 needs
# HLA-B as well as HLA-C, because B*46/B*73 carry it; Bw4 needs both B and A.
LIGAND_LOCI: Dict[str, tuple] = {
    "C1": ("C", "B"),
    "C2": ("C",),
    "Bw4-80I": ("B", "A"),
    "Bw4-80T": ("B", "A"),
    "A3/A11": ("A",),
}


def locus_dosage(dosage: pd.DataFrame) -> pd.DataFrame:
    """Total called dosage per class I locus per sample (2.0 when fully typed).

    Summed over *every* allele row at the locus, including null and
    ligand-irrelevant ones -- those are typed, they just carry no epitope.
    """
    rows = _hla_rows(dosage)
    if not rows:
        raise ValueError("no HLA-A/B/C allele rows found in the dosage matrix")
    loc = pd.Series({a: _locus_of(a) for a in rows}, name="locus")
    return dosage.loc[rows].groupby(loc).sum()


#: per-sample, per-locus call status. "excess" is not a biological state --
#: HLA-A/B/C have no copy-number variation, so a total above 2 is a typing
#: artifact (one haplotype emitting two calls, usually a duplicated or
#: mis-joined local assembly). It matters in the opposite direction from a
#: no-call: the extra copy pushes the sample out of the 0-ligand cell.
CALL_STATUS = ("none", "partial", "called", "excess")


def call_status(
    dosage: pd.DataFrame, *, ploidy: float = 2.0, tol: float = 0.01
) -> pd.DataFrame:
    """Locus x sample call status, one of CALL_STATUS."""
    tot = locus_dosage(dosage)
    out = pd.DataFrame("partial", index=tot.index, columns=tot.columns)
    out = out.where(tot > tol, "none")
    out = out.where(~(tot - ploidy).abs().le(tol), "called")
    return out.where(~(tot > ploidy + tol), "excess")


def call_summary(
    dosage: pd.DataFrame, *, ploidy: float = 2.0, tol: float = 0.01
) -> pd.DataFrame:
    """Per-locus counts of each call status -- run this before trusting a figure."""
    st = call_status(dosage, ploidy=ploidy, tol=tol)
    out = (
        st.apply(lambda r: r.value_counts(), axis=1)
        .reindex(columns=list(CALL_STATUS))
        .fillna(0)
        .astype(int)
    )
    out["n"] = st.shape[1]
    out["usable_pct"] = (out["called"] / out["n"] * 100).round(2)
    return out


def excess_calls(
    dosage: pd.DataFrame, *, ploidy: float = 2.0, tol: float = 0.01, top: int = 20
) -> pd.DataFrame:
    """Which allele rows turn up in over-called samples, ranked.

    An over-called locus is a bug upstream, not a property of the person, and it
    is rarely spread evenly: usually a handful of allele rows account for most of
    it, which points at the real cause (a 2-field and a 4-field row for the same
    call, a duplicated contig, paralogous GTF records). `enrichment` is how much
    more often an allele appears in over-called samples than in clean ones --
    anything well above 1 is a suspect.
    """
    st = call_status(dosage, ploidy=ploidy, tol=tol)
    rows = []
    for allele in _hla_rows(dosage):
        loc = _locus_of(allele)
        if loc not in st.index:
            continue
        bad = st.loc[loc].eq("excess").to_numpy()
        ok = st.loc[loc].eq("called").to_numpy()
        carried = (dosage.loc[allele] > 0).to_numpy()
        f_bad = carried[bad].mean() if bad.any() else float("nan")
        f_ok = carried[ok].mean() if ok.any() else float("nan")
        rows.append(
            {
                "allele": allele, "locus": loc,
                "n_excess_samples": int((carried & bad).sum()),
                "freq_in_excess": f_bad,
                "freq_in_called": f_ok,
                "enrichment": (f_bad / f_ok) if f_ok else float("nan"),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("n_excess_samples", ascending=False).head(int(top))


def ligand_called(
    dosage: pd.DataFrame, *, ploidy: float = 2.0, tol: float = 0.01
) -> pd.DataFrame:
    """Boolean ligand-group x sample: is a zero here a real zero?

    True when every locus feeding that ligand group carries exactly its diploid
    dosage. Less than that means a haplotype went uncalled, and a zero for the
    group cannot be told from an untyped allele. More than that is a typing
    artifact -- classical class I has no CNV -- and the phantom copy inflates
    ligand dosage, so it is excluded too rather than counted as a carrier.
    """
    st = call_status(dosage, ploidy=ploidy, tol=tol)
    full = st.eq("called")
    out = {}
    for g, loci in LIGAND_LOCI.items():
        need = [l for l in loci if l in full.index]
        if not need:                       # locus absent from the matrix entirely
            out[g] = pd.Series(False, index=dosage.columns)
        else:
            out[g] = full.loc[need].all(axis=0)
    return pd.DataFrame(out).T.reindex(
        [g for g in LIGAND_ORDER if g in LIGAND_LOCI]
    )


def kir_dosage(dosage: pd.DataFrame) -> pd.DataFrame:
    """KIR gene copy number per sample, summed over alleles of the same gene."""
    rows = [a for a in dosage.index if str(a).upper().startswith("KIR")]
    if not rows:
        raise ValueError("no KIR rows found (expected row names starting with 'KIR')")
    gene = pd.Series(
        {a: re.split(r"[*:\s]", str(a))[0].upper() for a in rows}, name="gene"
    )
    return dosage.loc[rows].groupby(gene).sum()


def ligand_dosage(
    dosage: pd.DataFrame,
    ligand_table: pd.DataFrame,
    verbose: bool = True,
    *,
    mask_uncalled: bool = False,
) -> pd.DataFrame:
    """Ligand-group dosage per sample: how many copies of C1, C2, Bw4-80I ... .

    A diploid individual has two class I haplotypes, so C1 + C2 dosage sums to 2
    across HLA-C; Bw4 dosage counts B and A copies together, which is the usual
    convention for KIR3DL1 ligand.

    With `mask_uncalled`, cells whose contributing loci are not fully typed come
    back as NaN instead of 0. Summing dosage cannot tell "carries no C1" from
    "HLA-C half-called", and the two mean opposite things for a missing-ligand
    count, so anything reporting one should ask for the mask.
    """
    hla = _hla_rows(dosage)
    if not hla:
        raise ValueError("no HLA-A/B/C allele rows found in the dosage matrix")
    groups = assign_ligands(hla, ligand_table, verbose=verbose)
    keep = groups[groups.isin(LIGAND_ORDER)]
    if keep.empty:
        raise ValueError("no HLA alleles mapped onto a KIR ligand group")
    out = dosage.loc[keep.index].groupby(keep).sum()
    out = out.reindex([g for g in LIGAND_ORDER if g in out.index])
    if mask_uncalled:
        ok = ligand_called(dosage).reindex(index=out.index, columns=out.columns)
        out = out.astype(float).where(ok.fillna(False).to_numpy())
    return out


# ------------------------------------------------------- receptor-ligand stats


def _fisher(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p for [[a, b], [c, d]]; nan if scipy is absent."""
    try:
        from scipy.stats import fisher_exact
    except ImportError:
        return float("nan")
    return float(fisher_exact([[a, b], [c, d]])[1])


def kir_ligand_association(
    dosage: pd.DataFrame,
    ligand_table: pd.DataFrame,
    population=None,
    *,
    groups: Optional[Sequence[str]] = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """Test every receptor-ligand pair for departure from independence.

    This is the number the Sankey draws but cannot show you: whether carrying a
    KIR gene and carrying its ligand co-occur more or less often than chance.

    KIR (chr19) and HLA (chr6) are unlinked, so within one panmictic population
    the null is strict independence and any real departure is selection on the
    receptor-ligand *pair* -- the coevolution signal. Across a mixed cohort it is
    not: ancestry drives both allele frequencies, so pooling populations
    manufactures association out of structure alone (Wahlund). Pass `population`
    and the result is stratified per group, which is the only form of this test
    worth reporting on an admixed cohort such as All of Us.

    Every count is over samples whose contributing HLA loci are fully typed, so
    `n` is per ligand group rather than per cohort and `n_hla_not_called` says
    how many were held out. A half-called locus reads as ligand-negative, which
    would otherwise land in the same cell as a true missing ligand.

    Returns one row per (population, ligand, KIR gene) with carrier counts, the
    observed and expected joint carrier frequency, their ratio, the odds ratio
    and a two-sided Fisher p.
    """
    dosage = _check_dosage(dosage)

    strata: list = [(None, None)]
    if population is not None:
        pop = _pop_series(population)
        shared = [c for c in dosage.columns if _as_sample_id(c) in pop.index]
        if not shared:
            _resolve_samples(dosage, population, "any")      # raises the good error
        labels = groups if groups is not None else sorted(
            pd.unique(pop.loc[[_as_sample_id(c) for c in shared]].astype(str))
        )
        strata = [(g, _resolve_samples(dosage, population, g)) for g in labels]

    rows = []
    for label, samples in strata:
        d = _subset(dosage, samples)
        n_total = d.shape[1]
        lig = ligand_dosage(d, ligand_table, verbose=verbose)
        called = ligand_called(d).reindex(index=lig.index, columns=lig.columns)
        kir = kir_dosage(d)
        for g in lig.index:
            # samples whose HLA typing can actually support a zero for this group
            ok = called.loc[g].fillna(False).to_numpy().astype(bool)
            if ok.sum() == 0:
                continue
            n = int(ok.sum())
            has_lig = (lig.loc[g] > 0).to_numpy()[ok]
            for gene in KIR_LIGAND_RECEPTORS.get(g, []):
                if gene not in kir.index:
                    continue
                has_kir = (kir.loc[gene] > 0).to_numpy()[ok]
                a = int((has_lig & has_kir).sum())
                b = int((~has_lig & has_kir).sum())
                c = int((has_lig & ~has_kir).sum())
                e = int((~has_lig & ~has_kir).sum())
                exp = float(has_lig.mean() * has_kir.mean())
                rows.append(
                    {
                        "population": label if label is not None else "ALL",
                        "ligand": g,
                        "kir_gene": gene,
                        "n": n,                      # HLA fully typed for this group
                        "n_total": n_total,
                        "n_hla_not_called": n_total - n,
                        "ligand_carriers": int(has_lig.sum()),
                        "kir_carriers": int(has_kir.sum()),
                        "both": a,
                        "kir_no_ligand": b,          # the missing-ligand genotype
                        "observed_freq": a / n if n else float("nan"),
                        "expected_freq": exp,
                        "obs_over_exp": (a / n / exp) if (n and exp > 0) else float("nan"),
                        "odds_ratio": ((a * e) / (b * c)) if (b and c) else float("nan"),
                        "fisher_p": _fisher(a, b, c, e),
                    }
                )
    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("no KIR gene in the matrix matches a ligand group present")
    return out


# ---------------------------------------------------------------- 1. Sankey


def plot_kir_ligand_sankey(
    dosage: pd.DataFrame,
    ligand_table: pd.DataFrame,
    population=None,
    group: Optional[str] = None,
    *,
    max_alleles: int = 12,
    min_carrier_freq: float = 0.01,
    theme: str = "light",
    node_width: float = 0.05,
    gap_budget: float = 0.30,
    figsize: tuple = (12.5, 8.0),
    title: Optional[str] = None,
    verbose: bool = True,
):
    """HLA allele -> KIR ligand group -> KIR gene.

    Left ribbons are each allele's carrier frequency. Right ribbons are the
    fraction of individuals carrying BOTH that ligand group and that KIR gene --
    the functional receptor-ligand pair frequency, which is what disease-association
    work keys on. Because the two loci are unlinked, that right-hand number is
    close to the product of the two marginals; a visible departure from it is the
    population-level signal Single et al. looked for.
    """
    if theme not in _THEME:
        raise ValueError("theme must be 'light' or 'dark'")
    th, pal = _THEME[theme], _PALETTE[theme]

    dosage = _check_dosage(dosage)
    samples = _resolve_samples(dosage, population, group)
    d = _subset(dosage, samples)
    n = d.shape[1]

    lig = ligand_dosage(d, ligand_table, verbose=verbose)
    kir = kir_dosage(d)

    # --- allele -> ligand, weighted by carrier frequency ---------------------
    # Every frequency here is over samples actually typed at the locus in
    # question, never over the whole cohort: a half-called HLA-C makes a sample
    # look C1-negative, and dividing by n would push every frequency down by
    # however much of the cohort failed to type.
    hla = _hla_rows(d)
    groups = assign_ligands(hla, ligand_table, verbose=False)
    full = locus_dosage(d) >= 2.0 - 0.01
    called = ligand_called(d).reindex(index=lig.index, columns=lig.columns)
    left = []
    for allele, g in groups.items():
        if g not in LIGAND_ORDER:
            continue
        loc = _locus_of(allele)
        ok = (full.loc[loc].to_numpy() if loc in full.index
              else np.ones(n, dtype=bool))
        if ok.sum() == 0:
            continue
        cf = float(((d.loc[allele] > 0).to_numpy() & ok).sum()) / int(ok.sum())
        if cf >= min_carrier_freq:
            left.append({"src": allele, "tgt": g, "w": cf})
    left = pd.DataFrame(left)
    if left.empty:
        raise ValueError(f"no alleles above min_carrier_freq={min_carrier_freq}")
    left = (
        left.sort_values("w", ascending=False)
        .groupby("tgt", sort=False)
        .head(int(max_alleles))
    )

    # --- ligand -> KIR, weighted by joint carrier frequency ------------------
    right = []
    for g in lig.index:
        ok = called.loc[g].fillna(False).to_numpy().astype(bool)
        if ok.sum() == 0:
            continue
        nk = int(ok.sum())
        has_lig = (lig.loc[g] > 0).to_numpy()[ok]
        for gene in KIR_LIGAND_RECEPTORS.get(g, []):
            if gene not in kir.index:
                continue
            has_kir = (kir.loc[gene] > 0).to_numpy()[ok]
            joint = float((has_lig & has_kir).sum()) / nk
            exp = float(has_lig.mean() * has_kir.mean())
            if joint > 0:
                right.append(
                    {
                        "src": g, "tgt": gene, "w": joint,
                        "n_called": nk, "n_total": n,
                        "expected": exp,
                        # >1 means receptor and ligand co-occur more often than
                        # two unlinked loci should; this is the whole point of
                        # the right-hand column, so carry it on the figure.
                        "obs_over_exp": joint / exp if exp > 0 else float("nan"),
                    }
                )
    right = pd.DataFrame(right)
    if right.empty:
        raise ValueError("no KIR gene in the matrix matches a ligand group present")

    # --- layout: three columns ------------------------------------------------
    cols = {
        0: list(dict.fromkeys(left["src"])),
        1: [g for g in LIGAND_ORDER if g in set(left["tgt"]) | set(right["src"])],
        2: list(dict.fromkeys(right["tgt"])),
    }
    size: Dict[tuple, float] = {}
    for a in cols[0]:
        size[(0, a)] = float(left.loc[left["src"] == a, "w"].sum())
    for g in cols[1]:
        size[(1, g)] = max(
            float(left.loc[left["tgt"] == g, "w"].sum()),
            float(right.loc[right["src"] == g, "w"].sum()),
        )
    for k in cols[2]:
        size[(2, k)] = float(right.loc[right["tgt"] == k, "w"].sum())

    # order column 0 by ligand group then weight, so ribbons do not cross
    rank = {g: i for i, g in enumerate(cols[1])}
    cols[0].sort(key=lambda a: (rank.get(left.loc[left["src"] == a, "tgt"].iloc[0], 9),
                                -size[(0, a)]))
    cols[2].sort(key=lambda k: (min(rank.get(g, 9)
                                    for g in right.loc[right["tgt"] == k, "src"]),
                                -size[(2, k)]))

    scale = np.inf
    for c, members in cols.items():
        tot = sum(size[(c, m)] for m in members)
        gap = min(0.03, gap_budget / max(len(members) - 1, 1))
        avail = 1.0 - gap * max(len(members) - 1, 0)
        if tot > 0:
            scale = min(scale, avail / tot)

    extents: Dict[tuple, tuple] = {}
    for c, members in cols.items():
        gap = min(0.03, gap_budget / max(len(members) - 1, 1))
        hs = [size[(c, m)] * scale for m in members]
        y = 0.5 + (sum(hs) + gap * max(len(hs) - 1, 0)) / 2
        for m, h in zip(members, hs):
            extents[(c, m)] = (y - h, y)
            y -= h + gap

    # --- draw ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=figsize, facecolor=th["surface"])
    ax.set_facecolor(th["surface"])
    lig_color = {g: pal[LIGAND_ORDER.index(g) % len(pal)] for g in cols[1]}

    cur_out: Dict[tuple, float] = {}
    cur_in: Dict[tuple, float] = {}
    for df, (cs, ct) in ((left, (0, 1)), (right, (1, 2))):
        for row in df.sort_values("w", ascending=False).itertuples():
            sk, tk = (cs, row.src), (ct, row.tgt)
            if sk not in extents or tk not in extents:
                continue
            h = row.w * scale
            st = cur_out.get(sk, extents[sk][1])
            tt = cur_in.get(tk, extents[tk][1])
            cur_out[sk], cur_in[tk] = st - h, tt - h
            g = row.tgt if cs == 0 else row.src
            ax.add_patch(
                plt.matplotlib.patches.PathPatch(
                    _bezier_ribbon(cs + node_width / 2, st, st - h,
                                   ct - node_width / 2, tt, tt - h),
                    facecolor=lig_color.get(g, th["ribbon"]), edgecolor="none",
                    alpha=float(np.clip(0.25 + 0.45 * row.w / df["w"].max(), 0, 1)),
                    zorder=1,
                )
            )

    for c, members in cols.items():
        for m in members:
            lo, hi = extents[(c, m)]
            col = (lig_color.get(m) if c == 1
                   else lig_color.get(
                       left.loc[left["src"] == m, "tgt"].iloc[0] if c == 0
                       else right.loc[right["tgt"] == m, "src"].iloc[0], th["ribbon"]))
            ax.add_patch(
                plt.Rectangle((c - node_width / 2, lo), node_width, max(hi - lo, 1e-4),
                              facecolor=col, edgecolor=th["surface"], lw=0.8, zorder=3)
            )
            label = str(m)   # keep the locus: "07:01" alone is ambiguous
            ha, dx = ("right", -node_width / 2 - 0.015) if c == 0 else \
                     ("left", node_width / 2 + 0.015)
            ax.text(c + dx, (lo + hi) / 2, label, ha=ha, va="center",
                    fontsize=8.5, color=th["primary"], zorder=4)

    for c, name in enumerate(["HLA allele", "KIR ligand", "KIR gene"]):
        ax.text(c, 1.05, name, ha="center", va="bottom", fontsize=12,
                fontweight="bold", color=th["primary"])

    if title is None:
        title = "KIR receptor–ligand structure"
        if group:
            title += f" — {group}"
        title += f"  (n = {n})"
    ax.set_title(title, fontsize=15, color=th["primary"], loc="left", pad=24)
    oe = right["obs_over_exp"].dropna()
    ax.text(0.0, -0.075,
            "left: allele carrier frequency   ·   right: fraction carrying both "
            "receptor and ligand\n"
            + (f"observed / expected under independence: "
               f"{oe.min():.2f}–{oe.max():.2f}   ·   "
               "full table in fig.pair_table; test with kir_ligand_association()"
               if len(oe) else ""),
            transform=ax.transAxes, fontsize=8.5, color=th["secondary"])
    ax.set_xlim(-0.75, 2.75)
    ax.set_ylim(-0.1, 1.14)
    ax.set_axis_off()
    fig.patch.set_facecolor(th["surface"])
    fig.tight_layout()
    fig.pair_table = right            # attach the numbers behind the right column
    return fig


# --------------------------------------------------------------- 2. heatmap


def plot_ligand_dosage_heatmap(
    dosage: pd.DataFrame,
    ligand_table: pd.DataFrame,
    population=None,
    group: Optional[str] = None,
    *,
    as_percent: bool = True,
    theme: str = "light",
    figsize: tuple = (9.5, 5.5),
    title: Optional[str] = None,
    verbose: bool = True,
):
    """KIR gene copy number against cognate-ligand dosage, per individual.

    One row per KIR gene that has a known ligand, one column per ligand dosage
    (0, 1, 2+ copies). The leftmost column is the missing-ligand genotype: the
    individual carries the receptor but nothing for it to bind. Cells are counts,
    or row percentages when as_percent.

    Receptor carriers whose contributing HLA loci are not fully typed are held
    out of the row entirely and reported separately on the right, because a zero
    ligand dosage from a half-called locus is a no-call, not a missing ligand --
    and the missing-ligand column is the one number this figure exists to show.
    Percentages are therefore over called carriers, not over all carriers.
    """
    th = _THEME[theme]
    dosage = _check_dosage(dosage)
    samples = _resolve_samples(dosage, population, group)
    d = _subset(dosage, samples)
    n = d.shape[1]

    lig = ligand_dosage(d, ligand_table, verbose=verbose)
    called = ligand_called(d).reindex(index=lig.index, columns=lig.columns)
    kir = kir_dosage(d)

    rows, labels, uncalled = [], [], []
    for g in lig.index:
        ok = called.loc[g].fillna(False).to_numpy().astype(bool)
        for gene in KIR_LIGAND_RECEPTORS.get(g, []):
            if gene not in kir.index:
                continue
            carriers = (kir.loc[gene] > 0).to_numpy()
            if carriers.sum() == 0:
                continue
            keep = carriers & ok
            if keep.sum() == 0:
                continue
            ld = lig.loc[g].to_numpy()[keep]
            counts = [int((ld == 0).sum()), int((ld == 1).sum()), int((ld >= 2).sum())]
            rows.append(counts)
            uncalled.append(int((carriers & ~ok).sum()))
            labels.append(f"{gene}  ·  {g}")
    if not rows:
        raise ValueError(
            "no KIR gene with a known ligand and a fully typed HLA locus found "
            "in the matrix"
        )

    M = np.array(rows, dtype=float)
    U = np.array(uncalled, dtype=float)
    shown = M / M.sum(axis=1, keepdims=True) * 100 if as_percent else M

    fig, ax = plt.subplots(figsize=figsize, facecolor=th["surface"])
    ax.set_facecolor(th["surface"])
    cmap = plt.matplotlib.colors.LinearSegmentedColormap.from_list("seq", _SEQ)
    im = ax.imshow(shown, cmap=cmap, aspect="auto",
                   vmin=0, vmax=shown.max() if shown.size else 1)

    ax.set_xticks(range(3))
    ax.set_xticklabels(["0 copies\nmissing ligand", "1 copy", "2+ copies"],
                       fontsize=9.5, color=th["primary"])
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=9.5, color=th["primary"])
    ax.set_xlabel("cognate ligand dosage in fully typed receptor carriers",
                  fontsize=10, color=th["secondary"], labelpad=8)

    # Held-out carriers, as a text column outside the ramp: they belong in the
    # figure (the reader needs to know how much of the cohort the row rests on)
    # but not in the colour scale, which is over called carriers only.
    if U.any():
        ax.set_xlim(-0.6, 3.75)
        ax.axvline(2.62, color=th["muted"], lw=0.8, alpha=0.5, zorder=4)
        ax.text(3.12, -0.72, "HLA not\nfully typed", ha="center", va="bottom",
                fontsize=8.5, color=th["secondary"])
        for i, u in enumerate(U):
            frac = u / (M[i].sum() + u) * 100
            ax.text(3.12, i, f"{int(u)}\n({frac:.0f}%)", ha="center", va="center",
                    fontsize=8.5,
                    color=th["muted"] if frac < 10 else "#e34948")

    # value labels, so the chart never rests on colour alone
    lim = shown.max() * 0.55 if shown.size else 1
    for i in range(shown.shape[0]):
        for j in range(shown.shape[1]):
            txt = f"{shown[i, j]:.0f}%" if as_percent else f"{int(M[i, j])}"
            ax.text(j, i, f"{txt}\n({int(M[i, j])})", ha="center", va="center",
                    fontsize=8.5,
                    color="#ffffff" if shown[i, j] > lim else th["primary"])

    # mark the missing-ligand column; the tick label carries the wording, so
    # nothing is drawn over the axis
    ax.add_patch(plt.Rectangle((-0.5, -0.5), 1, len(labels), fill=False,
                               edgecolor="#e34948", lw=2.0, zorder=5))
    ax.get_xticklabels()[0].set_color("#e34948")
    ax.get_xticklabels()[0].set_fontweight("bold")

    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    if title is None:
        title = "KIR receptor carriers by cognate ligand dosage"
        if group:
            title += f" — {group}"
        title += f"  (n = {n})"
    # the held-out-carrier header sits in the same band as the title
    ax.set_title(title, fontsize=14, color=th["primary"], loc="left",
                 pad=26 if U.any() else 14)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("% of called receptor carriers" if as_percent else "individuals",
                 fontsize=9, color=th["secondary"])
    cb.outline.set_visible(False)
    fig.patch.set_facecolor(th["surface"])
    fig.tight_layout()
    fig.counts = pd.DataFrame(
        np.column_stack([M, U]), index=labels,
        columns=["ligand_0", "ligand_1", "ligand_2plus", "hla_not_called"],
    )
    return fig
