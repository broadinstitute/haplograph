"""
allele_pair_ld — pairwise linkage of an allele x sample dosage matrix, on
presence and absence.

Notebook use
------------
    from allele_pair_ld import pair_ld

    pairs = pair_ld(dosage_df)                      # all samples together
    pairs = pair_ld(dosage_df, population=pops)     # per population, plus ALL
    pairs = pair_ld(dosage_df, population=pops, include_all=False)   # per population only

`dosage_df` is a DataFrame with alleles in the INDEX and samples in the COLUMNS
(pass `samples_in="rows"` if yours is the other way round). Values are dosages.
The allele names must really be in the index: read a saved matrix back with
`pd.read_csv(path, index_col=0)`, not a bare `read_csv`, or the index is row
numbers and every allele name in the output would be one. An allele-name column
left in place (`allele`, `marker`, `index`, …) is promoted with a warning; an
index of bare row numbers is an error rather than a table full of integers.
Missing calls remain unknown and are excluded pairwise. n_samples is the
number observed for both alleles; n_group_samples is the full group size.
Before any linkage statistic, each observed entry becomes presence (dosage > 0) or
absence (dosage == 0). A heterozygote and a homozygote are both carriers, and a
KIR allele at copy number 3 is the same as one at copy number 1. Zygosity and
copy number stay visible only in `allele_freq`, which is still counted from the
original dosages.

These are carrier associations, not haplotype LD or dosage correlations.
Estimating haplotype LD requires phased calls or a separate haplotype model.
Allele frequencies use complete locus totals; partially missing locus calls
are excluded from that dosage-frequency denominator.

Scope
-----
Every statistic is computed independently within one set of samples, and the
`population` column says which set. With no `population` argument there is one
set, labelled ALL. Given one, you get a row per allele pair per population, and
by default an ALL row as well. Nothing is compared across populations and no
statistic mixes them: to see whether a pair behaves the same way everywhere,
read its rows side by side.

Columns
-------
  carrier_r, carrier_r2 (aliases: r, r2)
                        Pearson correlation of the two presence/absence
                        indicators, and its square. For a 0/1 coding this r is
                        the phi coefficient, so it matches `Wn`.
  carrier_D, carrier_dmax, carrier_Dprime (aliases: D, dmax, Dprime)
                        Binary association coefficient D = p11 - p_a*p_b on the carrier
                        table, the feasible bound from the pair's carrier
                        frequencies, and D/Dmax. Because the bound uses the same
                        margins as D, |D'| stays within 1 apart from rounding.
  Wn                    signed Cramer's V (= phi) on the same 2x2 table, the
                        statistic used in Single et al. 2007 (Nat Genet 39:1114).
                        Equal to `r`. Undefined (NaN) where either allele's
                        carrier status is monomorphic in the group.
  n11/n10/n01/n00       raw carrier counts. `carrier_freq` is the share of
                        samples with the allele present (the margin of D and r).
                        `allele_freq` is the dosage frequency, copies of the
                        allele divided by copies at its locus, so homozygotes
                        count twice.
  p_chi2, p_exact       carrier-table test; p_exact (Fisher) is filled for sparse
                        tables when exact=True. q_value is BH-adjusted within
                        each population.

Same-locus pairs are dropped by default. Which alleles a sample carries at one
locus are not independent, so that correlation is not linkage between loci.
"""

from __future__ import annotations

import re
from typing import Iterable, Mapping, Optional, Union

import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2_dist
from scipy.stats import fisher_exact

__all__ = ["pair_ld", "make_demo_data"]

ALL_LABEL = "ALL"


# ------------------------------------------------------------------ coercion

# Column names that hold allele labels when they were never set as the index —
# what `reset_index()` or `read_csv` without `index_col=0` leaves behind.
_LABEL_COLUMNS = (
    "allele", "alleles", "marker", "markers", "variant", "id", "name",
    "index", "unnamed: 0", "allele_id", "marker_id",
)


def _is_positional(index: pd.Index) -> bool:
    """True when the index is row numbers rather than allele names."""
    if isinstance(index, pd.RangeIndex) or pd.api.types.is_integer_dtype(index):
        return True
    labels = pd.Index(index).astype(str)
    return len(labels) > 0 and bool(labels.str.fullmatch(r"\d+").all())


def _recover_allele_labels(df: pd.DataFrame) -> pd.DataFrame:
    """Promote an allele-name column to the index, or say why we can't.

    A matrix whose index is 0, 1, 2 … carries no allele names, and every
    downstream column would quietly report row numbers instead — including
    `locus_a`/`locus_b`, which would then make every pair look cross-locus.
    """
    labelish = [c for c in df.columns if not pd.api.types.is_numeric_dtype(df[c])]
    named = [c for c in labelish if str(c).strip().lower() in _LABEL_COLUMNS]
    pick = named[0] if named else (labelish[0] if len(labelish) == 1 else None)
    if pick is None:
        raise ValueError(
            "alleles are not named: the dosage matrix index is positional "
            "(0, 1, 2, …), so allele_a/allele_b and locus_a/locus_b would carry "
            "row numbers. Put the allele names in the index — "
            "df.set_index('allele'), or pd.read_csv(path, index_col=0) — or pass "
            "samples_in='rows' if it is the samples that are in the index."
            + (f" Non-numeric columns seen: {labelish[:5]}" if labelish else "")
        )
    print(f"[warn] alleles taken from column {pick!r}; the index was row numbers")
    return df.set_index(pick)


def _as_dosage(dosage: pd.DataFrame, samples_in: str) -> pd.DataFrame:
    if not isinstance(dosage, pd.DataFrame):
        raise TypeError("dosage must be a pandas DataFrame")
    if samples_in not in {"columns", "rows"}:
        raise ValueError("samples_in must be 'columns' or 'rows'")
    df = dosage.T.copy() if samples_in == "rows" else dosage.copy()

    if _is_positional(df.index):
        df = _recover_allele_labels(df)

    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    if df.index.has_duplicates:
        raise ValueError(
            f"duplicate allele names: {df.index[df.index.duplicated()].unique()[:5].tolist()}"
        )
    if df.columns.has_duplicates:
        raise ValueError(
            f"duplicate sample names: {df.columns[df.columns.duplicated()].unique()[:5].tolist()}"
        )

    non_numeric = df.select_dtypes(exclude=[np.number]).columns.tolist()
    if non_numeric:
        hint = ""
        if len(non_numeric) == df.shape[1]:
            # Every column non-numeric usually means one stray row of labels.
            hint = (
                " — every column is non-numeric, so the matrix probably still "
                "carries a label row; drop it before calling pair_ld"
            )
        raise TypeError(f"non-numeric sample columns: {non_numeric[:5]}{hint}")
    df = df.astype(np.float64)

    if np.isinf(df.values).any():
        raise ValueError("infinite dosages present")
    if (df.values < 0).any():
        raise ValueError("negative dosages present")
    return df


def _as_population(
    population: Union[pd.Series, pd.DataFrame, Mapping, Iterable, None],
    samples: np.ndarray,
) -> Optional[pd.Series]:
    """Normalise anything sensible into a Series of labels indexed by sample."""
    if population is None:
        return None

    if isinstance(population, pd.DataFrame):
        if population.shape[1] == 1:
            s = population.iloc[:, 0]
        elif population.shape[1] >= 2:
            cols = {str(c).lower(): c for c in population.columns}
            key = next((cols[k] for k in ("sample", "sample_id", "id", "iid") if k in cols), None)
            val = next(
                (cols[k] for k in ("population", "pop", "stratum", "strata", "group", "ancestry")
                 if k in cols),
                None,
            )
            if key is not None and val is not None:
                s = pd.Series(population[val].values, index=population[key].astype(str).values)
            else:
                s = pd.Series(
                    population.iloc[:, 1].values,
                    index=population.iloc[:, 0].astype(str).values,
                )
        else:
            raise ValueError("population DataFrame has no columns")
    elif isinstance(population, pd.Series):
        s = population
    elif isinstance(population, Mapping):
        s = pd.Series(population)
    else:
        arr = np.asarray(list(population), dtype=object)
        if arr.size != samples.size:
            raise ValueError(
                f"population has {arr.size} entries but there are {samples.size} samples"
            )
        s = pd.Series(arr, index=samples)

    s.index = s.index.astype(str)
    s = s[~s.index.duplicated()]
    missing = [x for x in samples if x not in s.index]
    if missing:
        raise ValueError(
            f"{len(missing)} samples have no population label, e.g. {missing[:5]}"
        )
    s = s.reindex(samples)
    # A null label is as missing as an absent one: left alone it matches no
    # group, so the sample would quietly vanish from every per-population block
    # while still counting towards ALL.
    blank = s.isna() | s.astype(str).str.strip().isin({"", "nan", "none", "NA", "N/A"})
    if blank.any():
        examples = s.index[blank][:5].tolist()
        raise ValueError(
            f"{int(blank.sum())} samples have a blank population label, "
            f"e.g. {examples}. Drop them from `dosage`, or label them."
        )
    s = s.astype(str)
    if (s == ALL_LABEL).any():
        raise ValueError(f"{ALL_LABEL!r} is reserved; rename that population label")
    return s


def _as_loci(
    locus: Union[pd.Series, Mapping, None], alleles: np.ndarray, locus_regex: str
) -> np.ndarray:
    if locus is not None:
        m = locus.to_dict() if isinstance(locus, pd.Series) else dict(locus)
        m = {str(k): str(v) for k, v in m.items()}
        missing = [a for a in alleles if a not in m]
        if missing:
            raise ValueError(
                f"{len(missing)} alleles absent from the locus map, e.g. {missing[:5]}"
            )
        return np.array([m[a] for a in alleles], dtype=object)

    pat = re.compile(locus_regex)
    out = []
    for a in alleles:
        hit = pat.match(a)
        out.append(hit.group(1) if hit and hit.lastindex else a)
    return np.array(out, dtype=object)


# ----------------------------------------------------------------- statistics


def _locus_totals(X: np.ndarray, codes: np.ndarray, n_loci: int) -> np.ndarray:
    """Per-sample total dosage at each allele's own locus -> D x n."""
    agg = np.zeros((n_loci, X.shape[1]), dtype=np.float64)
    np.add.at(agg, codes, X)
    return agg[codes]


def _corr(M: np.ndarray) -> np.ndarray:
    """Pearson correlation between rows, NaN where a row has no variance."""
    Mc = M - M.mean(axis=1, keepdims=True)
    SP = Mc @ Mc.T
    SS = (Mc**2).sum(axis=1)
    denom = np.sqrt(np.outer(SS, SS))
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(denom > 0, SP / denom, np.nan)
    return np.clip(r, -1.0, 1.0), SP / M.shape[1]


def _carrier_stats(C: np.ndarray, n: int):
    """Signed Cramer's V (= phi), its chi-square P, and the 2x2 counts."""
    n11 = C @ C.T
    carriers = C.sum(axis=1)
    n10 = carriers[:, None] - n11
    n01 = carriers[None, :] - n11
    n00 = n - n11 - n10 - n01

    marg = (
        carriers[:, None]
        * (n - carriers)[:, None]
        * carriers[None, :]
        * (n - carriers)[None, :]
    )
    num = n11 * n00 - n10 * n01
    with np.errstate(divide="ignore", invalid="ignore"):
        phi = np.where(marg > 0, num / np.sqrt(marg), np.nan)
    phi = np.clip(phi, -1.0, 1.0)
    p = np.where(np.isfinite(phi), chi2_dist.sf(n * np.nan_to_num(phi) ** 2, df=1), np.nan)
    return phi, p, n11, n10, n01, n00, carriers


def _dprime(delta: np.ndarray, p: np.ndarray):
    """Normalise binary carrier covariance by its feasible table bound."""
    pa, pb = p[:, None], p[None, :]
    pos = np.minimum(pa * (1 - pb), pb * (1 - pa))
    neg = np.maximum(-pa * pb, -(1 - pa) * (1 - pb))
    dmax = np.abs(np.where(delta >= 0, pos, neg))
    with np.errstate(divide="ignore", invalid="ignore"):
        dp = np.where(dmax > 0, delta / dmax, np.nan)
    return dp, dmax


def _bh_fdr(p: np.ndarray) -> np.ndarray:
    q = np.full(p.shape, np.nan, dtype=np.float64)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    m = ranked.size
    adj = np.minimum.accumulate((ranked * m / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


# ----------------------------------------------------------------- main entry


def pair_ld(
    dosage: pd.DataFrame,
    population: Union[pd.Series, pd.DataFrame, Mapping, Iterable, None] = None,
    locus: Union[pd.Series, Mapping, None] = None,
    *,
    samples_in: str = "columns",
    locus_regex: str = r"^([^*]+)\*",
    min_carriers: int = 5,
    min_freq: float = 0.0,
    include_same_locus: bool = False,
    include_all: bool = True,
    exact: bool = False,
    top: Optional[int] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """Presence/absence r2, D' and Wn for every allele pair.

    Missing calls are excluded pairwise. n_samples reports the complete-case
    count for each pair; n_group_samples reports group size. Explicit carrier_*
    columns and their legacy aliases describe carrier association, not haplotype LD.

    Dosages are turned into a presence/absence matrix first: an entry greater
    than 0 is a carrier. r, D' and Wn are then all computed on that matrix.

    Parameters
    ----------
    dosage : DataFrame
        Alleles in the index, samples in the columns (see `samples_in`).
    population : Series / DataFrame / dict / array, optional
        Sample -> population. Statistics are then computed separately within each
        population. Without it there is a single ALL group.
    locus : Series / dict, optional
        Allele -> locus. Without it the locus is parsed from the allele name with
        `locus_regex`: the text before the first '*' ("B*07:02" -> "B").
    min_carriers, min_freq : int, float
        Drop alleles carried by fewer than `min_carriers` samples, or below
        `min_freq` dosage allele frequency (copies of the allele / copies at the
        locus). Applied once over all samples, so the same allele set is used in
        every population and the rows stay comparable. The linkage coefficients
        themselves use carrier frequency, not this dosage frequency.
    include_same_locus : bool
        Keep within-locus pairs. Not recommended: the alleles one sample carries
        at a locus are not independent.
    include_all : bool
        With `population` given, also emit rows computed over all samples
        together, labelled ALL.
    exact : bool
        Fisher exact P for pairs whose carrier table is sparse (expected cell < 5).
    top : int, optional
        Within each population, keep only the N pairs with the largest |r|.

    Returns
    -------
    DataFrame with one row per allele pair per population, sorted by population
    then |r| descending.
    """

    def log(msg: str) -> None:
        if verbose:
            print(msg)

    mat = _as_dosage(dosage, samples_in)
    alleles_all = mat.index.to_numpy()
    samples = mat.columns.to_numpy()
    X_all = mat.to_numpy()
    log(f"[info] {len(alleles_all)} alleles x {len(samples)} samples")

    # ---- loci and per-locus totals ------------------------------------------
    loci_all = _as_loci(locus, alleles_all, locus_regex)
    uniq_loci, codes_all = np.unique(loci_all, return_inverse=True)
    if locus is None and uniq_loci.size == alleles_all.size and alleles_all.size > 1:
        # Nothing matched locus_regex, so each allele became its own locus and
        # include_same_locus has nothing to exclude.
        log(
            f"[warn] no locus parsed from any allele name with {locus_regex!r}; "
            "every allele is its own locus, so within-locus pairs cannot be "
            "dropped. Check the allele names, or pass locus=."
        )
    T_all = _locus_totals(X_all, codes_all, uniq_loci.size)

    # ---- filter alleles once, over all samples ------------------------------
    carriers_all = (X_all > 0).sum(axis=1)
    complete_totals = np.nansum(T_all, axis=1)
    complete_copies = np.where(np.isfinite(T_all), X_all, 0).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        freq_overall = np.where(complete_totals > 0,
                                complete_copies / complete_totals, 0.0)

    keep = (carriers_all >= min_carriers) & (freq_overall >= min_freq)
    if keep.sum() < 2:
        raise ValueError("fewer than two alleles survive filtering")
    log(
        f"[info] {int(keep.sum())} alleles retained "
        f"(min_carriers={min_carriers}, min_freq={min_freq})"
    )

    X = X_all[keep]
    T = T_all[keep]
    alleles = alleles_all[keep]
    loci = loci_all[keep]
    # Carriage, not copy number. Dosage 1 and dosage 2 are both present.
    P = np.where(np.isnan(X), np.nan, (X > 0).astype(np.float64))
    n_alleles = P.shape[0]

    # ---- pair index ----------------------------------------------------------
    ia, ib = np.triu_indices(n_alleles, k=1)
    same_locus = loci[ia] == loci[ib]
    if not include_same_locus:
        sel = ~same_locus
        if not sel.any():
            raise ValueError(
                "no cross-locus pairs; pass include_same_locus=True if intended"
            )
        log(
            f"[info] {int((~sel).sum())} within-locus pairs excluded "
            "(alleles at one locus are not independent in a sample)"
        )
        ia, ib, same_locus = ia[sel], ib[sel], same_locus[sel]

    # ---- groups to compute over ---------------------------------------------
    labels = _as_population(population, samples)
    groups = []
    if labels is None:
        groups.append((ALL_LABEL, np.arange(P.shape[1])))
    else:
        lab = labels.to_numpy()
        for g in pd.unique(lab):
            idx = np.where(lab == g)[0]
            if idx.size < 3:
                log(f"[warn] population {g!r} has {idx.size} samples, skipped")
                continue
            groups.append((str(g), idx))
        if not groups:
            raise ValueError("no population has at least 3 samples")
        log(f"[info] {len(groups)} populations: {', '.join(g for g, _ in groups)}")
        if include_all:
            groups.append((ALL_LABEL, np.arange(P.shape[1])))

    # ---- one independent block per group ------------------------------------
    blocks = []
    for name, idx in groups:
        Ps = P[:, idx]
        n = Ps.shape[1]

        # Pair-specific complete cases: unknown calls never become absences.
        observed = np.isfinite(Ps).astype(np.float64)
        present = np.nan_to_num(Ps)
        absent = observed - present
        pair_n = observed @ observed.T
        n11 = present @ present.T
        n10 = present @ absent.T
        n01 = absent @ present.T
        n00 = absent @ absent.T
        ca, cb = n11 + n10, n11 + n01
        with np.errstate(divide="ignore", invalid="ignore"):
            pa, pb = ca / pair_n, cb / pair_n
            Dcoef = n11 / pair_n - pa * pb
            denominator = np.sqrt(ca * (pair_n - ca) * cb * (pair_n - cb))
            phi = np.where(denominator > 0,
                           (n11 * n00 - n10 * n01) / denominator, np.nan)
            r_pa = np.clip(phi, -1.0, 1.0)
            dmax = np.where(Dcoef >= 0,
                            np.minimum(pa * (1 - pb), pb * (1 - pa)),
                            np.minimum(pa * pb, (1 - pa) * (1 - pb)))
            dp = np.where(dmax > 0, Dcoef / dmax, np.nan)
        p_chi2 = chi2_dist.sf(pair_n * r_pa**2, df=1)

        Xs = X[:, idx]
        totals = T[:, idx]
        # A partially missing locus has an unknown total copy count.
        valid = np.isfinite(Xs) & np.isfinite(totals)
        Ts = np.where(valid, totals, 0).sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            afreq = np.where(Ts > 0, np.where(valid, Xs, 0).sum(axis=1) / Ts, np.nan)

        blk = pd.DataFrame(
            {
                "population": name,
                "n_samples": pair_n[ia, ib].astype(int),
                "n_group_samples": n,
                "allele_a": alleles[ia],
                "allele_b": alleles[ib],
                "locus_a": loci[ia],
                "locus_b": loci[ib],
                "same_locus": same_locus,
                "n_carriers_a": ca[ia, ib].astype(int),
                "n_carriers_b": cb[ia, ib].astype(int),
                "n11": n11[ia, ib].astype(int),
                "n10": n10[ia, ib].astype(int),
                "n01": n01[ia, ib].astype(int),
                "n00": n00[ia, ib].astype(int),
                "carrier_freq_a": pa[ia, ib],
                "carrier_freq_b": pb[ia, ib],
                "allele_freq_a": afreq[ia],
                "allele_freq_b": afreq[ib],
                "r": r_pa[ia, ib],
                "r2": r_pa[ia, ib] ** 2,
                "D": Dcoef[ia, ib],
                "dmax": dmax[ia, ib],
                "Dprime": dp[ia, ib],
                "Wn": phi[ia, ib],
                "p_chi2": p_chi2[ia, ib],
            }
        )

        # Explicit names distinguish carrier association from haplotype LD.
        for legacy in ("r", "r2", "D", "dmax", "Dprime"):
            blk["carrier_" + legacy] = blk[legacy]

        # The carrier-frequency bound contains D. Anything past 1 is rounding.
        blk["dprime_out_of_bounds"] = blk["Dprime"].abs() > (1.0 + 1e-8)
        pair_sizes = blk["n_samples"].replace(0, np.nan)
        blk["min_expected_cell"] = np.minimum.reduce(
            [
                blk["n_carriers_a"] * blk["n_carriers_b"] / pair_sizes,
                blk["n_carriers_a"] * (blk["n_samples"] - blk["n_carriers_b"]) / pair_sizes,
                (blk["n_samples"] - blk["n_carriers_a"]) * blk["n_carriers_b"] / pair_sizes,
                (blk["n_samples"] - blk["n_carriers_a"]) * (blk["n_samples"] - blk["n_carriers_b"]) / pair_sizes,
            ]
        )
        blk["sparse_table"] = blk["min_expected_cell"] < 5

        if exact:
            pe = blk["p_chi2"].to_numpy().copy()
            sparse = np.flatnonzero(blk["sparse_table"].to_numpy())
            for k in sparse:
                tbl = np.array(
                    [
                        [blk.at[k, "n11"], blk.at[k, "n10"]],
                        [blk.at[k, "n01"], blk.at[k, "n00"]],
                    ],
                    dtype=np.int64,
                )
                if tbl.sum(axis=0).min() == 0 or tbl.sum(axis=1).min() == 0:
                    pe[k] = np.nan        # monomorphic carrier status, test undefined
                else:
                    pe[k] = fisher_exact(tbl)[1]
            blk["p_exact"] = pe
            log(f"[info] {name}: Fisher exact for {sparse.size} sparse pairs")

        pcol = "p_exact" if exact else "p_chi2"
        blk["q_value"] = _bh_fdr(blk[pcol].to_numpy())
        blk["p_used_for_q"] = pcol

        blk = blk.reindex(blk["r"].abs().sort_values(ascending=False).index)
        if top:
            blk = blk.head(int(top))
        blocks.append(blk)

        n_oob = int(blk["dprime_out_of_bounds"].sum())
        if n_oob:
            log(
                f"[warn] {name}: |D'| > 1 for {n_oob} pairs after rounding tolerance"
            )
        n_undef = int(blk["Wn"].isna().sum())
        if n_undef:
            log(
                f"[warn] {name}: Wn and r undefined for {n_undef} pairs — carrier "
                "status is monomorphic there"
            )

    out = pd.concat(blocks, ignore_index=True)
    return out


# ----------------------------------------------------------------------- demo


def make_demo_data(n: int = 400, seed: int = 0):
    """Synthetic (dosage, population) with HLA-B/C LD present in EUR only.

    B and C are in strong LD among the EUR samples and independent among the EAS
    samples, so the per-population rows should differ sharply for those pairs while
    the KIR pairs stay near zero everywhere.
    """
    rng = np.random.default_rng(seed)
    half = n // 2
    pops = np.array(["EUR"] * half + ["EAS"] * (n - half))

    pB = [0.40, 0.30, 0.20, 0.10]
    h1 = rng.choice(4, size=n, p=pB)
    h2 = rng.choice(4, size=n, p=pB)

    fidelity = np.where(pops == "EUR", 0.90, 0.0)   # LD in EUR, none in EAS

    def linked(h):
        return np.where(rng.random(n) < fidelity, h % 3, rng.integers(0, 3, n))

    c1, c2 = linked(h1), linked(h2)

    B = np.zeros((4, n)); C = np.zeros((3, n))
    for i in range(n):
        B[h1[i], i] += 1; B[h2[i], i] += 1
        C[c1[i], i] += 1; C[c2[i], i] += 1

    rows = {}
    for i in range(4):
        rows[f"B*{i + 1:02d}:01"] = B[i]
    for i in range(3):
        rows[f"C*{i + 1:02d}:01"] = C[i]
    for j, nm in enumerate(["KIR2DL1", "KIR3DS1", "KIR2DS2"]):
        base = np.where(pops == "EUR", 2, 1) + (1 if j == 1 else 0) * (pops == "EAS")
        rows[f"{nm}*0010101"] = np.clip(base + rng.integers(-1, 2, n), 0, 6).astype(float)
        rows[f"{nm}*0020101"] = np.clip(
            rng.integers(0, 3, n) + (pops == "EAS"), 0, 6
        ).astype(float)

    samples = [f"S{i:03d}" for i in range(n)]
    return pd.DataFrame(rows, index=samples).T, pd.Series(pops, index=samples, name="population")
