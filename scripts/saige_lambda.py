#!/usr/bin/env python3
"""
Genomic inflation per phecode from SAIGE step 2 single-variant output.

lambda      = median(chi2_1df) / 0.4549, from SAIGE's p.value column
lambda1000  = 1 + (lambda - 1) * 500 * (1/n_cases + 1/n_controls)
              i.e. lambda rescaled to a study of 1,000 cases and 1,000 controls,
              so phecodes with very different case counts can be compared
              (Freedman et al. 2004; de Bakker et al. 2008).

Also reports lambda within MAC bins (where variance-ratio or SPA problems show
up) and, if given, the step 1 variance ratios per MAC category.

Notebook use:
    from saige_lambda import read_saige_singlevar, inflation_table
    df = read_saige_singlevar("results/*_singlevar")
    tab = inflation_table(df)

Command line:
    python saige_lambda.py --step2 'results/*_singlevar' \
        [--vr 'results/*.varianceRatio.txt'] [--phecode-info phecodeX_info.csv] \
        --out inflation.tsv

Caveat: with only a few hundred HLA/KIR markers, which are strongly correlated
(LD, and alleles of one gene are mutually exclusive), lambda is noisy and is
pushed up by real MHC signal. Check n_markers, and prefer a genome-wide null
marker set run through the same model when you need a calibration number.
"""

import argparse
import glob
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

CHI2_MEDIAN_1DF = stats.chi2.ppf(0.5, 1)  # 0.454936...
STEP2_RE = re.compile(r"_step2Out_(?P<phecode>.+)_(?P<chrom>(?:chr)?[0-9XYMT]+)_singlevar$")
STEP1_VR_RE = re.compile(r"_step1Out_(?P<phecode>.+)\.varianceRatio\.txt$")


def lambda_gc(pvalues):
    """Median-based genomic inflation factor from p-values (1 df)."""
    p = np.asarray(pvalues, dtype=float)
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]  # SAIGE writes 0 only on underflow
    if p.size == 0:
        return np.nan
    chi2 = stats.chi2.isf(p, 1)
    return float(np.median(chi2) / CHI2_MEDIAN_1DF)


def lambda_1000(lam, n_cases, n_controls):
    """Rescale lambda to 1,000 cases and 1,000 controls."""
    if not (np.isfinite(lam) and n_cases > 0 and n_controls > 0):
        return np.nan
    return float(1 + (lam - 1) * 500 * (1.0 / n_cases + 1.0 / n_controls))


def read_saige_singlevar(paths):
    """Read SAIGE step 2 single-variant files; add phecode and MAC columns.

    `paths` is a glob string or a list of paths. The phecode is parsed from the
    SAIGE_phewas.wdl naming <prefix>_step2Out_<phecode>_<chrom>_singlevar.
    """
    files = sorted(glob.glob(paths)) if isinstance(paths, str) else list(paths)
    files = [f for f in files if not f.endswith(".index")]
    frames = []
    for f in files:
        m = STEP2_RE.search(os.path.basename(f))
        if not m or os.path.getsize(f) == 0:
            continue
        d = pd.read_csv(f, sep=r"\s+")
        if d.empty:
            continue
        d["phecode"] = m.group("phecode")
        frames.append(d)
    if not frames:
        raise FileNotFoundError(f"no SAIGE single-variant results matched {paths!r}")
    df = pd.concat(frames, ignore_index=True)

    if {"N_case", "N_ctrl"} <= set(df.columns):
        n = df["N_case"] + df["N_ctrl"]
    elif "N" in df.columns:
        n = df["N"]
    else:
        n = None
    if n is not None and "AC_Allele2" in df.columns:
        df["MAC"] = np.minimum(df["AC_Allele2"], 2 * n - df["AC_Allele2"])
    return df


def inflation_table(df, mac_bins=(20, 50, 100, 500, np.inf), pcol="p.value"):
    """One row per phecode: n, lambda, lambda1000 and lambda per MAC bin."""
    rows = []
    for phecode, g in df.groupby("phecode", sort=False):
        lam = lambda_gc(g[pcol])
        n_case = g["N_case"].median() if "N_case" in g else np.nan
        n_ctrl = g["N_ctrl"].median() if "N_ctrl" in g else np.nan
        row = {
            "phecode": phecode,
            "n_markers": int(g[pcol].notna().sum()),
            "n_cases": n_case,
            "n_controls": n_ctrl,
            "lambda": lam,
            "lambda1000": lambda_1000(lam, n_case, n_ctrl),
            "min_p": g[pcol].min(),
        }
        if "MAC" in g:
            edges = list(mac_bins)
            for lo, hi in zip(edges[:-1], edges[1:]):
                sub = g[(g["MAC"] >= lo) & (g["MAC"] < hi)]
                label = f"MAC{lo:g}-{hi:g}" if np.isfinite(hi) else f"MAC>={lo:g}"
                row[f"lambda_{label}"] = lambda_gc(sub[pcol])
                row[f"n_{label}"] = len(sub)
        rows.append(row)
    out = pd.DataFrame(rows)
    return out.sort_values("lambda1000", ascending=False, na_position="last")


def read_variance_ratios(paths):
    """Step 1 variance ratios, wide: one row per phecode, columns <type>_cat<k>.

    The file has rows `ratio  type  k` with type null / null_noXadj / sparse and
    k = MAC category (1 = lowest). SAIGE writes exactly 1 for a category it could
    not estimate.
    """
    files = sorted(glob.glob(paths)) if isinstance(paths, str) else list(paths)
    rows = []
    for f in files:
        m = STEP1_VR_RE.search(os.path.basename(f))
        if not m or os.path.getsize(f) == 0:
            continue
        # keep_default_na=False: pandas would otherwise read the type "null" as NaN
        vr = pd.read_csv(f, sep=r"\s+", header=None, names=["ratio", "type", "cat"],
                         keep_default_na=False, dtype={"type": str})
        rec = {"phecode": m.group("phecode")}
        for _, r in vr.iterrows():
            rec[f"{r['type']}_cat{int(r['cat'])}"] = float(r["ratio"])
        rows.append(rec)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--step2", required=True, help="glob for *_singlevar files (quote it)")
    ap.add_argument("--vr", help="glob for *.varianceRatio.txt files (quote it)")
    ap.add_argument("--phecode-info",
                    help="CSV/TSV with a 'phecode' column (e.g. phecodeX info with "
                         "'sex', 'category') to merge into the table")
    ap.add_argument("--out", required=True, help="output TSV")
    args = ap.parse_args()

    tab = inflation_table(read_saige_singlevar(args.step2))
    if args.vr:
        vr = read_variance_ratios(args.vr)
        if not vr.empty:
            tab = tab.merge(vr, on="phecode", how="left")
    if args.phecode_info:
        sep = "\t" if args.phecode_info.endswith((".tsv", ".txt")) else ","
        info = pd.read_csv(args.phecode_info, sep=sep, dtype={"phecode": str})
        tab = tab.merge(info, on="phecode", how="left")
    tab.to_csv(args.out, sep="\t", index=False, float_format="%.4g")
    print(tab.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
