#!/usr/bin/env python3
"""Bonferroni threshold from a PheWAS results table.

A PheWAS table has one row per phenotype tested (PheWAS R package `p`,
PheTK `p_value`, or SAIGE `p.value`). Phenotypes that were not tested have
a missing p-value and are not counted. The family-wise threshold is

    p < alpha / n_tests

which is the same line PheTK draws at -log10(alpha / n_tests).

Notebook use:
    from phewas_threshold import bonferroni_threshold
    thr = bonferroni_threshold(results)          # results: one row per phecode
    # thr["p_threshold"], thr["neglog10_threshold"]

If the table stacks several markers (one row per marker x phecode), pass
`group` so each marker is its own family of phenotype tests:

    bonferroni_threshold(saige, group="MarkerID")
"""

import numpy as np
import pandas as pd

# Column names written by the PheWAS R package, PheTK, and SAIGE step 2.
_P_COLUMNS = ("p", "p.value", "p_value", "pval", "pvalue", "P")


def _p_column(df, pcol):
    if pcol is not None:
        if pcol not in df.columns:
            raise KeyError(f"p-value column {pcol!r} not in {list(df.columns)}")
        return pcol
    found = [c for c in _P_COLUMNS if c in df.columns]
    if len(found) == 1:
        return found[0]
    if len(found) > 1:
        raise ValueError(
            f"several p-value columns {found}; pass pcol= to choose one"
        )
    raise KeyError(
        "no p-value column found; expected one of "
        + ", ".join(_P_COLUMNS)
        + f", got {list(df.columns)}"
    )


def _one_threshold(p, alpha):
    """Threshold from a 1-d array of p-values. Missing tests are dropped."""
    p = np.asarray(p, dtype=float)
    n = int(np.isfinite(p).sum())
    if n == 0:
        return {
            "n_tests": 0,
            "alpha": float(alpha),
            "p_threshold": np.nan,
            "neglog10_threshold": np.nan,
        }
    p_threshold = float(alpha) / n
    return {
        "n_tests": n,
        "alpha": float(alpha),
        "p_threshold": p_threshold,
        "neglog10_threshold": float(-np.log10(p_threshold)),
    }


def bonferroni_threshold(phewas, pcol=None, alpha=0.05, group=None):
    """Bonferroni threshold for the tests in a PheWAS table.

    Parameters
    ----------
    phewas : pandas.DataFrame
        One row per test. Rows with a missing p-value were not tested and
        are excluded from `n_tests`.
    pcol : str, optional
        P-value column. Inferred when the table has exactly one of `p`,
        `p.value`, `p_value`, `pval`, `pvalue`, `P`.
    alpha : float
        Family-wise error rate. Default 0.05, so the threshold is 0.05 / n.
    group : str, optional
        Column that splits the table into separate families (for example
        `MarkerID` when each allele was tested against the same phecodes).
        One threshold is computed per group. Omit it to treat every finite
        p-value in the table as one family.

    Returns
    -------
    dict
        ``n_tests``, ``alpha``, ``p_threshold``, ``neglog10_threshold``.
        A result is significant when its p-value is below ``p_threshold``,
        or when -log10(p) exceeds ``neglog10_threshold``.
    pandas.DataFrame
        Same fields, one row per group, when `group` is given.
    """
    if not 0 < alpha <= 1:
        raise ValueError(f"alpha must be in (0, 1]; got {alpha}")
    if not isinstance(phewas, pd.DataFrame):
        raise TypeError(f"phewas must be a DataFrame; got {type(phewas).__name__}")

    pcol = _p_column(phewas, pcol)
    if group is None:
        return _one_threshold(phewas[pcol], alpha)

    if group not in phewas.columns:
        raise KeyError(f"group column {group!r} not in {list(phewas.columns)}")
    rows = []
    for key, part in phewas.groupby(group, sort=False, dropna=False):
        rec = {group: key}
        rec.update(_one_threshold(part[pcol], alpha))
        rows.append(rec)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    demo = pd.DataFrame({"phecode": ["A", "B", "C", "D"], "p": [0.01, 0.2, np.nan, 1e-6]})
    thr = bonferroni_threshold(demo)
    # D was not tested, so n = 3 and the line is -log10(0.05 / 3).
    assert thr["n_tests"] == 3
    assert np.isclose(thr["p_threshold"], 0.05 / 3)
    assert np.isclose(thr["neglog10_threshold"], -np.log10(0.05 / 3))

    stacked = pd.DataFrame({
        "MarkerID": ["m1", "m1", "m2", "m2"],
        "p.value": [0.01, 0.2, np.nan, 0.4],
    })
    per_marker = bonferroni_threshold(stacked, group="MarkerID")
    assert list(per_marker["n_tests"]) == [2, 1]
    assert np.isclose(per_marker.loc[0, "p_threshold"], 0.05 / 2)
    print(thr)
