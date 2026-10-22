# Design: Publication-ready HLA dosage × ancestry ridgeline plot

Date: 2026-10-22  
Status: approved (Approach 1 — polish-in-place; pending user review of this file)

## Goal

Refine the existing seaborn `FacetGrid` ridgeline (joyplot) cell so ancestry labels are fully readable, the dosage axis focuses on the data range, and the figure is journal-ready (colorblind-safe palette, clean theme, PDF+PNG export).

## Non-goals

- Changing how `hla_dosage` or `ancestry_pred` are computed
- Switching to joypy / pure matplotlib / non-overlapping facets
- Building a reusable Python package/theme
- Wiring into WDL (optional later)

## User choices (locked)

| Decision | Choice |
|----------|--------|
| Approach | Polish-in-place (keep overlapping FacetGrid ridgeline) |
| Label bug | Fix clipping so full `ancestry_pred` names are visible |
| X-axis | Fixed `xlim=(60, 80)` |
| Palette | Colorblind-safe discrete (Okabe–Ito), keyed by ancestry |

## Approach

**Polish-in-place:** keep `sns.FacetGrid` + dual `kdeplot` + `refline` + overlap; fix label placement/margins, xlim, palette, labels/theme, and export.

## Behavior

### Inputs

- `df_meta` with at least `ancestry_pred` (categorical/string) and `hla_dosage` (numeric).
- Optional: `row_order` / `hue_order` if a fixed ancestry order is desired; otherwise seaborn categorical order.

### Label fix (root cause)

Current `ax.text(0, .2, label, …, transform=ax.transAxes)` places text at the left edge of each axes; with overlap and default margins, names are clipped (only fragments / `"oth"` visible).

1. Draw labels with `clip_on=False`.
2. Place text slightly inside or just left of the axes in axes coords (e.g. `x=0.01` or `x=-0.02`) with `ha="left"`, `va="center"`, `fontweight="bold"`, color from the hue mapping.
3. Reserve left figure margin so labels are not cut by the figure boundary: e.g. `g.figure.subplots_adjust(left=0.12, hspace=-0.25)` (tune `left` if longest ancestry name needs more room). Prefer adjusting `left` over relying on `tight_layout`, which fights negative `hspace`.

### X-axis

1. After mapping KDEs: `g.set(xlim=(60, 80))`.
2. X label: `"HLA dosage"` (human-readable; not raw `hla_dosage`).
3. Keep bottom spine for the shared x-axis; no y ticks/labels.

### Densities

1. Filled KDE: `bw_adjust=0.5`, `clip_on=False`, `fill=True`, `alpha=1`, `linewidth=1.5`.
2. White outline KDE on top (same `bw_adjust`).
3. `g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)`.
4. Overlap: `hspace=-0.25` (unchanged unless labels collide; then ease toward `-0.15`).

### Palette & order

1. Okabe–Ito discrete colors, recycled if more ancestries than colors:
   `#E69F00`, `#56B4E9`, `#009E73`, `#F0E442`, `#0072B2`, `#D55E00`, `#CC79A7`, `#999999`.
2. Pass as `palette` to `FacetGrid` with `hue="ancestry_pred"` (same as row).
3. Optional canonical order when present (e.g. AFR, AMR, EAS, EUR, MID, SAS, other) — append any unseen labels alphabetically. If order is unknown in the notebook, omit and use data order.

### Theme

1. `g.set_titles("")`.
2. `g.set(yticks=[], ylabel="")`.
3. `g.despine(bottom=True, left=True)`.
4. Optional title: `"HLA dosage by ancestry"` (omit if empty).
5. No color legend (row labels carry the identity).

### Export

1. Basename: `"HLA_dosage_ridgeline"` (optionally prefix with allele/locus string if the notebook defines one).
2. Write **PDF** (vector) and **PNG** at 300 dpi.
3. Size: width ~10 in; height ≈ `0.5 * n_ancestries + 1.5` in (or FacetGrid `height=0.5`, `aspect=15` as now, then `savefig` from `g.figure`).

## Deliverable

Drop-in replacement for the user’s notebook cell (same inputs: `df_meta` with `ancestry_pred`, `hla_dosage`). Optionally also save as `haplograph/scripts/plot_hla_dosage_ridgeline.py` for reuse — only if requested; default is notebook cell only.

## Out of scope / later

- Sample-size annotations on each ridge
- Per-allele faceting or multi-panel dosage plots
- Port into WDL figure tasks

## Spec self-review

- No TBD placeholders.
- Label fix (clip + margin), xlim `(60, 80)`, Okabe–Ito, and export are explicit.
- Scope is one notebook cell; matches Approach A polish-in-place.
- Ambiguity resolved: human x-label `"HLA dosage"`; optional ancestry order is best-effort, not required for the first cell.
