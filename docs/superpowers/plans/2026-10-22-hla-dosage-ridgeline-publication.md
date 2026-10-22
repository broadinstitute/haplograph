# HLA Dosage Ridgeline Publication Plot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Deliver a publication-ready drop-in seaborn ridgeline cell matching the approved spec (label fix, xlim 60–80, Okabe–Ito, PDF+PNG).

**Architecture:** Single Python script / notebook cell: FacetGrid ridgeline by `ancestry_pred`, dual kdeplot, fixed left margin for labels, `xlim=(60, 80)`, export from `g.figure`.

**Tech Stack:** Python, seaborn, matplotlib, pandas

## Global Constraints

- Labels: `clip_on=False`; left margin reserved; full `ancestry_pred` visible
- X-axis: `xlim=(60, 80)`; label `"HLA dosage"`
- Colors: Okabe–Ito by ancestry; no legend
- Overlap: `hspace=-0.25`
- Export: PDF + PNG @ 300 dpi; basename `HLA_dosage_ridgeline`

---

### Task 1: Write drop-in Python cell

**Files:**
- Create: `haplograph/scripts/plot_hla_dosage_ridgeline.py` (canonical cell body; paste into notebook as needed)

**Interfaces:**
- Consumes: `df_meta` with columns `ancestry_pred`, `hla_dosage`
- Produces: `g` (FacetGrid); files `HLA_dosage_ridgeline.pdf` / `.png`

- [x] **Step 1: Implement ridgeline plot per spec** → `scripts/plot_hla_dosage_ridgeline.py`
- [x] **Step 2: Smoke-check with synthetic data** (`MPLBACKEND=Agg python scripts/plot_hla_dosage_ridgeline.py` → OK, xlim 60–80)
- [x] **Step 3: Commit** script + plan
## Spec coverage self-review

- Label fix, xlim, Okabe–Ito, theme, export → Task 1
- No placeholders
