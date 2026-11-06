#!/usr/bin/env python3
"""
Extract complete per-haplotype gene structures, block by block, from one
Immuannot annotation table (e.g. hprc.immuannot.csv).

This is `collect_haplotype_features` from hap_structure_plot.py repurposed: instead
of grouping haplotypes into ranked structures for plotting, it keeps every
(sample, haplotype) and reports, for each of eight gene blocks, the gene order of
that block on that haplotype -- or 0 when the block is not completely assembled.

Definition of complete
----------------------
A block is complete on a haplotype when ONE contig carries both of the block's
anchor genes. The anchors are always-present (copy-invariant) genes at the two
ends of the block, so "both anchors on one contig" means the whole block was
traversed by a single contig. Only genes lying between the anchors on that contig
are reported, and the structure is oriented so the first anchor reads on the left.

Anchors are not always the first/last gene of the block's gene list: the first
gene of the DR block (HLA-DRB5) varies in copy number and so cannot witness
completeness; DR is anchored on DRB1..DRA instead.

Haplotypes with no such contig get structure "0" and a reason
(anchor_missing:<gene> / anchors_on_different_contigs).

RCCX is the exception: it has "anchors": None. The block is nothing but
copy-variable C4 genes, so no gene's presence can witness a traversal. Every C4
copy on the haplotype is reported for inspection, across contigs if need be.
These observations remain incomplete unless --assume-complete is supplied.
Multiple contigs cannot establish gene order and remain incomplete regardless.
Such a block is listed by contig then position, flipped when most copies read
on the minus strand.

Resolution
----------
Structures are reported twice per haplotype:
  structure         gene names only, comma-joined left to right
                    e.g. "HLA-F,HLA-V,HLA-P,HLA-G,HLA-A,HLA-J"
  structure_allele  the same order at the requested allele resolution
                    e.g. "HLA-F*01,HLA-V*01,...", KIR as "KIR2DL1*001"
HLA resolution is `--hla-digits` (0, 2, 4, 6, or 8; default 2), two digits per
colon field: 0 -> HLA-A, 2 -> HLA-A*01, 4 -> HLA-A*01:01. KIR resolution is
`--kir-digits` (0, 3, 5, or 7; default 3): 0 -> KIR2DL1, 3 -> KIR2DL1*001,
5 -> KIR2DL1*00101. Digit 0 is the gene-level haplotype: the allele is dropped
and `structure_allele` matches `structure` for that gene family. A call shorter
than a positive resolution keeps the bare gene, so a 2-digit HLA name does not
appear inside a 4-digit structure. IMGT expression suffixes such as N are kept;
Immuannot's "new" novel tag is dropped.

Input
-----
The same table hap_structure_plot.py takes, with columns
    Sample  Hap  Gene  Strand  Consensus_Allele  Contig  Spos  Epos  [Tool]
`Hap` identifies the haplotype within a sample; pass hap_col=None to treat each
contig as its own unit instead. `Gene` may be a bare gene name or an allele
("HLA-A*01:01"); anything after "*" is stripped.

Library use
-----------
    import pandas as pd
    from hap_block_haplotypes import collect_block_haplotypes, wide_structures

    annot = pd.read_csv("hprc.immuannot.csv", index_col=0)
    long = collect_block_haplotypes(annot)                 # one row per hap x block
    wide = wide_structures(long)                           # one row per hap, 8 columns
    alpha = long[(long["block"] == "alpha") & long["complete"].eq(1)]

CLI
---
    python hap_block_haplotypes.py --annot hprc.immuannot.csv --out hprc_blocks
writes <out>.block_haplotypes.tsv (long), <out>.block_structures_wide.tsv
(gene order), <out>.block_haplotypes_wide.tsv (the resolution above),
<out>.block_structure_counts.tsv, and one dosage VCF of that resolved
haplotype. The VCF matches allele_dosage_vcf.py: DS is a float dosage and
the haplotype is the HI field. The ID is the block name plus its index
(DR_1, KIR_2, ...).
Additive dosage (0, 1, 2) is the default; --dominant writes carrier status.
Unresolved dosages are '.'. RCCX observed structures remain partial unless
--assume-complete certifies the annotations; multiple-contig RCCX order remains
unresolved even with that assumption.
"""
from __future__ import annotations

import argparse
import sys
from datetime import date

import numpy as np
import pandas as pd

from allele_dosage_vcf import candidate_alleles, collapse_copies

# --------------------------------------------------------------------------- #
# Blocks                                                                       #
# --------------------------------------------------------------------------- #
# genes   : the genes reported in the structure (as requested, gene level)
# anchors : (left, right) copy-invariant genes that must share one contig for the
#           block to count as complete. The structure is oriented so `left` comes
#           first, so the pair also fixes the reading direction: each one runs in
#           the direction of the gene list it belongs to (alpha HLA-F -> HLA-J,
#           beta HLA-C/B -> MICB, DR HLA-DRB1 -> HLA-DRA, ...).
BLOCKS = {
    "alpha": {
        "genes": ["HLA-F", "HLA-V", "HLA-P", "HLA-G", "HLA-H", "HLA-T", "HLA-K",
                  "HLA-U", "HLA-A", "HLA-Y", "HLA-W", "HLA-J"],
        "anchors": ("HLA-F", "HLA-J"),
    },
    "beta": {
        "genes": ["HLA-B", "HLA-C", "HLA-S", "MICA", "HLA-X", "MICB"],
        "anchors": ("HLA-C", "MICB"),
    },
    "RCCX": {
        # No anchors: the block is nothing but copy-variable C4 genes, so there is
        # no gene whose presence could witness a complete traversal. Every C4 copy
        # on the haplotype counts; a haplotype with no C4 at all is 0.
        "genes": ["C4AL", "C4BL", "C4AS", "C4BS"],
        "anchors": None,
    },
    "DR": {
        "genes": ["HLA-DRB5", "HLA-DRB3", "HLA-DRB4", "HLA-DRB1", "HLA-DRA"],
        "anchors": ("HLA-DRB1", "HLA-DRA"),
    },
    "delta": {
        "genes": ["HLA-DOB", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2"],
        "anchors": ("HLA-DQA1", "HLA-DOB"),
    },
    "epsilon": {
        "genes": ["HLA-DMB", "HLA-DPB2", "HLA-DPA2", "HLA-DPB1", "HLA-DOA",
                  "HLA-DPA1", "HLA-DMA"],
        "anchors": ("HLA-DMB", "HLA-DPB2"),
    },
    "kappa": {
        "genes": ["HLA-L", "HLA-N", "HLA-E"],
        "anchors": ("HLA-L", "HLA-E"),
    },
    "KIR": {
        "genes": ["KIR3DP1", "KIR2DL4", "KIR3DS1", "KIR2DS1", "KIR3DL2",
                  "KIR2DL5B", "KIR2DL3", "KIR2DL5A", "KIR2DS2", "KIR2DS5",
                  "KIR2DS4", "KIR2DL2", "KIR2DS3", "KIR2DP1", "KIR3DL1",
                  "KIR3DL3", "KIR2DL1"],
        "anchors": ("KIR3DL3", "KIR3DL2"),
    },
}

BLOCK_ORDER = list(BLOCKS)

# Single-letter IMGT expression suffixes (not Immuannot's "new" novel tag).
IMGT_EXPR_SUFFIX = frozenset("NLSQCA")

# Same resolution knobs as allele_dosage_vcf.py. A haplotype is one string, so
# KIR keeps a single length here rather than the nested 3/5/7 markers that
# script emits for one gene copy.
HLA_DIGITS = 2
KIR_DIGITS = 3
# 0 is gene level: the allele is dropped and the haplotype is gene names only.
HLA_DIGIT_CHOICES = (0, 2, 4, 6, 8)
KIR_DIGIT_CHOICES = (0, 3, 5, 7)

# Immuannot calls the column Consensus_Allele. GTF exports of the same
# annotation call it Alleles or template_allele.
ALLELE_COL_CANDIDATES = (
    "Consensus_Allele", "Allele", "Alleles", "allele",
    "template_allele", "Template_Allele", "consensus",
)

MISSING_VALUE = "0"  # structure of a haplotype that does not span the block

# Left-anchor coordinate used only as the VCF position label. Later structures
# of the same block take the next free base, then a per-chromosome pass keeps
# every POS unique.
BLOCK_ANCHOR = {
    "alpha": ("chr6", 29723434),      # HLA-F
    "kappa": ("chr6", 30259548),      # HLA-L
    "beta": ("chr6", 31268749),       # HLA-C
    "RCCX": ("chr6", 31982057),       # C4A
    "DR": ("chr6", 32439872),         # HLA-DRA
    "delta": ("chr6", 32637406),      # HLA-DQA1
    "epsilon": ("chr6", 32934629),    # HLA-DMB
    "KIR": ("chr19", 54724442),       # KIR3DL3
}
CONTIG_LENGTH = {
    "chr6": 170805979,
    "chr11": 135086622,
    "chr19": 58617616,
}


# --------------------------------------------------------------------------- #
# Allele names                                                                 #
# --------------------------------------------------------------------------- #
def trim_allele(gene, allele, hla_digits=HLA_DIGITS, kir_digits=KIR_DIGITS):
    """`gene`*<allele at the requested resolution>, or just `gene`.

    hla_digits : 0, 2, 4, 6, or 8. Two digits are one colon field
                 (0 -> HLA-A, 2 -> HLA-A*01, 4 -> HLA-A*01:01,
                 8 -> HLA-A*01:01:01:01).
    kir_digits : 0, 3, 5, or 7 leading digits of a KIR name
                 (0 -> KIR2DL1, 3 -> KIR2DL1*001, 7 -> KIR2DL1*0010101).
    KIR genes take `kir_digits`; every other gene takes `hla_digits`. Digit 0
    drops the allele. A call with fewer fields or digits than a positive
    resolution contributes no allele. IMGT expression suffixes are kept and
    Immuannot's "new" tag is dropped.
    """
    if hla_digits not in HLA_DIGIT_CHOICES:
        raise ValueError(
            f"hla_digits must be one of {HLA_DIGIT_CHOICES}, got {hla_digits}")
    if kir_digits not in KIR_DIGIT_CHOICES:
        raise ValueError(
            f"kir_digits must be one of {KIR_DIGIT_CHOICES}, got {kir_digits}")
    if str(gene).startswith("C4"):
        return gene
    if str(gene).startswith("KIR"):
        if kir_digits == 0:
            return gene
    elif hla_digits == 0:
        return gene
    if allele is None or (isinstance(allele, float) and np.isnan(allele)):
        return gene
    name = str(allele).strip()
    if not name or name.lower() in ("na", "nan", "none", "."):
        return gene
    name = name.split("*", 1)[1] if "*" in name else name
    if name.endswith(":new"):
        name = name[:-4]
    elif name.endswith("new"):
        name = name[:-3]
    if not name:
        return gene
    suffix = name[-1] if name[-1] in IMGT_EXPR_SUFFIX else ""
    core = name[:-1] if suffix else name
    if not core:
        return gene
    if str(gene).startswith("KIR"):
        if len(core) < kir_digits:
            return gene
        trimmed = core[:kir_digits]
    else:
        # Filename-safe exports write 01:02:01:01 as 01_02_01_01. At 4 digits
        # a single underscore token has one field and used to fall back to
        # the bare gene, so the whole haplotype stayed at gene level.
        if ":" not in core and "_" in core:
            core = core.replace("_", ":")
        n_fields = hla_digits // 2
        fields = core.split(":")
        if len(fields) < n_fields or not all(fields[:n_fields]):
            return gene
        trimmed = ":".join(fields[:n_fields])
    return f"{gene}*{trimmed}{suffix}"


# --------------------------------------------------------------------------- #
# Extraction                                                                   #
# --------------------------------------------------------------------------- #
ANNOT_COLS = ("Sample", "Gene", "Contig", "Spos", "Epos")


def _allele_series(annot, allele_col, gene_col):
    """The allele column, or None when the table has no allele calls.

    A missing `Consensus_Allele` is looked up under the GTF names (`Alleles`,
    `template_allele`). Alleles written into the gene column (`HLA-A*01:01`)
    are used when no allele column exists. An absent call is reported, because
    every haplotype is then a gene name no matter what `--hla-digits` is.
    """
    by_lower = {str(c).lower(): c for c in annot.columns}
    chosen = None
    if allele_col and allele_col in annot.columns:
        chosen = allele_col
    elif allele_col and str(allele_col).lower() in by_lower:
        chosen = by_lower[str(allele_col).lower()]
    else:
        for cand in ALLELE_COL_CANDIDATES:
            if cand.lower() in by_lower:
                chosen = by_lower[cand.lower()]
                print(f"allele column {allele_col!r} is not in the table; "
                      f"using {chosen!r}", file=sys.stderr)
                break
    if chosen is not None:
        series = annot[chosen]
        called = series.notna() & ~series.astype(str).str.strip().str.lower().isin(
            ["", "na", "nan", "none", "."])
        if called.any():
            return series
        print(f"allele column {chosen!r} has no allele calls", file=sys.stderr)
    if annot[gene_col].astype(str).str.contains(r"\*").any():
        print(f"reading alleles from {gene_col!r}", file=sys.stderr)
        return annot[gene_col]
    print("no allele calls found; every haplotype is a gene name",
          file=sys.stderr)
    return None


def _normalise(annot, sample_col, hap_col, gene_col, allele_col,
               contig_col, start_col, end_col, strand_col):
    """Rename the annotation table to the internal column names and type it."""
    need = [sample_col, gene_col, contig_col, start_col, end_col]
    if hap_col is not None:
        need.append(hap_col)
    missing = [c for c in need if c not in annot.columns]
    if missing:
        raise ValueError(f"annotation table is missing columns: {missing}")

    d = pd.DataFrame({
        "sample": annot[sample_col].astype(str).str.strip(),
        "contig": annot[contig_col].astype(str).str.strip(),
        "spos": annot[start_col].astype(np.int64),
        "epos": annot[end_col].astype(np.int64),
    }, index=annot.index)
    d["hap"] = (annot[hap_col].astype(str).str.strip() if hap_col is not None
                else d["contig"])
    # "HLA-A*01:01" -> "HLA-A"; a bare gene name is left alone
    d["gene"] = annot[gene_col].astype(str).str.strip().str.split("*").str[0]
    d["allele"] = _allele_series(annot, allele_col, gene_col)
    if "_candidate_alleles" in annot:
        d["_candidate_alleles"] = annot["_candidate_alleles"]
    d["strand"] = (annot[strand_col].astype(str).str.strip()
                   if strand_col in annot.columns else ".")
    d["hap_id"] = d["sample"] + "\t" + d["hap"]
    return d


def _block_rows_anchorless(d, genes, hla_digits, kir_digits, orient_strand="+",
                           assume_complete=False):
    """One record per haplotype for a block with no anchor genes (RCCX).

    Nothing can witness a complete traversal here, so every copy found on the
    haplotype is reported, across contigs if need be. Completeness must be
    explicitly certified and order requires one contig. Order is by contig then
    position, flipped when most copies read on the opposite strand (C4 is + in GRCh38).
    """
    sub = d[d["gene"].isin(set(genes))].copy()
    records = []
    if len(sub):
        f = (sub.groupby("hap_id")["strand"]
             .apply(lambda s: (s == orient_strand).sum() < (s != orient_strand).sum()))
        sub["flipped"] = sub["hap_id"].map(f)
        sub["key"] = np.where(sub["flipped"].to_numpy(bool), -sub["epos"], sub["spos"])
        sub["label_allele"] = sub["resolved_allele"]
        sub = sub.sort_values(["hap_id", "contig", "key", "epos"])
        agg = sub.groupby("hap_id").agg(
            structure=("gene", ",".join),
            structure_allele=("label_allele", ",".join),
            n_genes=("gene", "size"),
            allele_ok=("allele_ok", "min"),
            contigs=("contig", lambda c: ",".join(sorted(set(c)))),
            flipped=("flipped", "first"),
            lo=("spos", "min"), hi=("epos", "max"))
    else:
        agg = pd.DataFrame()

    for hid in sorted(d["hap_id"].unique()):
        sample, hap = hid.split("\t", 1)
        if hid in agg.index:
            a = agg.loc[hid]
            one_contig = "," not in a["contigs"]
            records.append(dict(
                sample=sample, hap=hap, complete=int(assume_complete and one_contig),
                allele_complete=int(a["allele_ok"] and assume_complete and one_contig),
                contig=a["contigs"],
                n_genes=int(a["n_genes"]), structure=a["structure"],
                structure_allele=a["structure_allele"],
                span_bp=int(a["hi"] - a["lo"]) if one_contig else 0,
                flipped=bool(a["flipped"]),
                reason=("assumed_complete" if assume_complete else "completeness_unknown")
                       if one_contig else "order_unknown:multiple_contigs"))
        else:
            records.append(dict(
                sample=sample, hap=hap, complete=int(assume_complete),
                allele_complete=int(assume_complete), contig="", n_genes=0,
                structure="" if assume_complete else MISSING_VALUE,
                structure_allele="" if assume_complete else MISSING_VALUE,
                span_bp=0, flipped=False, reason="no_gene_of_block_found"))
    return records


def _block_rows(d, genes, anchors, hla_digits, kir_digits):
    """One record per haplotype for a single block (see collect_block_haplotypes)."""
    a_left, a_right = anchors
    gene_set = set(genes)
    wanted = gene_set | {a_left, a_right}
    sub = d[d["gene"].isin(wanted)]

    # anchor spans per (hap_id, contig)
    def span(gene):
        g = sub[sub["gene"] == gene].groupby(["hap_id", "contig"])
        return pd.concat({f"{gene}_s": g["spos"].min(), f"{gene}_e": g["epos"].max()},
                         axis=1)

    left, right = span(a_left), span(a_right)
    both = left.join(right, how="inner").reset_index()

    # how many of the block's own genes sit on each candidate contig
    own = (sub[sub["gene"].isin(gene_set)]
           .groupby(["hap_id", "contig"])["gene"].nunique().rename("n_block_genes"))
    both = both.merge(own, on=["hap_id", "contig"], how="left")
    both["n_block_genes"] = both["n_block_genes"].fillna(0).astype(int)
    # one contig per haplotype: the one carrying most of the block
    both = (both.sort_values(["hap_id", "n_block_genes", "contig"])
            .drop_duplicates("hap_id", keep="last").set_index("hap_id"))

    ls, le = both[f"{a_left}_s"], both[f"{a_left}_e"]
    rs, re_ = both[f"{a_right}_s"], both[f"{a_right}_e"]
    fwd = ls <= rs                     # left anchor already on the left
    both["lo"] = np.where(fwd, ls, rs)
    both["hi"] = np.where(fwd, re_, le)
    both["flipped"] = ~fwd

    # block genes inside the anchor span on the chosen contig
    g = sub[sub["gene"].isin(gene_set)].merge(
        both[["contig", "lo", "hi", "flipped"]].reset_index(),
        on=["hap_id", "contig"], how="inner")
    g = g[(g["spos"] >= g["lo"]) & (g["epos"] <= g["hi"])].copy()
    f = g["flipped"].to_numpy(bool)
    g["key"] = np.where(f, -g["epos"], g["spos"])
    g["rel_start"] = np.where(f, g["hi"] - g["epos"], g["spos"] - g["lo"])
    g["label_gene"] = g["gene"]
    g["label_allele"] = g["resolved_allele"]
    g = g.sort_values(["hap_id", "key", "epos"])

    agg = g.groupby("hap_id").agg(
        structure=("label_gene", ",".join),
        structure_allele=("label_allele", ",".join),
        n_genes=("label_gene", "size"),
        allele_ok=("allele_ok", "min"))

    # reasons for the haplotypes that are not complete
    has = {a: set(sub.loc[sub["gene"] == a, "hap_id"]) for a in (a_left, a_right)}
    contigs_of = (sub[sub["gene"].isin((a_left, a_right))]
                  .groupby("hap_id")["gene"].apply(set))

    def reason(hid):
        absent = [a for a in (a_left, a_right) if hid not in has[a]]
        if absent:
            return "anchor_missing:" + "+".join(absent)
        if len(contigs_of.get(hid, set())) == 2:
            return "anchors_on_different_contigs"
        return "anchor_missing:unknown"

    records = []
    for hid in sorted(d["hap_id"].unique()):
        sample, hap = hid.split("\t", 1)
        if hid in both.index:
            row, a = both.loc[hid], agg.loc[hid] if hid in agg.index else None
            records.append(dict(
                sample=sample, hap=hap, complete=1, contig=row["contig"],
                allele_complete=int(a is not None and a["allele_ok"]),
                n_genes=int(a["n_genes"]) if a is not None else 0,
                structure=a["structure"] if a is not None else "",
                structure_allele=a["structure_allele"] if a is not None else "",
                span_bp=int(row["hi"] - row["lo"]),
                flipped=bool(row["flipped"]), reason="spanned"))
        else:
            records.append(dict(
                sample=sample, hap=hap, complete=0, contig="",
                allele_complete=0,
                n_genes=0, structure=MISSING_VALUE,
                structure_allele=MISSING_VALUE, span_bp=0,
                flipped=False, reason=reason(hid)))
    return records


def collect_block_haplotypes(annot, blocks=None, *, sample_col="Sample",
                             hap_col="Hap", gene_col="Gene",
                             allele_col="Consensus_Allele", contig_col="Contig",
                             start_col="Spos", end_col="Epos", strand_col="Strand",
                             hla_digits=HLA_DIGITS, kir_digits=KIR_DIGITS,
                             verbose=True, assume_complete=False):
    """Per (sample, haplotype, block) gene structure from one annotation table.

    annot      : Immuannot-style table (Sample, Hap, Gene, Consensus_Allele,
                 Contig, Spos, Epos, Strand).
    blocks     : {name: {"genes": [...], "anchors": (left, right)}}; default BLOCKS.
    hap_col    : column holding the haplotype id within a sample; None = each
                 contig is its own haplotype.
    hla_digits / kir_digits : allele resolution for `structure_allele`
                 (0/2/4/6/8 HLA digits, 0/3/5/7 KIR digits; 0 = gene level).

    Returns a DataFrame with one row per haplotype per block:
      sample, hap, block, complete (1/0), allele_complete (1/0), structure, structure_allele, n_genes,
      contig, span_bp, flipped, reason
    Missing anchored blocks have complete=0 and structure "0". Anchorless
    partial observations retain their structure for inspection. allele_complete
    additionally requires usable, unambiguous calls at the requested resolution.
    """
    blocks = BLOCKS if blocks is None else blocks
    d = _normalise(annot, sample_col, hap_col, gene_col, allele_col,
                   contig_col, start_col, end_col, strand_col)
    d = collapse_copies(d, sample="sample", hap="hap", gene="gene",
                        allele="allele", contig="contig", start="spos", end="epos",
                        verbose=verbose)
    def resolve(row):
        gene = row["gene"]
        labels = {trim_allele(gene, call, hla_digits, kir_digits)
                  for call in candidate_alleles(row, "allele")}
        gene_level = (gene.startswith("C4") or
                      (kir_digits == 0 if gene.startswith("KIR") else hla_digits == 0))
        ok = len(labels) == 1 and (gene_level or gene not in labels)
        return (next(iter(labels)) if len(labels) == 1 else gene, ok)
    resolved = [resolve(r) for _, r in d.iterrows()]
    d["resolved_allele"] = [r[0] for r in resolved]
    d["allele_ok"] = [r[1] for r in resolved]
    n_haps = d["hap_id"].nunique()

    out = []
    for name, spec in blocks.items():
        anchors = spec.get("anchors")
        if anchors:
            rows = _block_rows(d, spec["genes"], tuple(anchors), hla_digits, kir_digits)
            how = f"anchors {anchors[0]}..{anchors[1]}"
        else:
            rows = _block_rows_anchorless(d, spec["genes"], hla_digits, kir_digits,
                                           assume_complete=assume_complete)
            how = "no anchors, completeness assumed" if assume_complete else "no anchors, partial observation"
        for r in rows:
            r["block"] = name
            if str(r["hap"]).strip().lower() in ("", "na", "nan", "none", "."):
                r["complete"] = r["allele_complete"] = 0
                r["reason"] = "haplotype_id_missing"
        out.extend(rows)
        if verbose:
            n_ok = sum(r["complete"] for r in rows)
            percent = 100 * n_ok / n_haps if n_haps else 0
            print(f"{name:8s} {how}: {n_ok}/{n_haps} haplotypes complete "
                  f"({percent:.1f}%)")

    cols = ["sample", "hap", "block", "complete", "allele_complete", "structure", "structure_allele",
            "n_genes", "contig", "span_bp", "flipped", "reason"]
    res = pd.DataFrame.from_records(out, columns=cols)
    res["block"] = pd.Categorical(res["block"], list(blocks), ordered=True)
    return res.sort_values(["sample", "hap", "block"]).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Reshaping                                                                    #
# --------------------------------------------------------------------------- #
def wide_structures(long, value="structure"):
    """One row per (sample, hap), one column per block ("0" where incomplete)."""
    w = long.pivot_table(index=["sample", "hap"], columns="block", values=value,
                         aggfunc="first", observed=True)
    w.columns.name = None
    return w.reset_index()


def _chrom_key(chrom):
    rest = chrom[3:] if chrom.startswith("chr") else chrom
    try:
        return (0, int(rest))
    except ValueError:
        return (1, rest)


def _safe_token(text):
    """VCF INFO value: drop the characters that split fields."""
    return "".join(ch if ch not in " \t;,=" else "_" for ch in str(text))


def encode_dosage(n, additive=True):
    """Dosage written to the VCF.

    Dominant maps counts 0, 1, 2+ to 0, 1, 1. Additive keeps the raw count,
    which is 0, 1, or 2 for a diploid sample.
    """
    n = int(n)
    if additive or n <= 1:
        return n
    return 1


def dominant_carriers(long, value="structure"):
    """How many complete haplotypes of each sample carry each structure.

    Returns {block: {structure: {sample: n_haplotypes}}}. Incomplete haplotypes
    (structure "0") are not counted. `encode_dosage` turns the count into the
    VCF dosage.
    """
    ok = long[(long["complete"] == 1)
              & long[value].ne(MISSING_VALUE)
              & long[value].astype(str).str.len().gt(0)]
    if value == "structure_allele" and "allele_complete" in ok:
        ok = ok[ok["allele_complete"].eq(1)]
    carriers = {}
    for (block, structure), grp in ok.groupby(["block", value], observed=True):
        counts = grp.groupby(grp["sample"].astype(str)).size()
        carriers.setdefault(str(block), {})[structure] = counts.to_dict()
    return carriers


def _place_structures(carriers):
    """(chrom, pos, block, structure) with a unique POS on each chromosome."""
    placed = []
    for block, structures in carriers.items():
        chrom, start = BLOCK_ANCHOR.get(block, ("chr6", 1))
        for i, structure in enumerate(sorted(structures)):
            placed.append((chrom, start + i, block, structure))
    placed.sort(key=lambda row: (_chrom_key(row[0]), row[1], row[2], row[3]))
    last = {}
    unique = []
    for chrom, pos, block, structure in placed:
        if chrom in last and pos <= last[chrom]:
            pos = last[chrom] + 1
        last[chrom] = pos
        unique.append((chrom, pos, block, structure))
    return unique


def write_structure_vcf(long, samples, path, value="structure", min_dosage=0,
                        exclude_samples=(), source="hap_block_haplotypes",
                        additive=True):
    """DS-dosage VCF of block haplotypes, one row per haplotype.

    The record layout matches allele_dosage_vcf.py. The ID is the block name
    plus a 1-based index within that block (DR_1). HI is the haplotype, with
    the same characters removed that allele_dosage_vcf.py removes from an
    allele name. Dosage is additive (0, 1, 2) by default. Dominant is explicit.
    Exact dosage requires two distinct complete haplotypes; unresolved dosage
    is '.', except that a confirmed dominant carrier is always 1.
    Samples in exclude_samples are
    dropped. A marker whose dosages sum to less than min_dosage is dropped;
    0 keeps every marker.
    """
    excluded = set(exclude_samples)
    sample_list = [s for s in samples if s not in excluded]
    if long.duplicated(["sample", "hap", "block"]).any():
        raise ValueError("expected one row per sample, haplotype, and block")
    carriers = dominant_carriers(long, value)
    callable_rows = long["complete"].eq(1)
    if value == "structure_allele" and "allele_complete" in long:
        callable_rows &= long["allele_complete"].eq(1)
    state = long.assign(_callable=callable_rows).groupby(["sample", "block"], observed=True)
    complete = {key: len(g) == 2 and g["hap"].nunique() == 2 and g["_callable"].all()
                for key, g in state}
    placed = _place_structures(carriers)
    unknown = sorted(set(carriers) - set(BLOCK_ANCHOR))
    if unknown:
        print(f"{len(unknown)} blocks have no anchor coordinate and were "
              f"placed at chr6:1: {', '.join(unknown)}", file=sys.stderr)

    contigs = list(CONTIG_LENGTH)
    for chrom, _pos, _block, _structure in placed:
        if chrom not in contigs:
            contigs.append(chrom)
    contigs.sort(key=_chrom_key)

    header = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
              "INFO", "FORMAT", *sample_list]
    meta = [
        "##fileformat=VCFv4.3",
        f"##fileDate={date.today():%Y%m%d}",
        f"##source={source}",
        "##reference=GRCh38",
        "##dosageMode=" + ("additive" if additive else "dominant"),
    ]
    for chrom in contigs:
        length = CONTIG_LENGTH.get(chrom)
        meta.append(f"##contig=<ID={chrom},length={length}>" if length
                    else f"##contig=<ID={chrom}>")
    meta.append('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">')
    meta.append('##INFO=<ID=HI,Number=1,Type=String,'
                'Description="Haplotype Information">')

    rank = {}
    n_written = 0
    with open(path, "w", buffering=1 << 20) as fp:
        fp.write("".join(line + "\n" for line in meta))
        fp.write("#" + "\t".join(header) + "\n")
        for chrom, pos, block, structure in placed:
            counts = carriers[block][structure]
            row = []
            total = 0
            for sample in sample_list:
                count = counts.get(sample, 0)
                if not additive and count:
                    dosage = 1
                elif complete.get((sample, block), False):
                    dosage = encode_dosage(count, additive)
                else:
                    dosage = None
                if dosage is not None:
                    total += dosage
                row.append("." if dosage is None else str(dosage))
            if min_dosage and total < min_dosage:
                continue
            rank[block] = rank.get(block, 0) + 1
            vid = f"{block}_{rank[block]}"
            fp.write(f"{chrom}\t{pos}\t{vid}\tA\tG\t.\t.\t"
                     f"HI={_safe_token(structure)}\tDS\t"
                     + "\t".join(row) + "\n")
            n_written += 1
    return n_written


def structure_counts(long, value="structure"):
    """How often each distinct block structure occurs across haplotypes."""
    ok = long[long["complete"] == 1]
    if value == "structure_allele" and "allele_complete" in ok:
        ok = ok[ok["allele_complete"].eq(1)]
    c = (ok.groupby(["block", value], observed=True).size()
         .rename("n_haplotypes").reset_index())
    tot = c.groupby("block", observed=True)["n_haplotypes"].transform("sum")
    c["freq"] = c["n_haplotypes"] / tot
    return c.sort_values(["block", "n_haplotypes"], ascending=[True, False]) \
            .reset_index(drop=True)


# --------------------------------------------------------------------------- #
# CLI                                                                          #
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True,
                    help="Immuannot table (.csv / .tsv, optionally gzipped)")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--sep", default=None,
                    help="field separator (default: sniffed; the Immuannot '.csv' "
                         "tables are in fact tab-separated)")
    ap.add_argument("--sample-col", default="Sample")
    ap.add_argument("--hap-col", default="Hap",
                    help="haplotype column; 'none' = each contig is a haplotype")
    ap.add_argument("--gene-col", default="Gene")
    ap.add_argument("--allele-col", default="Consensus_Allele")
    ap.add_argument("--hla-digits", type=int, default=HLA_DIGITS,
                    choices=HLA_DIGIT_CHOICES,
                    help="HLA resolution in structure_allele: "
                         "0 -> gene name, 2 -> A*01, 4 -> A*01:01, "
                         "6 -> A*01:01:01, 8 -> A*01:01:01:01 (default 2)")
    ap.add_argument("--kir-digits", type=int, default=KIR_DIGITS,
                    choices=KIR_DIGIT_CHOICES,
                    help="KIR resolution in structure_allele: "
                         "0 -> gene name, 3 -> KIR2DL1*001, "
                         "5 -> KIR2DL1*00101, 7 -> KIR2DL1*0010101 (default 3)")
    coding = ap.add_mutually_exclusive_group()
    coding.add_argument("--additive", dest="additive", action="store_true",
                        help="raw haplotype count 0, 1, 2 (default)")
    coding.add_argument("--dominant", dest="additive", action="store_false",
                        help="carrier coding 0, 1, 1")
    ap.set_defaults(additive=True)
    ap.add_argument("--assume-complete", action="store_true",
                    help="certify complete anchorless RCCX annotation on each haplotype; "
                         "multi-contig order and uncalled alleles remain unresolved")
    ap.add_argument("--min-dosage", type=int, default=0,
                    help="drop a marker whose dosages sum to less than this")
    ap.add_argument("--exclude-samples", default=None,
                    help="file of sample IDs to leave out of the VCF, one per line")
    ap.add_argument("--no-vcf", action="store_true",
                    help="do not write the dosage VCF")
    args = ap.parse_args()

    # sniff the separator: the Immuannot ".csv" tables are tab-separated
    annot = (pd.read_csv(args.annot, sep=args.sep) if args.sep
             else pd.read_csv(args.annot, sep=None, engine="python"))

    hap_col = None if str(args.hap_col).lower() == "none" else args.hap_col
    long = collect_block_haplotypes(
        annot, sample_col=args.sample_col, hap_col=hap_col, gene_col=args.gene_col,
        allele_col=args.allele_col, hla_digits=args.hla_digits,
        kir_digits=args.kir_digits, assume_complete=args.assume_complete)

    long.to_csv(f"{args.out}.block_haplotypes.tsv", sep="\t", index=False)
    wide_structures(long, "structure").to_csv(
        f"{args.out}.block_structures_wide.tsv", sep="\t", index=False)
    wide_structures(long, "structure_allele").to_csv(
        f"{args.out}.block_haplotypes_wide.tsv", sep="\t", index=False)
    structure_counts(long, "structure_allele").to_csv(
        f"{args.out}.block_structure_counts.tsv", sep="\t", index=False)
    ok = long[long["complete"] == 1]
    n_ok = len(ok)
    n_allele = int(ok["structure_allele"].astype(str).str.contains(r"\*").sum()) if n_ok else 0
    print(f"Wrote {args.out}.block_haplotypes.tsv, "
          f"{args.out}.block_structures_wide.tsv, "
          f"{args.out}.block_haplotypes_wide.tsv, "
          f"{args.out}.block_structure_counts.tsv")
    print(f"{n_allele}/{n_ok} complete block haplotypes carry an allele name "
          f"(hla-digits {args.hla_digits}, kir-digits {args.kir_digits})")
    if n_ok and n_allele == 0 and (args.hla_digits or args.kir_digits):
        print("haplotypes are gene names: the table has no usable allele calls "
              "at this resolution", file=sys.stderr)

    if not args.no_vcf:
        excluded = set()
        if args.exclude_samples:
            with open(args.exclude_samples) as fh:
                excluded = {line.strip() for line in fh if line.strip()}
        samples = sorted(long["sample"].astype(str).unique())
        path = f"{args.out}.block_haplotype.dosage.vcf"
        n = write_structure_vcf(
            long, samples, path, value="structure_allele",
            min_dosage=args.min_dosage, exclude_samples=excluded,
            source=f"hap_block_haplotypes_hla{args.hla_digits}_kir{args.kir_digits}",
            additive=args.additive)
        print(f"{n} haplotype markers -> {path}")


if __name__ == "__main__":
    main()
