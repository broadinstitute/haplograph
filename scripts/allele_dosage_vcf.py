#!/usr/bin/env python3
"""
HLA/KIR/C4 allele dosages from a per-haplotype annotation table.

Every annotated gene copy contributes one count to its allele at one HLA
resolution (KIR still nests coarse to fine), so a marker's dosage is the
number of copies the sample carries that are compatible with it:

    HLA-A*01:01:01:01  ->  HLA-A*01                 (default: 2-digit allele only)
    KIR2DL1*0020101    ->  KIR2DL1*002, KIR2DL1*00201, KIR2DL1*0020101
    C4AL               ->  C4AL              (counted by gene; allele ignored)

HLA keeps one resolution (`--hla-digits`, default 2). A finer call is trimmed to
that resolution, so HLA-A*01:01 and HLA-A*01:02 both count as HLA-A*01. The bare
gene is not a marker. KIR stays nested, coarse to fine, and likewise omits the
bare gene.

Classical HLA calls are resolved within each haplotype, preserving homozygotes.
Competing calls are usable only at resolutions where they agree.

What IS an error is counting annotation rows instead of gene copies. When one
locus yields several rows -- overlapping annotations, two candidate consensus
calls, the same contig supplied twice -- the extra rows usually differ only in
the deep fields, so they merge at low resolution and split at high resolution.
Coarse allele dosages inflate while the deep markers still look plausible.
`collapse_copies` removes that by merging overlapping annotations of one gene on
one contig into a single copy before anything is counted; genuinely duplicated
genes (C4, KIR2DL4/3DP1/3DS1) sit at non-overlapping coordinates and survive.

Haplotype-aware by design: the (sample, hap) pair is part of the locus key, so a
homozygote keeps both copies.

Additive dosage is the default. Unresolved dosage is '.' in VCF and NA in the
matrix. Exact dosage requires two observed haplotypes and usable locus calls;
variable genes additionally require a spanned block. C4 has no completeness
anchors: use --assume-complete only for independently verified complete input.
Counters hold confirmed observed copies; use dosage_matrix or dosage_value to
include uncertainty rather than interpreting a Counter's absent entry as zero.

Library use
-----------
    import pandas as pd
    from allele_dosage_vcf import collapse_copies, allele_dosages, dosage_matrix

    df = pd.read_csv("hprc.immuannot.csv", sep=None, engine="python")
    Allele_dict = allele_dosages(collapse_copies(df))   # {marker: {sample: n}}
    mat = dosage_matrix(df)                             # markers x samples frame
    write_dosage_vcf(Allele_dict, sorted(df["Sample"].unique()),
                     "out.dosage.vcf")

    duplicate_loci(df).head(20)       # how many rows the collapse would remove

CLI
---
    python allele_dosage_vcf.py --annot hprc.immuannot.csv --out hprc
writes <out>.dosage_matrix.tsv, <out>.duplicate_loci.tsv, and <out>.dosage.vcf.
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from datetime import date

import pandas as pd

# Genes whose copy number genuinely varies between haplotypes; everything else
# is expected exactly once per haplotype, so a second copy is an artefact.
# Same model as hla_kir_dosage_vcf.py's VARIABLE_REGIONS.
VARIABLE_COPY_GENES = frozenset({
    "HLA-DRB2", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DRB6", "HLA-DRB7",
    "HLA-DRB8", "HLA-DRB9", "C4A", "C4B", "C4AL", "C4BL", "C4AS", "C4BS",
    "MICA", "HLA-Y", "HLA-H", "HLA-K", "HLA-T", "HLA-U",
    "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B",
    "KIR2DP1", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5",
    "KIR3DL1", "KIR3DP1", "KIR3DS1",
})

# Classical loci on alternate contigs are representations of one haplotype.
CLASSIC_HLA_GENES = frozenset({
    "HLA-A", "HLA-B", "HLA-C",
    "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
})

HLA_DIGITS = 2            # 2, 4, 6, or 8: the one HLA resolution emitted
KIR_DIGITS = (3, 5, 7)    # KIR resolutions emitted, coarse to fine
NOVEL_TAGS = (":new", "new")   # Immuannot's novel-allele marks, end-anchored

# Markers that are not alleles. Dropped rather than written as variants.
SKIP_MARKERS = frozenset({"", "NA", "undetermined", "Undetermined", "nan", "None"})

# GRCh38 primary contig lengths (the chr11 header is kept so a coordinate
# table that uses it still validates; no default gene sits on chr11).
CONTIG_LENGTH = {
    "chr6": 170805979,
    "chr11": 135086622,
    "chr19": 58617616,
}

# Gene -> (chrom, pos, end). pos is the BED start, used as the VCF position
# label. Intervals are HLA_gene_intervals_hg38_v2.bed, plus HLA/KIR genes from
# the primary assembly that bed omits. DRB3/4, HLA-Y, and the KIR genes that
# are not on the GRCh38 reference haplotype are labels in the right interval,
# not transcription starts. C4 length forms are filled in from C4A/C4B.
_COORD_TABLE = """
chr6 26087226 26098343 HLA-HFE
chr6 29723434 29727296 HLA-F
chr6 29790954 29805501 HLA-V
chr6 29800415 29802425 HLA-P
chr6 29827825 29831021 HLA-G
chr6 29887752 29890482 HLA-H
chr6 29896654 29897786 HLA-T
chr6 29926459 29929232 HLA-K
chr6 29934101 29934286 HLA-U
chr6 29942532 29945870 HLA-A
chr6 29956596 29958570 HLA-W
chr6 29980000 29981000 HLA-Y
chr6 30006606 30009539 HLA-J
chr6 30259548 30293016 HLA-L
chr6 30351416 30351550 HLA-N
chr6 30489509 30494194 HLA-E
chr6 31268749 31272130 HLA-C
chr6 31353867 31367067 HLA-B
chr6 31382074 31382288 HLA-S
chr6 31403595 31415315 MICA
chr6 31494881 31511124 MICB
chr6 31982057 32002681 C4A
chr6 32014795 32035418 C4B
chr6 32439872 32445047 HLA-DRA
chr6 32517345 32530294 HLA-DRB5
chr6 32536000 32537000 HLA-DRB9
chr6 32542000 32543000 HLA-DRB4
chr6 32548000 32549000 HLA-DRB3
chr6 32554000 32555000 HLA-DRB2
chr6 32560000 32561000 HLA-DRB6
chr6 32566000 32567000 HLA-DRB7
chr6 32572000 32573000 HLA-DRB8
chr6 32578775 32589848 HLA-DRB1
chr6 32637406 32643684 HLA-DQA1
chr6 32659467 32668383 HLA-DQB1
chr6 32741391 32747198 HLA-DQA2
chr6 32756094 32763543 HLA-DQB2
chr6 32812763 32820466 HLA-DOB
chr6 32825415 32838739 TAP2
chr6 32845209 32853704 TAP1
chr6 32934629 32941028 HLA-DMB
chr6 32948613 32969094 HLA-DMA
chr6 33004182 33009591 HLA-DOA
chr6 33064569 33080781 HLA-DPA1
chr6 33073786 33089696 HLA-DPB1
chr6 33091753 33097295 HLA-DPA2
chr6 33112516 33129084 HLA-DPB2
chr19 54724442 54736632 KIR3DL3
chr19 54736700 54737700 KIR2DS2
chr19 54738513 54753052 KIR2DL3
chr19 54753200 54754200 KIR2DL2
chr19 54753600 54754600 KIR2DL5B
chr19 54754000 54755000 KIR2DS3
chr19 54754400 54755400 KIR2DS5
chr19 54755023 54767371 KIR2DP1
chr19 54769793 54784322 KIR2DL1
chr19 54786355 54790318 KIR3DP1
chr19 54803610 54814517 KIR2DL4
chr19 54815000 54816000 KIR2DL5A
chr19 54816468 54830778 KIR3DL1
chr19 54831200 54832200 KIR3DS1
chr19 54832676 54848569 KIR2DS4
chr19 54849000 54850000 KIR2DS1
chr19 54850443 54867207 KIR3DL2
"""


def builtin_coordinates():
    """{gene: (chrom, pos, end)} for the MHC and KIR genes this script emits."""
    coords = {}
    for line in _COORD_TABLE.strip().splitlines():
        chrom, start, end, name = line.split()
        coords[name] = (chrom, int(start), int(end))
    for form, parent in (("C4AL", "C4A"), ("C4AS", "C4A"),
                         ("C4BL", "C4B"), ("C4BS", "C4B")):
        coords[form] = coords[parent]
    return coords


GENE_COORDINATES = builtin_coordinates()


def load_coordinates(path):
    """Read a coordinate table over the builtin one.

    A line is either BED (`chrom start end gene`) or `gene chrom start end`.
    `start` is used as the VCF POS, same as the builtin table. Genes in the
    file replace the builtin entry; genes absent from the file keep theirs.
    """
    coords = builtin_coordinates()
    named = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split()
            if len(parts) < 4:
                raise ValueError(f"coordinate line needs 4 fields: {line}")
            if parts[0].lower().startswith("chr"):
                chrom, start, end, name = parts[:4]
            else:
                name, chrom, start, end = parts[:4]
            name = name.strip('"')
            coords[name] = (chrom, int(start), int(end))
            named.add(name)
    for form, parent in (("C4AL", "C4A"), ("C4AS", "C4A"),
                         ("C4BL", "C4B"), ("C4BS", "C4B")):
        if form not in named and parent in coords:
            coords[form] = coords[parent]
    return coords

# column name -> candidates tried in order, case-insensitively
COLUMN_CANDIDATES = {
    "sample": ("Sample", "sample", "research_id"),
    "hap": ("Hap", "hap", "Haplotype", "haplotype"),
    "gene": ("Gene", "gene"),
    "allele": ("Consensus_Allele", "Allele", "Alleles", "template_allele", "consensus"),
    "contig": ("Contig", "contig"),
    "start": ("Spos", "start", "Start"),
    "end": ("Epos", "end", "End"),
}


def resolve_columns(df, **overrides):
    """Map the roles above onto real column names; missing roles come back None."""
    lower = {str(c).lower(): c for c in df.columns}
    out = {}
    for role, candidates in COLUMN_CANDIDATES.items():
        if overrides.get(role):
            name = overrides[role]
            if name not in df.columns:
                raise ValueError(f"{role} column {name!r} is not in the table")
            out[role] = name
            continue
        out[role] = next((lower[c.lower()] for c in candidates
                          if c.lower() in lower), None)
    for role in ("sample", "gene"):
        if out[role] is None:
            raise ValueError(f"no {role} column found; tried "
                             f"{list(COLUMN_CANDIDATES[role])}")
    return out


# --------------------------------------------------------------------------- #
# Marker names                                                                 #
# --------------------------------------------------------------------------- #
def strip_novel(allele):
    """Drop Immuannot's novel tag, anchored at the end: ':new' or 'new'."""
    a = str(allele).strip()
    for tag in NOVEL_TAGS:
        if a.endswith(tag):
            return a[:-len(tag)]
    return a


def marker_keys(gene, allele, hla_digits=HLA_DIGITS, kir_digits=KIR_DIGITS):
    """Markers one annotated gene copy contributes.

    HLA alleles are kept at a single resolution, `hla_digits`:
        2 -> HLA-A*01, 4 -> HLA-A*01:01, 6 -> HLA-A*01:01:01,
        8 -> HLA-A*01:01:01:01.
    A deeper call is trimmed to that many fields. A call shorter than requested
    contributes nothing, so neither a 2-digit name nor the bare gene appears in
    a 4-digit matrix. KIR still emits every length in `kir_digits`, and omits
    the bare gene. C4 is counted by gene name (C4AL / C4BL / C4AS / C4BS); the
    allele adds nothing. A copy with no allele call, or a name with no '*',
    contributes nothing outside C4.
    """
    gene = str(gene)
    if gene.startswith("C4"):
        return [gene]
    if allele is None:
        return []
    a = strip_novel(allele)
    if "*" not in a:
        return []
    g, _, name = a.partition("*")
    if not name or not g:
        return []
    if g.startswith("KIR"):
        keys = [f"{g}*{name[:n]}" for n in kir_digits if len(name) >= n]
    else:
        if hla_digits not in (2, 4, 6, 8):
            raise ValueError(f"hla_digits must be 2, 4, 6, or 8, got {hla_digits}")
        n_fields = hla_digits // 2
        fields = name.split(":")
        keys = []
        if len(fields) >= n_fields and all(fields[:n_fields]):
            keys.append(f"{g}*" + ":".join(fields[:n_fields]))
    # a short KIR name yields the same prefix more than once; it is one marker
    return list(dict.fromkeys(keys))


# --------------------------------------------------------------------------- #
# One row per gene copy                                                        #
# --------------------------------------------------------------------------- #
def _locus_ids(d, grp, starts, ends):
    """Merge overlapping intervals within each group; returns a locus id per row."""
    ids, run_end, run_id = [0] * len(d), None, -1
    for i in range(len(d)):
        if i == 0 or grp[i] != grp[i - 1] or starts[i] > run_end:
            run_id += 1
            run_end = ends[i]
        else:
            run_end = max(run_end, ends[i])
        ids[i] = run_id
    return ids


def candidate_alleles(row, allele_col):
    """All competing calls retained for one physical copy, including no-call."""
    if "_candidate_alleles" in row:
        return row["_candidate_alleles"]
    value = row[allele_col] if allele_col is not None else None
    return (None if pd.isna(value) or str(value).strip().lower() in
            ("", "na", "nan", "none", ".", "undetermined") else str(value).strip(),)


def collapse_copies(df, verbose=True, cross_contig="keep", one_per_hap=None,
                    **overrides):
    """Collapse redundant annotation rows so each row is one gene copy.

    Rows for the same gene on the same (sample, haplotype, contig) whose
    intervals overlap are one locus and keep one representative: a row carrying
    an allele is preferred, then the longest annotation. All usable competing
    calls survive in _candidate_alleles for resolution-aware counting.
    Non-overlapping copies of a gene are left alone, so real duplications survive.

    cross_contig : what to do when one (sample, haplotype, gene) has copies on
        several contigs -- usually redundant/alternate contigs of the same
        region, which coordinate overlap cannot detect because the coordinates
        live in different frames.
        "keep"    (default) leave them; every contig's copy is counted.
        "largest" keep only the contig carrying the most copies of that gene,
                  as hap_structure_plot.py's fragmented="largest" does. This
                  also discards genuine duplications split across contigs, so
                  check copy_number_outliers() before turning it on.

    one_per_hap : cap a haplotype at one copy of a gene, for genes that cannot
        legitimately be duplicated. Catches what coordinates cannot: two calls
        for one locus that neither overlap nor sit on different contigs, and
        unphased tables where both of a heterozygote's alleles were written to
        both haplotypes (HLA-A hap1 *34:02 + *68:01, hap2 the same two).
        None (default) off; "auto" = every gene outside VARIABLE_COPY_GENES;
        or pass your own collection of gene names. The surviving copy is the
        one with an allele call, then the longest annotation -- so which of two
        competing alleles wins is arbitrary. Prefer fixing the input where you
        can; use this when the table is what it is.

    Needs contig and start/end columns. Without them the table cannot be
    deduplicated safely -- collapsing on (sample, gene, allele) would merge a
    homozygote's two real copies -- so the frame is returned unchanged with a
    warning.
    """
    if cross_contig not in ("keep", "largest"):
        raise ValueError("cross_contig must be 'keep' or 'largest'")
    c = resolve_columns(df, **overrides)
    needed = [c["contig"], c["start"], c["end"]]
    if any(x is None for x in needed):
        if verbose:
            print("collapse_copies: no contig/start/end columns, returning the "
                  "table unchanged (duplicate loci will inflate coarse markers)")
        return df

    key = [c["sample"], c["gene"], c["contig"]]
    if c["hap"] is not None:
        key.insert(1, c["hap"])   # a homozygote must keep both copies
    elif verbose:
        print("collapse_copies: no haplotype column; copies of one gene on one "
              "contig are merged by position only")

    d = df.sort_values(key + [c["start"], c["end"]]).copy()
    candidates = [candidate_alleles(r, c["allele"]) for _, r in d.iterrows()]
    d["_candidate_alleles"] = candidates
    d["_locus"] = _locus_ids(d, d.groupby(key, sort=False).ngroup().to_numpy(),
                             d[c["start"]].to_numpy(), d[c["end"]].to_numpy())
    d["_has"] = (d[c["allele"]].notna() if c["allele"] is not None
                 else False)
    d["_len"] = d[c["end"]] - d[c["start"]]
    calls = d.groupby("_locus")["_candidate_alleles"].agg(
        lambda rows: tuple(dict.fromkeys(a for calls in rows for a in calls)))
    # Missing annotations add no alternative when a usable call exists at the
    # same locus. Two different usable calls remain competing candidates.
    calls = calls.map(lambda values: tuple(v for v in values if v is not None) or (None,))
    d["_candidate_alleles"] = d["_locus"].map(calls)
    out = (d.sort_values(["_locus", "_has", "_len"])
             .drop_duplicates("_locus", keep="last")
             .drop(columns=["_locus", "_has", "_len"])
             .sort_index())
    n_overlap = len(df) - len(out)

    n_cross = 0
    if cross_contig == "largest":
        gkey = [k for k in key if k != c["contig"]]
        k = out.groupby(gkey + [c["contig"]]).size().rename("k").reset_index()
        best = k.sort_values(["k", c["contig"]]).drop_duplicates(gkey, keep="last")
        allowed = set(map(tuple, best[gkey + [c["contig"]]].to_numpy()))
        mask = [tuple(r) in allowed
                for r in out[gkey + [c["contig"]]].to_numpy()]
        n_cross = len(out) - sum(mask)
        out = out[mask]

    n_capped = 0
    if one_per_hap is not None:
        gkey = [k for k in key if k != c["contig"]]   # sample [, hap], gene
        o = out.copy()
        o["_has"] = (o[c["allele"]].notna() if c["allele"] is not None else False)
        o["_len"] = o[c["end"]] - o[c["start"]]
        gene_str = o[c["gene"]].astype(str)
        exempt = (gene_str.isin(VARIABLE_COPY_GENES) if one_per_hap == "auto"
                  else ~gene_str.isin(set(one_per_hap)))
        keep = (o[~exempt].sort_values(["_has", "_len"])
                .groupby(gkey, sort=False).tail(1))
        n_capped = int((~exempt).sum()) - len(keep)
        out = out[out.index.isin(set(o.index[exempt]) | set(keep.index))]

    if verbose and (n_overlap or n_cross or n_capped):
        extra = f", {n_cross} on redundant contigs" if n_cross else ""
        extra += (f", {n_capped} beyond one copy per haplotype"
                  if n_capped else "")
        print(f"collapse_copies: {n_overlap} overlapping annotation rows "
              f"removed{extra} ({len(df)} -> {len(out)} gene copies)")
    return out


def duplicate_loci(df, **overrides):
    """Rows per (sample, hap, gene) above one, worst first -- the collapse preview.

    Overlapping entries are redundant annotations; non-overlapping ones are real
    duplications (C4, KIR2DL4/3DP1/3DS1) and are left alone by collapse_copies.
    """
    c = resolve_columns(df, **overrides)
    key = [c[k] for k in ("sample", "hap", "gene") if c[k] is not None]
    n = df.groupby(key).size().rename("n_rows")
    kept = collapse_copies(df, verbose=False, **overrides).groupby(key).size() \
        .rename("n_copies")
    out = pd.concat([n, kept], axis=1).fillna(0).astype(int)
    out["n_redundant"] = out["n_rows"] - out["n_copies"]
    return (out[out["n_rows"] > 1]
            .sort_values(["n_redundant", "n_rows"], ascending=False)
            .reset_index())


def copy_number_outliers(df, expected=2, genes=None, collapse=True, **overrides):
    """(sample, gene) pairs carrying more copies than a diploid can, worst first.

    Any gene that is not copy-variable should top out at `expected` copies per
    sample, so HLA-A at 4 means four rows survived for one sample. The columns
    say which explanation applies:

      n_haps    > 2  the sample has more than two haplotypes in the table
      n_contigs > n_haps  copies sit on several contigs of one haplotype
                          (redundant/alternate contigs -- cross_contig="largest")
      otherwise       several rows at one locus that collapse_copies did not
                      merge, e.g. annotations that abut rather than overlap

    collapse=True (default) measures what is left AFTER collapse_copies, which
    is what actually reaches the dosages.
    """
    c = resolve_columns(df, **overrides)
    if collapse:
        df = collapse_copies(df, verbose=False, **overrides)
    if genes is not None:
        df = df[df[c["gene"]].astype(str).isin(set(genes))]
    g = df.groupby([c["sample"], c["gene"]])
    cols = {"n_copies": g.size()}
    if c["hap"] is not None:
        cols["n_haps"] = g[c["hap"]].nunique()
    if c["contig"] is not None:
        cols["n_contigs"] = g[c["contig"]].nunique()
    if c["allele"] is not None:
        cols["alleles"] = g[c["allele"]].apply(
            lambda s: ",".join(sorted(set(s.dropna().astype(str)))))
    out = pd.DataFrame(cols)
    return (out[out["n_copies"] > expected]
            .sort_values("n_copies", ascending=False).reset_index())


# --------------------------------------------------------------------------- #
# Counting                                                                     #
# --------------------------------------------------------------------------- #
class AlleleDosages(defaultdict):
    """Confirmed copy counters plus samples whose gene dosage is unresolved."""
    def __init__(self):
        super().__init__(Counter)
        self.uncertain = defaultdict(set)
        self.uncertain_resolution = defaultdict(set)
        self.samples = []


def marker_resolution(marker):
    gene, _, allele = marker.partition("*")
    return len(allele) if gene.startswith("KIR") else len(allele.split(":")) if allele else 0


def dosage_value(dosages, marker, sample, additive=True):
    """Exact dosage or None; a confirmed carrier is known in dominant mode."""
    if hasattr(dosages, "samples") and sample not in dosages.samples:
        return None
    count = int(dosages.get(marker, {}).get(sample, 0))
    if not additive and count:
        return 1
    if sample in getattr(dosages, "uncertain", {}).get(_marker_gene(marker), set()):
        return None
    resolution_key = (_marker_gene(marker), marker_resolution(marker))
    if sample in getattr(dosages, "uncertain_resolution", {}).get(resolution_key, set()):
        return None
    return count if additive else int(count > 0)


def allele_dosages(df, hla_digits=HLA_DIGITS, kir_digits=KIR_DIGITS,
                   collapse=True, verbose=True, assume_complete=False, **overrides):
    """Confirmed observed counts with completeness and resolution uncertainty.

    Collapse physical copies by default. Classical HLA copy counts must agree
    across alternate contig representations of each haplotype. Distinct loci
    within a contig retain their multiplicity. No sample-wide set operation is
    used. A diploid sample needs exactly two nonmissing hap IDs.
    assume_complete certifies absence/completeness, not unknown allele identity.
    """
    c = resolve_columns(df, **overrides)
    d = df.copy()
    for role in ("sample", "hap", "gene", "contig"):
        if c[role] is not None:
            d[c[role]] = d[c[role]].map(lambda v: str(v).strip() if pd.notna(v) else "")
    d[c["gene"]] = d[c["gene"]].str.split("*").str[0]
    if collapse:
        d = collapse_copies(d, verbose=verbose, **overrides)
    out = AlleleDosages()
    out.samples = sorted(d[c["sample"]].unique())
    hap_col = c["hap"]
    invalid_hap = {"", "nan", "None", "NA", "."}
    haps = {s: set(g[hap_col]) - invalid_hap if hap_col else set()
            for s, g in d.groupby(c["sample"])}
    # Variable-gene absence is known only when the corresponding block spans.
    coverage, covered_contig, gene_block = {}, {}, {}
    if all(c[k] is not None for k in ("hap", "contig", "start", "end")) and len(d):
        from hap_block_haplotypes import BLOCKS, collect_block_haplotypes
        gene_block = {g: name for name, spec in BLOCKS.items() for g in spec["genes"]}
        blocks = collect_block_haplotypes(
            d, sample_col=c["sample"], hap_col=hap_col, gene_col=c["gene"],
            allele_col=c["allele"], contig_col=c["contig"], start_col=c["start"],
            end_col=c["end"], assume_complete=assume_complete, verbose=False)
        coverage = {(r.sample, r.hap, str(r.block)): bool(r.complete)
                    for r in blocks.itertuples()}
        covered_contig = {(r.sample, r.hap, str(r.block)): r.contig
                          for r in blocks.itertuples()}
    group_cols = [c["sample"], c["gene"]] + ([hap_col] if hap_col else [])
    observed = set()
    for group, frame in d.groupby(group_cols, dropna=False):
        sample, gene = group[:2]
        hap = group[2] if hap_col else ""
        observed.add((sample, hap, gene))
        # A spanned contig does not certify additional representations of the
        # same variable gene on alternate contigs as extra biological copies.
        if (not assume_complete and c["contig"] is not None and
                (gene in VARIABLE_COPY_GENES or gene.startswith(("KIR", "C4")))):
            chosen = covered_contig.get((sample, hap, gene_block.get(gene)))
            if chosen is not None and set(frame[c["contig"]]) != {chosen}:
                out.uncertain[gene].add(sample)
        representations = ([g for _, g in frame.groupby(c["contig"])]
                           if gene in CLASSIC_HLA_GENES and hap_col and c["contig"] else [frame])
        resolutions = ((0,) if gene.startswith("C4") else kir_digits
                       if gene.startswith("KIR") else (hla_digits // 2,))
        rep_counts = []
        for representation in representations:
            counts = Counter()
            for _, locus in representation.iterrows():
                alternatives = [set(marker_keys(gene, v, hla_digits, kir_digits))
                                for v in candidate_alleles(locus, c["allele"])]
                possible = set.union(*alternatives)
                certain = set.intersection(*alternatives)
                for marker in possible:
                    out[marker]  # retain ambiguous-only markers
                counts.update(certain)
                for resolution in resolutions:
                    choices = {tuple(sorted(m for m in markers if marker_resolution(m) == resolution))
                               for markers in alternatives}
                    if len(choices) != 1 or not next(iter(choices)):
                        out.uncertain_resolution[(gene, resolution)].add(sample)
            rep_counts.append(counts)
        markers = set().union(*(set(counts) for counts in rep_counts))
        for marker in markers:
            out[marker][sample] += min(counts[marker] for counts in rep_counts)
        for resolution in resolutions:
            signatures = {tuple(sorted((m, n) for m, n in counts.items()
                                       if marker_resolution(m) == resolution))
                          for counts in rep_counts}
            if len(signatures) > 1:
                out.uncertain_resolution[(gene, resolution)].add(sample)
    for gene in {_marker_gene(m) for m in out}:
        variable = gene in VARIABLE_COPY_GENES or gene.startswith(("KIR", "C4"))
        for sample in out.samples:
            if len(haps[sample]) != 2:
                out.uncertain[gene].add(sample)
            for hap in haps[sample]:
                if variable:
                    if not assume_complete and not coverage.get((sample, hap, gene_block.get(gene)), False):
                        out.uncertain[gene].add(sample)
                elif (sample, hap, gene) not in observed:
                    out.uncertain[gene].add(sample)
    return out


def dosage_matrix(df=None, samples=None, dosages=None, **kwargs):
    """Markers x samples nullable integers; unresolved dosage stays missing."""
    d = allele_dosages(df, **kwargs) if dosages is None else dosages
    samples = list(samples) if samples is not None else getattr(
        d, "samples", sorted({s for counts in d.values() for s in counts}))
    mat = pd.DataFrame([[dosage_value(d, m, s) for s in samples] for m in d],
                       index=list(d), columns=samples, dtype="Int64")
    mat = mat.sort_index(axis=0).sort_index(axis=1)
    mat.index.name = "marker"
    return mat


# --------------------------------------------------------------------------- #
# Dosage VCF                                                                   #
# --------------------------------------------------------------------------- #
def _chrom_key(chrom):
    """chr6 < chr11 < chr19, which plain string order does not give."""
    rest = chrom[3:] if chrom.startswith("chr") else chrom
    try:
        return (0, int(rest))
    except ValueError:
        return (1, rest)


def _marker_gene(marker):
    """HLA-A*01 -> HLA-A; C4AL stays C4AL."""
    return marker.split("*", 1)[0]


def keeps_copy_count(gene):
    """True for genes whose dosage is the copy count, never the dominant code.

    KIR, C4 (C4A/C4B and the long/short forms), TAP1/TAP2, and MICA vary in
    copy number between haplotypes.
    """
    g = str(gene)
    return g.startswith(("KIR", "C4", "TAP")) or g == "MICA"


def _safe_token(text):
    """VCF ID / INFO value: drop the characters that split fields."""
    return "".join(ch if ch not in " \t;,=" else "_" for ch in text)


def place_markers(markers, coordinates):
    """(chrom, pos, marker, gene) in VCF order, one unique position each.

    Every allele of a gene starts at that gene's coordinate and then takes
    the next free base, so two alleles never share a POS. Markers with no
    coordinate are returned separately and are not written.
    """
    grouped = defaultdict(list)
    missing = []
    for marker in markers:
        if marker is None or marker in SKIP_MARKERS:
            continue
        gene = _marker_gene(marker)
        if gene not in coordinates:
            missing.append(marker)
            continue
        grouped[gene].append(marker)

    placed = []
    for gene, names in grouped.items():
        chrom, start, _end = coordinates[gene]
        for i, name in enumerate(sorted(names)):
            placed.append((chrom, start + i, name, gene))
    placed.sort(key=lambda row: (_chrom_key(row[0]), row[1], row[2]))

    last = {}
    unique = []
    for chrom, pos, name, gene in placed:
        if chrom in last and pos <= last[chrom]:
            pos = last[chrom] + 1
        last[chrom] = pos
        unique.append((chrom, pos, name, gene))
    return unique, missing


def write_dosage_vcf(allele_dict, samples, path, coordinates=None,
                     mode="additive", min_dosage=0, exclude_samples=(),
                     source="HLA_KIR_Allele_digit2_VCF"):
    """Write a DS-dosage VCF, one row per marker, one column per sample.

    Default "additive" keeps copy counts. Unresolved dosage is '.'.
    mode "dominant" maps copy counts 0, 1, 2+ to dosages 0, 1, 1. KIR, C4,
    TAP, and MICA are copy-number variable, so they keep the raw count in
    either mode. "additive" keeps the raw count for every marker. Samples in
    exclude_samples are dropped. A marker whose dosages sum to less than
    min_dosage is dropped; 0 keeps every marker.
    """
    if mode not in ("dominant", "additive"):
        raise ValueError("mode must be 'dominant' or 'additive'")
    if coordinates is None:
        coordinates = GENE_COORDINATES
    excluded = set(exclude_samples)
    sample_list = [s for s in samples if s not in excluded]
    placed, missing = place_markers(allele_dict, coordinates)
    if missing:
        print(f"{len(missing)} markers have no coordinate and were skipped: "
              + ", ".join(sorted(set(missing))[:8]), file=sys.stderr)

    contigs = list(CONTIG_LENGTH)
    for chrom, _pos, _name, _gene in placed:
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
        "##dosageMode=" + mode,
    ]
    for chrom in contigs:
        length = CONTIG_LENGTH.get(chrom)
        meta.append(f"##contig=<ID={chrom},length={length}>" if length
                    else f"##contig=<ID={chrom}>")
    meta.append('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">')
    meta.append('##INFO=<ID=HI,Number=1,Type=String,'
                'Description="Allele Digit 2 Information">')

    n_written = 0
    with open(path, "w", buffering=1 << 20) as fp:
        fp.write("".join(line + "\n" for line in meta))
        fp.write("#" + "\t".join(header) + "\n")
        for chrom, pos, marker, gene in placed:
            # KIR, C4, TAP, and MICA vary in copy number, so a second copy is
            # a real dosage even when every other gene is coded dominant.
            raw = keeps_copy_count(gene) or mode == "additive"
            row = []
            total = 0
            for sample in sample_list:
                dosage = dosage_value(allele_dict, marker, sample, additive=raw)
                if dosage is not None:
                    total += dosage
                row.append("." if dosage is None else str(dosage))
            if min_dosage and total < min_dosage:
                continue
            token = _safe_token(marker)
            fp.write(f"{chrom}\t{pos}\t{token}\tA\tG\t.\t.\tHI={token}\tDS\t"
                     + "\t".join(row) + "\n")
            n_written += 1
    return n_written


# --------------------------------------------------------------------------- #
# CLI                                                                          #
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot", required=True,
                    help="annotation table (.csv / .tsv, optionally gzipped); "
                         "the separator is sniffed")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--no-collapse", action="store_true",
                    help="count annotation rows as-is (reproduces the inflated "
                         "gene / 1-field dosages; for comparison only)")
    ap.add_argument("--hla-digits", type=int, default=HLA_DIGITS, choices=(2, 4, 6, 8),
                    help="HLA resolution to keep: 2 -> A*01, 4 -> A*01:01, "
                         "6 -> A*01:01:01, 8 -> A*01:01:01:01")
    ap.add_argument("--kir-digits", type=int, nargs="+", default=list(KIR_DIGITS))
    coding = ap.add_mutually_exclusive_group()
    coding.add_argument("--additive", action="store_true", help="raw copy dosage (default)")
    coding.add_argument("--dominant", action="store_true",
                        help="carrier coding; KIR, C4, TAP, and MICA retain raw dosage")
    ap.add_argument("--assume-complete", action="store_true",
                    help="certify complete gene annotation on each observed haplotype; "
                         "does not resolve uncalled or conflicting alleles")
    ap.add_argument("--mode", choices=("dominant", "additive"), default=None,
                    help="same choice as --additive: dominant (0, 1, 1) or "
                         "additive (0, 1, 2). KIR, C4, TAP, and MICA stay raw")
    ap.add_argument("--min-dosage", type=int, default=0,
                    help="drop a marker whose dosages sum to less than this")
    ap.add_argument("--exclude-samples", default=None,
                    help="file of sample IDs to leave out of the VCF, one per line")
    ap.add_argument("--coords", default=None,
                    help="BED (chrom start end gene) or 'gene chrom start end'; "
                         "replaces builtin coordinates for genes it names")
    ap.add_argument("--no-vcf", action="store_true",
                    help="do not write <out>.dosage.vcf")
    for role in COLUMN_CANDIDATES:
        ap.add_argument(f"--col-{role}", default=None,
                        help=f"name of the {role} column")
    args = ap.parse_args()
    if (args.additive and args.mode == "dominant") or (args.dominant and args.mode == "additive"):
        ap.error("dosage mode flags disagree")

    overrides = {role: getattr(args, f"col_{role}") for role in COLUMN_CANDIDATES}
    df = pd.read_csv(args.annot, sep=None, engine="python")

    dup = duplicate_loci(df, **overrides)
    dup.to_csv(f"{args.out}.duplicate_loci.tsv", sep="\t", index=False)
    print(f"{len(dup)} sample x hap x gene loci have more than one row; "
          f"{int(dup['n_redundant'].sum()) if len(dup) else 0} rows are redundant")

    dosages = allele_dosages(df, collapse=not args.no_collapse,
                             hla_digits=args.hla_digits,
                             kir_digits=tuple(args.kir_digits),
                             assume_complete=args.assume_complete, **overrides)
    mat = dosage_matrix(dosages=dosages)
    mat.to_csv(f"{args.out}.dosage_matrix.tsv", sep="\t", na_rep="NA")
    print(f"{mat.shape[0]} markers x {mat.shape[1]} samples -> "
          f"{args.out}.dosage_matrix.tsv")

    if not args.no_vcf:
        all_samples = dosages.samples
        excluded = set()
        if args.exclude_samples:
            with open(args.exclude_samples) as fh:
                excluded = {line.strip() for line in fh if line.strip()}
        sample_list = [s for s in all_samples if s not in excluded]
        print(f"{len(sample_list)} samples in the VCF, {len(all_samples)} in "
              f"the table, {len(all_samples) - len(sample_list)} excluded")
        coordinates = (load_coordinates(args.coords) if args.coords
                       else GENE_COORDINATES)
        mode = "dominant" if args.dominant or args.mode == "dominant" else "additive"
        n_markers = write_dosage_vcf(
            dosages, sample_list, f"{args.out}.dosage.vcf",
            coordinates=coordinates, mode=mode,
            min_dosage=args.min_dosage, exclude_samples=excluded,
            source=f"HLA_KIR_Allele_digit{args.hla_digits}_VCF")
        print(f"{n_markers} markers -> {args.out}.dosage.vcf")


if __name__ == "__main__":
    main()
