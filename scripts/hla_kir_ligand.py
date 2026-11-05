"""
hla_kir_ligand — derive KIR ligand groups from the IPD-IMGT/HLA protein alignments.

Rather than hand-copying a lookup table, this reads the official alignments and
reads the motif off the sequence itself, so the assignment is reproducible and
carries the IPD-IMGT/HLA release it came from.

The rules it applies, all on mature-protein numbering:

  HLA-C  residue 80     Asn -> C1   (ligand for KIR2DL2 / KIR2DL3 / KIR2DS2)
                        Lys -> C2   (ligand for KIR2DL1 / KIR2DS1)
  HLA-B  residue 83     Arg -> Bw4, Gly -> Bw6. Bw6 is not a KIR ligand.
         residue 80     within Bw4: Ile -> Bw4-80I (strong KIR3DL1 ligand, and the
                        putative KIR3DS1 ligand), Thr -> Bw4-80T (weak).
  HLA-A  residue 83     Arg -> Bw4, subdivided at 80 as above (A*23/24/32 are the
                        common Bw4 carriers).
  HLA-A  A*03, A*11     ligand for KIR3DL2, assigned by allele group, since that
                        specificity is not a single-residue motif.
  HLA-B  B*46, B*73     carry the C1 epitope by gene conversion from HLA-C, so
                        they are KIR2DL2/2DL3 ligands despite being Bw6 at 83.
                        Confirmed on sequence: both read Val76 + Asn80, where
                        every other common HLA-B reads Glu76.

An allele whose motif residues are not sequenced (the alignments mark those '*')
is reported as "unassigned" rather than guessed at. An allele with a null
expression suffix (B*15:01:01:02N) is reported as "null": the motif is in the
sequence but the molecule never reaches the cell surface, so it presents no
ligand to KIR.

Usage
-----
    from hla_kir_ligand import build_ligand_table, assign_ligands

    tbl = build_ligand_table()              # downloads + parses, returns DataFrame
    tbl.to_csv("kir_ligand_table.tsv", sep="\\t", index=False)

    groups = assign_ligands(dosage_df.index, tbl)   # allele -> ligand group
"""

from __future__ import annotations

import re
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, Optional

import pandas as pd

__all__ = [
    "build_ligand_table",
    "assign_ligands",
    "parse_alignment",
    "KIR_LIGAND_RECEPTORS",
    "IMGT_BASE",
]

IMGT_BASE = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments"

# Which KIR genes each ligand group feeds. Inhibitory first, then activating.
KIR_LIGAND_RECEPTORS: Dict[str, list] = {
    "C1": ["KIR2DL2", "KIR2DL3", "KIR2DS2"],
    "C2": ["KIR2DL1", "KIR2DS1"],
    "Bw4-80I": ["KIR3DL1", "KIR3DS1"],
    "Bw4-80T": ["KIR3DL1"],
    "A3/A11": ["KIR3DL2"],
}


# ------------------------------------------------------------------- download


def _fetch(locus: str, cache_dir: Optional[Path]) -> str:
    name = f"{locus}_prot.txt"
    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        path = cache_dir / name
        if path.exists() and path.stat().st_size > 1000:
            return path.read_text(errors="replace")
    url = f"{IMGT_BASE}/{name}"
    with urllib.request.urlopen(url, timeout=180) as r:
        text = r.read().decode("utf-8", errors="replace")
    if cache_dir is not None:
        (cache_dir / name).write_text(text)
    return text


# --------------------------------------------------------------------- parser


def parse_alignment(text: str) -> tuple:
    """Parse an IPD-IMGT/HLA *_prot.txt alignment.

    Returns (sequences, position_index, version) where `sequences` maps allele ->
    aligned residue string with '-' already resolved against the reference, and
    `position_index` maps mature-protein position -> offset into that string.

    Numbering is walked along the REFERENCE allele (the first in the file),
    counting only real residues so that alignment indels ('.') do not consume a
    position, and stepping -1 -> 1 with no zero.
    """
    version = ""
    m = re.search(r"^# version:\s*(.+)$", text, re.M)
    if m:
        version = m.group(1).strip()

    lines = text.splitlines()
    n = len(lines)
    starts = [i for i, l in enumerate(lines) if l.lstrip().startswith("Prot")]
    if not starts:
        raise ValueError("no 'Prot' blocks found; is this an IPD-IMGT/HLA alignment?")

    seqs: Dict[str, list] = {}
    reference: Optional[str] = None
    pos1_offset: Optional[int] = None    # index of mature residue 1 in the ref string
    ref_len_so_far = 0

    row_re = re.compile(r"^\s(\S+)(?:\s|$)")

    for b, prot_i in enumerate(starts):
        end = starts[b + 1] if b + 1 < len(starts) else n
        rows = [k for k in range(prot_i + 2, end) if row_re.match(lines[k])]
        if not rows:
            continue

        # Label columns on the Prot line. These are LABEL anchors and sit to the
        # LEFT of where residues actually start, so they cannot be used to slice.
        labels = {m.group(): m.start() for m in re.finditer(r"-?\d+", lines[prot_i])}

        # the block's first populated row fixes the residue column
        ref_line = next(
            (lines[k] for k in rows if len(lines[k].split()) > 1), lines[rows[0]]
        )
        name = ref_line.split()[0]
        seq_start = len(ref_line) - len(
            ref_line[ref_line.index(name) + len(name):].lstrip()
        )
        ref_seg = ref_line[seq_start:].replace(" ", "")
        seg_len = len(ref_seg)

        if pos1_offset is None and "1" in labels:
            pos1_offset = ref_len_so_far + len(
                ref_line[seq_start: labels["1"]].replace(" ", "")
            )

        seen = set()
        for k in rows:
            allele = row_re.match(lines[k]).group(1)
            # A row may carry no sequence at all in this block (null alleles such
            # as B*07:67N). Pad it rather than treating it as the end of the block.
            seg = lines[k][seq_start:].replace(" ", "") if len(lines[k]) > seq_start else ""
            seg = seg.ljust(seg_len, "*")[:seg_len]
            if allele not in seqs:
                seqs[allele] = ["*" * ref_len_so_far] if ref_len_so_far else []
                if reference is None:
                    reference = allele
            seqs[allele].append(seg)
            seen.add(allele)

        for allele in seqs:                      # absent from this block entirely
            if allele not in seen:
                seqs[allele].append("*" * seg_len)
        ref_len_so_far += seg_len

    if reference is None:
        raise ValueError("no sequences parsed from alignment")
    if pos1_offset is None:
        raise ValueError("could not locate the mature position-1 anchor")

    joined = {a: "".join(parts) for a, parts in seqs.items()}
    ref = joined[reference]

    # resolve '-' against the reference
    resolved: Dict[str, str] = {}
    for allele, s in joined.items():
        s = s.ljust(len(ref), "*")[: len(ref)]
        resolved[allele] = "".join(
            ref[k] if ch == "-" else ch for k, ch in enumerate(s)
        )

    # Walk outward from the position-1 anchor, so numbering never depends on how
    # much leader the alignment happens to show. Indels ('.') consume no position.
    position_index: Dict[int, int] = {}
    pos = 1
    for k in range(pos1_offset, len(ref)):
        if ref[k] == ".":
            continue
        position_index[pos] = k
        pos += 1
    pos = -1
    for k in range(pos1_offset - 1, -1, -1):
        if ref[k] == ".":
            continue
        position_index[pos] = k
        pos -= 1
    return resolved, position_index, version


def _residues(seq: str, index: Dict[int, int], positions: Iterable[int]) -> Dict[int, str]:
    out = {}
    for p in positions:
        k = index.get(p)
        out[p] = seq[k] if (k is not None and k < len(seq)) else "?"
    return out


# ---------------------------------------------------------------- assignment


def _classify(locus: str, res: Dict[int, str], allele: str) -> tuple:
    """Return (ligand_group, note)."""
    r76, r80, r83 = res.get(76, "?"), res.get(80, "?"), res.get(83, "?")

    # A null allele is not transcribed into a surface molecule, so whatever the
    # motif says, there is no epitope for KIR to see. Counting one as a ligand
    # carrier inflates every ligand frequency and every missing-ligand count.
    if re.search(r"N$", allele):
        return "null", "null expression suffix; not presented"

    if locus == "C":
        if r80 == "N":
            return "C1", ""
        if r80 == "K":
            return "C2", ""
        if r80 in "*?X":
            return "unassigned", "residue 80 not sequenced"
        return "unassigned", f"unexpected residue 80 = {r80}"

    if locus in {"B", "A"}:
        # A*03 / A*11 are the KIR3DL2 ligands, by allele group not by motif
        if locus == "A" and re.match(r"^A\*(03|11):", allele):
            return "A3/A11", ""
        if r83 == "R":
            if r80 == "I":
                return "Bw4-80I", ""
            if r80 == "T":
                return "Bw4-80T", ""
            if r80 in "*?X":
                return "unassigned", "Bw4 but residue 80 not sequenced"
            return "Bw4-80T", f"unexpected residue 80 = {r80}, treated as 80T"
        if r83 == "G":
            if locus == "B":
                # B*46:01 / B*73:01 took the C1 motif from HLA-C by gene
                # conversion and are genuine KIR2DL2/2DL3 ligands. Asn80 alone
                # does NOT identify them -- B*07:02 is Ser77+Asn80 too and binds
                # nothing. Val76 is the discriminator: HLA-C and these two read
                # Val76, essentially every other HLA-B reads Glu76.
                if r76 == "V" and r80 == "N":
                    if re.match(r"^B\*(46|73):", allele):
                        return "C1", "C1 epitope on HLA-B (gene conversion)"
                    return "Bw6", (
                        "Val76+Asn80 C1-like motif outside B*46/B*73; "
                        "not an established KIR2DL2/3 ligand"
                    )
                return "Bw6", ""
            return "none", ""
        if r83 in "*?X":
            return "unassigned", "residue 83 not sequenced"
        return "unassigned", f"unexpected residue 83 = {r83}"

    return "unassigned", "unsupported locus"


def build_ligand_table(
    loci: Iterable[str] = ("A", "B", "C"),
    cache_dir: Optional[str] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """Download and parse the IPD-IMGT/HLA alignments into a ligand table.

    Returns one row per allele with its motif residues, the assigned ligand group,
    the KIR genes that group feeds, and the IPD-IMGT/HLA version it came from.
    """
    rows = []
    version = ""
    for locus in loci:
        text = _fetch(locus, Path(cache_dir) if cache_dir else None)
        seqs, index, version = parse_alignment(text)
        if verbose:
            print(f"[info] HLA-{locus}: {len(seqs)} alleles  ({version})")
        for allele, seq in seqs.items():
            res = _residues(seq, index, (76, 77, 80, 81, 82, 83))
            group, note = _classify(locus, res, allele)
            rows.append(
                {
                    "allele": allele,
                    "locus": locus,
                    "motif_77_83": "".join(res[p] for p in (77, 80, 81, 82, 83)),
                    "res76": res[76], "res77": res[77], "res80": res[80],
                    "res81": res[81], "res82": res[82], "res83": res[83],
                    "ligand_group": group,
                    "kir_receptors": ",".join(KIR_LIGAND_RECEPTORS.get(group, [])),
                    "note": note,
                    "ipd_version": version,
                }
            )
    tbl = pd.DataFrame(rows)
    if verbose:
        n_un = int((tbl.ligand_group == "unassigned").sum())
        print(
            f"[info] {len(tbl)} alleles; "
            + ", ".join(
                f"{g}={n}" for g, n in tbl.ligand_group.value_counts().items()
            )
        )
        if n_un:
            print(f"[warn] {n_un} unassigned (motif residues not sequenced)")
    return tbl


# ------------------------------------------------------------------ self-test

# Alleles whose ligand group is well established, used to check that the
# alignment parse still lines up after an IPD-IMGT/HLA release.
# NOTE A*25:01 is Bw4-80I: it reads ESLRIALR, Ile at 80, same as A*32:01.
GROUND_TRUTH: Dict[str, str] = {
    "C*01:02": "C1", "C*03:04": "C1", "C*07:01": "C1", "C*07:02": "C1",
    "C*08:02": "C1", "C*12:03": "C1", "C*14:02": "C1", "C*16:01": "C1",
    "C*02:02": "C2", "C*04:01": "C2", "C*05:01": "C2", "C*06:02": "C2",
    "C*15:02": "C2", "C*17:01": "C2",
    "B*57:01": "Bw4-80I", "B*58:01": "Bw4-80I", "B*51:01": "Bw4-80I",
    "B*52:01": "Bw4-80I", "B*53:01": "Bw4-80I", "B*38:01": "Bw4-80I",
    "B*44:02": "Bw4-80T", "B*44:03": "Bw4-80T", "B*27:05": "Bw4-80T",
    "B*13:02": "Bw4-80T", "B*37:01": "Bw4-80T",
    "B*07:02": "Bw6", "B*08:01": "Bw6", "B*35:01": "Bw6",
    "B*40:01": "Bw6", "B*15:01": "Bw6",
    # Bw6 at residue 83, but C1 by the Val76+Asn80 motif they picked up from
    # HLA-C; the two HLA-B alleles that really do feed KIR2DL2/2DL3.
    "B*46:01": "C1", "B*73:01": "C1",
    "A*23:01": "Bw4-80I", "A*24:02": "Bw4-80I", "A*32:01": "Bw4-80I",
    "A*25:01": "Bw4-80I",
    "A*03:01": "A3/A11", "A*11:01": "A3/A11",
    "A*02:01": "none", "A*01:01": "none",
    # null alleles carry the motif but never reach the surface
    "A*03:01:01:02N": "null", "B*15:01:01:02N": "null",
}


def validate_table(table: pd.DataFrame, verbose: bool = True) -> pd.DataFrame:
    """Check the derived table against alleles whose group is well established.

    Run this after an IPD-IMGT/HLA release: if the alignment format shifts, the
    residue numbering silently moves and every assignment goes with it. Returns
    the rows that disagree (empty is a pass).
    """
    idx = table.set_index("allele")
    bad = []
    for allele, expect in GROUND_TRUTH.items():
        hit = idx[(idx.index == allele) | idx.index.str.startswith(allele + ":")]
        # A two-field name like A*03:01 also prefix-matches A*03:01:01:02N, which
        # is now "null"; take the expressed representative unless the expectation
        # was written against a suffixed name in the first place.
        if len(hit) > 1 and not re.search(r"[A-Z]$", allele):
            expressed = hit[~hit.index.str.contains(r"[A-Z]$", regex=True)]
            if not expressed.empty:
                hit = expressed
        got = hit["ligand_group"].iloc[0] if len(hit) else "NOT FOUND"
        motif = hit["motif_77_83"].iloc[0] if len(hit) else ""
        if got != expect:
            bad.append({"allele": allele, "expected": expect, "got": got, "motif": motif})
    out = pd.DataFrame(bad)
    if verbose:
        n = len(GROUND_TRUTH) - len(out)
        print(f"[{'ok' if out.empty else 'FAIL'}] ligand table: {n}/{len(GROUND_TRUTH)} known alleles correct")
        for r in bad:
            print(f"   {r['allele']}: expected {r['expected']}, got {r['got']} (motif {r['motif']})")
    return out


# ------------------------------------------------- mapping onto typed alleles


def _keys(allele: str):
    """Progressively shorter forms of an allele name, for fallback matching."""
    a = str(allele).strip()
    a = re.sub(r"[A-Z]$", "", a)                       # trailing N/L/S/Q expression tag
    yield a
    if "*" in a:
        locus, fields = a.split("*", 1)
        parts = fields.split(":")
        for n in range(len(parts) - 1, 0, -1):
            yield f"{locus}*{':'.join(parts[:n])}"


def assign_ligands(
    alleles: Iterable[str], table: pd.DataFrame, verbose: bool = True
) -> pd.Series:
    """Map typed allele names onto ligand groups.

    Matches the full name first, then falls back to progressively fewer fields
    (so C*07:01:01:14 finds C*07:01:01 or C*07:01). A group is only used as a
    fallback when every allele sharing that prefix agrees on it — otherwise the
    allele is left unassigned rather than being given a coin-flip.
    """
    exact = dict(zip(table["allele"], table["ligand_group"]))

    # prefix -> group, but only where the prefix is unambiguous
    prefix: Dict[str, str] = {}
    pref_rows = []
    for a, g in zip(table["allele"], table["ligand_group"]):
        for k in list(_keys(a))[1:]:
            pref_rows.append((k, g))
    if pref_rows:
        pf = pd.DataFrame(pref_rows, columns=["key", "group"])
        agree = pf.groupby("key")["group"].nunique() == 1
        first = pf.groupby("key")["group"].first()
        prefix = first[agree].to_dict()

    out, miss, ambig = {}, [], []
    for allele in alleles:
        hit = None
        for k in _keys(allele):
            if k in exact:
                hit = exact[k]
                break
            if k in prefix:
                hit = prefix[k]
                break
        if hit is None:
            # distinguish "never heard of it" from "prefix disagrees"
            if any(k in set(prefix) for k in _keys(allele)):
                ambig.append(allele)
            else:
                miss.append(allele)
            hit = "unassigned"
        out[allele] = hit

    if verbose:
        s = pd.Series(out)
        named = s[s != "unassigned"]
        print(f"[info] assigned {len(named)}/{len(s)} alleles")
        if miss:
            print(f"[warn] {len(miss)} not found in IPD-IMGT/HLA, e.g. {miss[:5]}")
        if ambig:
            print(f"[warn] {len(ambig)} ambiguous at this resolution, e.g. {ambig[:5]}")
    return pd.Series(out, name="ligand_group")
