#!/usr/bin/env python3
"""Fixed-effect inverse-variance meta-analysis of SAIGE step 2 runs.

Pairs phecodes across ancestries (AFR, AMR, and EUR), concatenates
chromosomes, and runs METAL SCHEME STDERR. A phecode is pooled when at least
two ancestries have results. SAIGE's BETA is the effect of Allele2. METAL
reports its Effect for Allele1 after lowercasing alleles and ordering them
lexicographically; META_BETA is that estimate flipped back onto Allele2.

The studies must be disjoint. SCHEME STDERR cannot correct for sample
overlap, and the sample-size scheme is not used. Direction follows the
ancestries actually processed for that phecode, in AFR, AMR, EUR order;
the studies column names those characters.

  saige_metal.py index --study1-list a.list --study1-label AFR --study2-list b.list --study2-label AMR ...
  saige_metal.py meta  --bundle phecode.tar.gz --out-dir .
  saige_metal.py self-test
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import math
import os
import re
import shutil
import subprocess
import sys
import tarfile
from collections import Counter, defaultdict


STEP2_RE = re.compile(
    r"(?P<phecode>.+)_(?P<chrom>(?:chr)?[0-9XYMT]+)_singlevar(?:\.txt)?(?:\.gz)?$",
    re.IGNORECASE,
)
MARKER_CANDIDATES = ("MarkerID", "SNPID", "SNP", "ID")
P_CANDIDATES = ("p.value", "P", "pval", "PVALUE", "p_value")
AF_CANDIDATES = ("AF_Allele2", "AF")
N_PAIR_CANDIDATES = (
    ("N_case", "N_ctrl"),
    ("N_cases", "N_ctrls"),
    ("N.Cases", "N.Controls"),
    ("N.case", "N.ctrl"),
)
LABEL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
MISSING = {"", "NA", "NaN", "nan", "NAN", ".", "None"}


def parse_saige_basename(name):
    """Return (phecode, chrom) from a SAIGE step2 single-variant basename.

    The split is the rightmost `_step2Out_`, so an output prefix that itself
    contains that token does not swallow the phecode.
    """
    lower = name.lower()
    key = "_step2out_"
    i = lower.rfind(key)
    if i < 0:
        return None
    m = STEP2_RE.match(name[i + len(key) :])
    if not m:
        return None
    return m.group("phecode"), m.group("chrom")


def chrom_sort_key(chrom):
    c = chrom.lower()
    if c.startswith("chr"):
        c = c[3:]
    special = {"x": 100, "y": 101, "m": 102, "mt": 102}
    if c in special:
        return (special[c], c)
    if c.isdigit():
        return (int(c), c)
    return (1000, c)


def file_sort_key(path):
    base = os.path.basename(path)
    parsed = parse_saige_basename(base)
    if parsed is None and "__" in base:
        parts = base.split("__", 2)
        if len(parts) == 3:
            parsed = parse_saige_basename(parts[2])
            if parsed is None:
                return (chrom_sort_key(parts[1]), base)
    if parsed is not None:
        return (chrom_sort_key(parsed[1]), base)
    return ((1000, ""), base)


def safe_phecode(phecode):
    safe = re.sub(r"[^A-Za-z0-9._+-]+", "_", phecode).strip("._")
    if not safe or safe in {".", ".."}:
        raise ValueError("phecode %r cannot be used as a file name" % phecode)
    return safe


def read_list(path):
    paths = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                paths.append(line)
    return paths


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "rt", newline="")


def parse_float(token):
    if token is None or token in MISSING:
        return None
    try:
        value = float(token)
    except ValueError:
        return None
    if value != value or value in (float("inf"), float("-inf")):
        return None
    return value


def format_n(value):
    if abs(value - round(value)) < 1e-6:
        return str(int(round(value)))
    return format(value, ".6g")


def negate_printed(text):
    """Flip the sign of a printed decimal without changing its digits."""
    if text.startswith("-"):
        return text[1:]
    if text.startswith("+"):
        return "-" + text[1:]
    return "-" + text


def complement_af(text):
    value = parse_float(text)
    if value is None:
        return ""
    return "%.4f" % (1.0 - value)


def resolve_columns(header, path):
    idx = {name: i for i, name in enumerate(header)}
    if len(idx) != len(header):
        raise ValueError("duplicate column names in %s" % path)

    def need(candidates, what):
        for name in candidates:
            if name in idx:
                return idx[name]
        raise ValueError(
            "%s has no %s column (looked for %s); header is %s"
            % (path, what, ", ".join(candidates), " ".join(header))
        )

    cols = {
        "marker": need(MARKER_CANDIDATES, "marker id"),
        "allele1": need(("Allele1",), "Allele1"),
        "allele2": need(("Allele2",), "Allele2"),
        "beta": need(("BETA",), "BETA"),
        "se": need(("SE",), "SE"),
        "p": need(P_CANDIDATES, "p-value"),
        "af": need(AF_CANDIDATES, "allele frequency"),
        "chr": idx.get("CHR"),
        "pos": idx.get("POS"),
        "n": idx.get("N"),
        "n_case": None,
        "n_ctrl": None,
    }
    if cols["n"] is None:
        for case_name, ctrl_name in N_PAIR_CANDIDATES:
            if case_name in idx and ctrl_name in idx:
                cols["n_case"] = idx[case_name]
                cols["n_ctrl"] = idx[ctrl_name]
                break
        if cols["n_case"] is None:
            raise ValueError(
                "%s has no sample-size column (N, or N_case and N_ctrl)" % path
            )
    return cols


def load_study(paths):
    """Read SAIGE single-variant files into one marker table.

    BETA is kept as the effect of Allele2. The first valid row wins when a
    marker is repeated across chromosomes.
    """
    rows = {}
    n_in = 0
    n_dup = 0
    dropped = Counter()
    for path in sorted(paths, key=file_sort_key):
        with open_text(path) as fh:
            header = None
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                header = line.rstrip("\r\n").split()
                break
            if not header:
                raise ValueError("no header in %s" % path)
            file_cols = resolve_columns(header, path)
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                tok = line.rstrip("\r\n").split()
                n_in += 1
                if len(tok) <= max(i for i in file_cols.values() if i is not None):
                    dropped["short_row"] += 1
                    continue
                marker = tok[file_cols["marker"]]
                allele1 = tok[file_cols["allele1"]]
                allele2 = tok[file_cols["allele2"]]
                if marker in MISSING:
                    dropped["marker"] += 1
                    continue
                if allele1 in MISSING or allele2 in MISSING or allele1 == allele2:
                    dropped["allele"] += 1
                    continue
                beta = parse_float(tok[file_cols["beta"]])
                se = parse_float(tok[file_cols["se"]])
                p = parse_float(tok[file_cols["p"]])
                af = parse_float(tok[file_cols["af"]])
                if beta is None:
                    dropped["beta"] += 1
                    continue
                if se is None or se <= 0:
                    dropped["se"] += 1
                    continue
                if p is None or p < 0 or p > 1:
                    dropped["p"] += 1
                    continue
                if af is None or af < 0 or af > 1:
                    dropped["af"] += 1
                    continue
                if file_cols["n"] is not None:
                    n = parse_float(tok[file_cols["n"]])
                else:
                    n_case = parse_float(tok[file_cols["n_case"]])
                    n_ctrl = parse_float(tok[file_cols["n_ctrl"]])
                    n = None if n_case is None or n_ctrl is None else n_case + n_ctrl
                if n is None or n <= 0:
                    dropped["n"] += 1
                    continue
                if marker in rows:
                    n_dup += 1
                    continue
                p_tok = "1e-300" if p == 0 else tok[file_cols["p"]]
                rows[marker] = {
                    "marker": marker,
                    "allele1": allele1,
                    "allele2": allele2,
                    "beta": tok[file_cols["beta"]],
                    "se": tok[file_cols["se"]],
                    "p": p_tok,
                    "af": tok[file_cols["af"]],
                    "n": format_n(n),
                    "chr": tok[file_cols["chr"]] if file_cols["chr"] is not None else "",
                    "pos": tok[file_cols["pos"]] if file_cols["pos"] is not None else "",
                }
    return {
        "rows": rows,
        "n_in": n_in,
        "n_used": len(rows),
        "n_dup": n_dup,
        "dropped": dropped,
    }


def write_metal_input(table, path):
    """METAL input. EFFECT_ALLELE is SAIGE Allele2, the allele BETA refers to."""
    fields = ("MARKER", "EFFECT_ALLELE", "OTHER_ALLELE", "BETA", "SE", "P", "AF", "N")
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE)
        writer.writerow(fields)
        for rec in table["rows"].values():
            writer.writerow(
                (
                    rec["marker"],
                    rec["allele2"],
                    rec["allele1"],
                    rec["beta"],
                    rec["se"],
                    rec["p"],
                    rec["af"],
                    rec["n"],
                )
            )


def write_metal_script(path, labels):
    """SCHEME STDERR is inverse-variance weighting. SAMPLESIZE is not used.

    PROCESS order is Direction order. One line per ancestry that has results
    for this phecode.
    """
    direction = "\n".join(
        "# Direction character %d = %s" % (i, label) for i, label in enumerate(labels, 1)
    )
    processes = "\n".join("PROCESS %s.metal.tsv" % label for label in labels)
    script = """# Fixed-effect inverse-variance meta-analysis.
# EFFECT_ALLELE is SAIGE Allele2 (the allele BETA refers to).
%s
# The studies are disjoint. OVERLAP cannot be combined with SCHEME STDERR.
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ OFF
GENOMICCONTROL OFF
SEPARATOR TAB
COLUMNCOUNTING LENIENT
MAXWARNINGS 1000
CUSTOMVARIABLE N
MARKER MARKER
ALLELE EFFECT_ALLELE OTHER_ALLELE
EFFECT BETA
STDERR SE
PVALUE P
FREQ AF
OUTFILE meta .tbl
%s
ANALYZE HETEROGENEITY
QUIT
""" % (
        direction,
        processes,
    )
    with open(path, "w") as fh:
        fh.write(script)


def align_to_allele(rec, target_allele):
    """Return (beta, se, af) for target_allele, or None if it is not in this row."""
    if rec is None:
        return None
    beta = parse_float(rec["beta"])
    se = parse_float(rec["se"])
    af = parse_float(rec["af"])
    if beta is None or se is None or se <= 0 or af is None:
        return None
    if rec["allele2"].lower() == target_allele.lower():
        return beta, se, af
    if rec["allele1"].lower() == target_allele.lower():
        return -beta, se, 1.0 - af
    return None


def inverse_variance(parts):
    """Fixed-effect inverse-variance pool of (beta, se, af) on one allele."""
    weight_sum = 0.0
    beta_sum = 0.0
    af_sum = 0.0
    for beta, se, af in parts:
        weight = 1.0 / (se * se)
        weight_sum += weight
        beta_sum += weight * beta
        af_sum += weight * af
    return beta_sum / weight_sum, 1.0 / math.sqrt(weight_sum), af_sum / weight_sum


def orient_effect(meta_allele1, effect, freq1, rec):
    """Return (META_BETA, META_AF) on SAIGE Allele2, or (None, None)."""
    if rec is None or effect == "":
        return None, None
    meta = meta_allele1.lower()
    if meta == rec["allele2"].lower():
        return effect, freq1
    if meta == rec["allele1"].lower():
        af = complement_af(freq1) if freq1 != "" else ""
        return negate_printed(effect), af
    return None, None


def n_studies(direction):
    return sum(ch != "?" for ch in direction)


def result_columns(labels):
    per_study = []
    for label in labels:
        for suffix in ("BETA", "SE", "P", "AF", "N"):
            per_study.append("%s_%s" % (label, suffix))
    return [
        "phecode",
        "CHR",
        "POS",
        "MarkerName",
        "SAIGE_Allele1",
        "SAIGE_Allele2",
        "allele_source",
        "META_BETA",
        "META_SE",
        "P-value",
        "Direction",
        "studies",
        "n_studies",
        "META_AF",
        "N",
        "HetISq",
        "HetChiSq",
        "HetDf",
        "HetPVal",
        "Allele1",
        "Allele2",
        "Effect",
        "StdErr",
        "Freq1",
        "FreqSE",
    ] + per_study


def _study_fields(label, rec):
    if rec is None:
        return {("%s_%s" % (label, suffix)): "" for suffix in ("BETA", "SE", "P", "AF", "N")}
    return {
        "%s_BETA" % label: rec["beta"],
        "%s_SE" % label: rec["se"],
        "%s_P" % label: rec["p"],
        "%s_AF" % label: rec["af"],
        "%s_N" % label: rec["n"],
    }


def build_result(tbl_path, phecode, labels, tables, all_labels, both_studies_only, out_path):
    """Write the result table.

    `labels` is the PROCESS order (Direction order) for this phecode.
    `all_labels` is the workflow order, used so every phecode file has the
    same ancestry columns. `tables` maps a label to a load_study() result.
    """
    columns = result_columns(all_labels)
    studies = ",".join(labels)
    n_meta = 0
    n_complete = 0
    n_partial = 0
    n_unoriented = 0
    n_discordant = 0
    n_written = 0
    min_p = None
    with open(tbl_path, newline="") as src, open(out_path, "w", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        missing = [c for c in ("MarkerName", "Allele1", "Effect", "StdErr", "P-value", "Direction") if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError("METAL table %s is missing %s" % (tbl_path, ", ".join(missing)))
        writer = csv.DictWriter(
            dst,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in reader:
            n_meta += 1
            direction = row.get("Direction") or ""
            n_here = n_studies(direction)
            complete = "?" not in direction and n_here >= 2 and n_here == len(labels)
            if complete:
                n_complete += 1
            else:
                n_partial += 1
            if both_studies_only and not complete:
                continue
            marker = row["MarkerName"]
            present = []
            for label in labels:
                rec_l = tables[label]["rows"].get(marker)
                if rec_l is not None:
                    present.append((label, rec_l))
            if present:
                source, rec = present[0]
            else:
                source, rec = "", None
            beta, _af = orient_effect(row.get("Allele1") or "", row.get("Effect") or "", row.get("Freq1") or "", rec)
            parts = []
            if rec is not None:
                for _label, candidate in present:
                    aligned = align_to_allele(candidate, rec["allele2"])
                    if aligned is not None:
                        parts.append(aligned)
            if parts:
                pooled_beta, pooled_se, pooled_af = inverse_variance(parts)
                meta_beta = format(pooled_beta, ".12g")
                meta_se = format(pooled_se, ".12g")
                meta_af = format(pooled_af, ".6g")
                if beta is not None and abs(float(beta) - pooled_beta) > 5e-4:
                    n_discordant += 1
            else:
                n_unoriented += 1
                meta_beta = ""
                meta_se = ""
                meta_af = ""
            p = parse_float(row.get("P-value"))
            if p is not None and (min_p is None or p < min_p):
                min_p = p
            out = {key: "" for key in columns}
            out.update(
                {
                    "phecode": phecode,
                    "CHR": "" if rec is None else rec["chr"],
                    "POS": "" if rec is None else rec["pos"],
                    "MarkerName": marker,
                    "SAIGE_Allele1": "" if rec is None else rec["allele1"],
                    "SAIGE_Allele2": "" if rec is None else rec["allele2"],
                    "allele_source": source,
                    "META_BETA": meta_beta,
                    "META_SE": meta_se,
                    "P-value": row.get("P-value") or "",
                    "Direction": direction,
                    "studies": studies,
                    "n_studies": str(n_here),
                    "META_AF": meta_af,
                    "N": row.get("N") or "",
                    "HetISq": row.get("HetISq") or "",
                    "HetChiSq": row.get("HetChiSq") or "",
                    "HetDf": row.get("HetDf") or "",
                    "HetPVal": row.get("HetPVal") or "",
                    "Allele1": row.get("Allele1") or "",
                    "Allele2": row.get("Allele2") or "",
                    "Effect": row.get("Effect") or "",
                    "StdErr": row.get("StdErr") or "",
                    "Freq1": row.get("Freq1") or "",
                    "FreqSE": row.get("FreqSE") or "",
                }
            )
            for label in all_labels:
                table = tables.get(label)
                rec_l = None if table is None else table["rows"].get(marker)
                out.update(_study_fields(label, rec_l))
            writer.writerow(out)
            n_written += 1
    return {
        "n_meta": n_meta,
        "n_complete": n_complete,
        "n_partial": n_partial,
        "n_unoriented": n_unoriented,
        "n_discordant": n_discordant,
        "n_written": n_written,
        "min_p": "" if min_p is None else format(min_p, ".6g"),
    }


def append_info(path, phecode, labels, ran):
    with open(path, "a") as fh:
        fh.write("\n")
        fh.write("# phecode: %s\n" % phecode)
        for i, label in enumerate(labels, 1):
            fh.write("# Direction character %d = %s\n" % (i, label))
        fh.write("# Ancestries with no results for this phecode are omitted from Direction.\n")
        fh.write("# The studies column lists Direction in order.\n")
        fh.write("# META_BETA and META_SE are the fixed-effect inverse-variance estimate on SAIGE Allele2.\n")
        fh.write("# They are pooled from the SAIGE betas as sum(beta/SE^2) / sum(1/SE^2), not from METAL's rounded Effect.\n")
        fh.write("# Effect and StdErr are METAL's printed values; Effect is for Allele1 after lexicographic allele ordering.\n")
        fh.write("# N is the sum of the study sample sizes. The studies are disjoint.\n")
        fh.write("# Freq1 is the inverse-variance-weighted frequency of Allele1.\n")
        fh.write("# GENOMICCONTROL was OFF. Weighting is 1/SE^2 (SCHEME STDERR).\n")
        if not ran:
            fh.write("# METAL was not run: fewer than two ancestries had usable variants.\n")


def write_qc(path, row):
    fields = list(row.keys())
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def reason_string(counter):
    if not counter:
        return ""
    return ";".join("%s=%d" % (key, counter[key]) for key in sorted(counter))


def qc_row(phecode, all_labels, included, status, tables, file_counts, counts):
    row = {
        "phecode": phecode,
        "studies": ",".join(included),
        "status": status,
    }
    for label in all_labels:
        table = tables.get(label)
        row["n_%s_files" % label] = file_counts.get(label, 0)
        row["n_%s_rows_in" % label] = 0 if table is None else table["n_in"]
        row["n_%s_used" % label] = 0 if table is None else table["n_used"]
        row["n_%s_dup" % label] = 0 if table is None else table["n_dup"]
        row["drop_reasons_%s" % label] = "" if table is None else reason_string(table["dropped"])
    row["n_meta"] = counts.get("n_meta", 0)
    row["n_complete"] = counts.get("n_complete", 0)
    row["n_partial"] = counts.get("n_partial", 0)
    row["n_written"] = counts.get("n_written", 0)
    row["n_unoriented"] = counts.get("n_unoriented", 0)
    row["n_discordant"] = counts.get("n_discordant", 0)
    row["min_p"] = counts.get("min_p", "")
    return row


def write_empty_result(path, labels):
    columns = result_columns(labels)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()


def check_labels(labels):
    if len(labels) < 2:
        raise SystemExit("at least two studies are required")
    for label in labels:
        if not LABEL_RE.match(label):
            raise SystemExit(
                "study label %r must match %s (it is used as a column name and a directory name)"
                % (label, LABEL_RE.pattern)
            )
    if len(set(labels)) != len(labels):
        raise SystemExit("study labels must be unique, got %s" % ", ".join(labels))


def bundle_member_ok(name):
    if name in {"phecode.txt", "studies.txt", "all_labels.txt"}:
        return True
    parts = name.split("/")
    if len(parts) != 2 or not parts[1] or parts[1] in {".", ".."}:
        return False
    return LABEL_RE.match(parts[0]) is not None


def extract_bundle(bundle, dest):
    os.makedirs(dest, exist_ok=True)
    with tarfile.open(bundle, "r:*") as tar:
        for member in tar.getmembers():
            name = member.name
            if name.startswith("/") or ".." in name.split("/"):
                raise ValueError("unsafe path in bundle: %s" % name)
            if not bundle_member_ok(name):
                raise ValueError("unexpected bundle member: %s" % name)
            tar.extract(member, dest)


def add_file(tar, arcname, path):
    info = tarfile.TarInfo(name=arcname)
    info.size = os.path.getsize(path)
    with open(path, "rb") as fh:
        tar.addfile(info, fh)


def add_text(tar, arcname, text):
    data = text.encode()
    info = tarfile.TarInfo(name=arcname)
    info.size = len(data)
    tar.addfile(info, io.BytesIO(data))


def read_label_file(path):
    labels = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                labels.append(line)
    check_labels(labels) if len(labels) >= 2 else None
    return labels


def parse_study_args(items):
    """--study LABEL=file_list, repeated in Direction order."""
    studies = []
    for item in items:
        if "=" not in item:
            raise SystemExit("--study must be LABEL=path, got %r" % item)
        label, path = item.split("=", 1)
        studies.append((label, path))
    check_labels([label for label, _path in studies])
    return studies


def collect_studies(args):
    """Studies in Direction order.

    Repeated --study LABEL=list, or --studyN-list with --studyN-label for
    N = 1..5. The numbered flags are what the WDL passes. --study is not a
    prefix of those names when abbreviation is off; with abbreviation on,
    argparse rejects --study as ambiguous.
    """
    repeated = getattr(args, "study", None) or []
    numbered = []
    gap = False
    for i in range(1, 6):
        path = getattr(args, "study%d_list" % i, None)
        label = getattr(args, "study%d_label" % i, None)
        if path is None and label is None:
            gap = True
            continue
        if gap:
            raise SystemExit("--study%d-list is set but a preceding study was omitted" % i)
        if path is None or label is None:
            raise SystemExit("--study%d-list and --study%d-label must be set together" % (i, i))
        numbered.append((label, path))
    if repeated and numbered:
        raise SystemExit("pass either repeated --study LABEL=list or --studyN-list/--studyN-label, not both")
    if repeated:
        return parse_study_args(repeated)
    check_labels([label for label, _path in numbered])
    return numbered


def cmd_index(args):
    studies = collect_studies(args)
    labels = [label for label, _path in studies]
    owner = {}
    lists = {}
    for label, list_path in studies:
        paths = []
        seen = set()
        for path in read_list(list_path):
            if path in seen:
                continue
            seen.add(path)
            if path in owner:
                raise SystemExit("%s is listed for both %s and %s" % (path, owner[path], label))
            owner[path] = label
            paths.append(path)
        lists[label] = paths

    grouped = {label: defaultdict(lambda: defaultdict(list)) for label in labels}
    skipped = []
    for label, paths in lists.items():
        for path in paths:
            if not os.path.isfile(path):
                skipped.append((path, "missing"))
                continue
            if os.path.getsize(path) == 0:
                skipped.append((path, "empty"))
                continue
            parsed = parse_saige_basename(os.path.basename(path))
            if parsed is None:
                skipped.append((path, "unrecognized_name"))
                continue
            phecode, chrom = parsed
            grouped[label][phecode][chrom].append(path)

    phecodes = sorted(set().union(*[set(grouped[label]) for label in labels]))
    os.makedirs(args.bundle_dir, exist_ok=True)
    pairing = []
    safe_used = {}
    n_paired = 0
    for phecode in phecodes:
        by_label = {label: grouped[label].get(phecode, {}) for label in labels}
        present = [label for label in labels if by_label[label]]
        dup = [chrom for label in present for chrom, paths in by_label[label].items() if len(paths) > 1]
        if dup:
            status = "duplicate_chrom"
        elif len(present) >= 2:
            status = "+".join(present)
        else:
            status = "only_%s" % present[0]
        row = {"phecode": phecode, "status": status, "studies": ",".join(present)}
        for label in labels:
            chroms = by_label[label]
            row["n_%s_files" % label] = sum(len(v) for v in chroms.values())
            row["%s_chroms" % label] = ",".join(sorted(chroms, key=chrom_sort_key))
        pairing.append(row)
        if len(present) < 2 or dup:
            continue
        safe = safe_phecode(phecode)
        if safe in safe_used:
            raise SystemExit("phecode file names collide: %s and %s" % (safe_used[safe], phecode))
        safe_used[safe] = phecode
        bundle_path = os.path.join(args.bundle_dir, safe + ".tar.gz")
        with tarfile.open(bundle_path, "w:gz") as tar:
            add_text(tar, "phecode.txt", phecode + "\n")
            add_text(tar, "studies.txt", "".join(label + "\n" for label in present))
            add_text(tar, "all_labels.txt", "".join(label + "\n" for label in labels))
            for label in present:
                items = []
                for chrom in sorted(by_label[label], key=chrom_sort_key):
                    items.append((chrom, by_label[label][chrom][0]))
                for i, (chrom, path) in enumerate(items):
                    arc = "%s/%04d__%s__%s" % (label, i, chrom, os.path.basename(path))
                    add_file(tar, arc, path)
        n_paired += 1

    fields = ["phecode", "status", "studies"]
    for label in labels:
        fields.append("n_%s_files" % label)
        fields.append("%s_chroms" % label)
    with open(args.pairing, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(pairing)
    with open(args.direction_order, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(("position", "study"))
        for i, label in enumerate(labels, 1):
            writer.writerow((i, label))
    with open(args.skipped, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(("path", "reason"))
        writer.writerows(skipped)

    print(
        "paired %d phecode(s); %d unpaired or duplicate; %d skipped file(s)"
        % (n_paired, len(phecodes) - n_paired, len(skipped))
    )
    if n_paired == 0:
        raise SystemExit("no phecode is present in at least two studies")



def cmd_meta(args):
    out = args.out_dir
    os.makedirs(out, exist_ok=True)
    dest = os.path.join(out, "inbundle")
    if os.path.exists(dest):
        shutil.rmtree(dest)
    extract_bundle(args.bundle, dest)
    with open(os.path.join(dest, "phecode.txt")) as fh:
        phecode = fh.read().strip()
    if not phecode:
        raise SystemExit("bundle has an empty phecode.txt")
    included = read_label_file(os.path.join(dest, "studies.txt"))
    all_path = os.path.join(dest, "all_labels.txt")
    all_labels = read_label_file(all_path) if os.path.isfile(all_path) else list(included)
    if [label for label in all_labels if label in included] != included:
        raise SystemExit("studies.txt is not in all_labels.txt order")

    def members(study):
        directory = os.path.join(dest, study)
        if not os.path.isdir(directory):
            return []
        paths = [
            os.path.join(directory, name)
            for name in os.listdir(directory)
            if os.path.isfile(os.path.join(directory, name))
        ]
        return sorted(paths, key=file_sort_key)

    tables = {}
    file_counts = {}
    for label in included:
        files = members(label)
        file_counts[label] = len(files)
        tables[label] = load_study(files)
    usable = [label for label in included if tables[label]["n_used"] > 0]
    result_path = os.path.join(out, "metal_result.tsv")
    info_path = os.path.join(out, "metal_result.info")
    log_path = os.path.join(out, "metal.log")
    script_path = os.path.join(out, "metal_script.txt")
    qc_path = os.path.join(out, "metal_qc.tsv")

    if len(usable) < 2:
        write_empty_result(result_path, all_labels)
        with open(info_path, "w") as fh:
            fh.write("# METAL was not run\n")
        append_info(info_path, phecode, usable, ran=False)
        with open(log_path, "w") as fh:
            fh.write("METAL was not run: fewer than two ancestries had usable variants\n")
        with open(script_path, "w") as fh:
            fh.write("# METAL was not run: fewer than two ancestries had usable variants\n")
        write_qc(
            qc_path,
            qc_row(phecode, all_labels, usable, "no_variants", tables, file_counts, {}),
        )
        print("%s\tno_variants" % phecode)
        return

    # METAL reads these names from the working directory.
    cwd = os.getcwd()
    try:
        os.chdir(out)
        for label in usable:
            write_metal_input(tables[label], "%s.metal.tsv" % label)
        write_metal_script("metal_script.txt", usable)
        with open("metal.log", "w") as log:
            proc = subprocess.run(
                ["metal", "metal_script.txt"],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        log_text = open("metal.log").read()
    finally:
        os.chdir(cwd)
    if proc.returncode != 0 or "Completed meta-analysis" not in log_text:
        raise SystemExit("METAL failed for phecode %s; see metal.log" % phecode)
    tbl_path = os.path.join(out, "meta1.tbl")
    if not os.path.isfile(tbl_path):
        raise SystemExit("METAL did not write meta1.tbl for phecode %s" % phecode)
    used_tables = {label: tables[label] for label in usable}
    counts = build_result(
        tbl_path,
        phecode,
        usable,
        used_tables,
        all_labels,
        args.both_studies_only,
        result_path,
    )
    info_src = os.path.join(out, "meta1.tbl.info")
    if os.path.isfile(info_src):
        shutil.copyfile(info_src, info_path)
    else:
        with open(info_path, "w") as fh:
            fh.write("# METAL info file was missing\n")
    append_info(info_path, phecode, usable, ran=True)
    write_qc(
        qc_path,
        qc_row(phecode, all_labels, usable, "ok", tables, file_counts, counts),
    )
    print(
        "%s\tok\twritten=%d\tcomplete=%d\tunoriented=%d\tdiscordant=%d"
        % (phecode, counts["n_written"], counts["n_complete"], counts["n_unoriented"], counts["n_discordant"])
    )

def self_test():
    # Allele orientation. METAL's Effect is for its printed Allele1.
    assert negate_printed("-0.200000000000") == "0.200000000000"
    assert negate_printed("0.200000000000") == "-0.200000000000"
    rec = {"allele1": "A", "allele2": "T"}
    # Printed allele is A, so Effect is -BETA and must be flipped back onto T.
    beta, af = orient_effect("a", "-0.200000000000", "0.9500", rec)
    assert beta == "0.200000000000", beta
    assert af == "0.0500", af
    beta, af = orient_effect("t", "0.200000000000", "0.0500", rec)
    assert beta == "0.200000000000" and af == "0.0500"
    # Two studies, both already on Allele2=T: 0.2 (SE 0.1) and 0.4 (SE 0.2).
    pooled_beta, pooled_se, pooled_af = inverse_variance(
        [(0.2, 0.1, 0.05), (0.4, 0.2, 0.08)]
    )
    assert abs(pooled_beta - 0.24) < 1e-12, pooled_beta
    assert abs(pooled_se - 1.0 / math.sqrt(125.0)) < 1e-12, pooled_se
    assert abs(pooled_af - (100 * 0.05 + 25 * 0.08) / 125.0) < 1e-12
    # AFR 0.2/0.1, AMR 0.4/0.2, EUR 0.1/0.1.
    three_beta, three_se, _three_af = inverse_variance(
        [(0.2, 0.1, 0.05), (0.4, 0.2, 0.08), (0.1, 0.1, 0.02)]
    )
    assert abs(three_beta - 40.0 / 225.0) < 1e-12, three_beta
    assert abs(three_se - 1.0 / 15.0) < 1e-12, three_se
    hla = {"allele1": "N", "allele2": "HLA-A*02:01"}
    beta, _af = orient_effect("hla-a*02:01", "0.400000000000", "0.1000", hla)
    assert beta == "0.400000000000"
    assert orient_effect("g", "0.1", "0.2", rec) == (None, None)

    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        s1_chr6 = os.path.join(tmp, "afr_step2Out_GI_519.2_chr6_singlevar")
        s1_chr10 = os.path.join(tmp, "afr_step2Out_GI_519.2_chr10_singlevar")
        s1_only = os.path.join(tmp, "afr_step2Out_ONLY.1_chr6_singlevar")
        s2 = os.path.join(tmp, "amr_step2Out_GI_519.2_6_singlevar.gz")
        s3 = os.path.join(tmp, "eur_step2Out_GI_519.2_chr6_singlevar")
        s1_ae = os.path.join(tmp, "afr_step2Out_AE.1_chr6_singlevar")
        s3_ae = os.path.join(tmp, "eur_step2Out_AE.1_chr6_singlevar")
        s3_only = os.path.join(tmp, "eur_step2Out_EUR.9_chr6_singlevar")
        unrecognized = os.path.join(tmp, "failed_step2_singlevar.txt")
        header = "CHR POS MarkerID Allele1 Allele2 AC_Allele2 AF_Allele2 MissingRate BETA SE Tstat var p.value N\n"
        with open(s1_chr6, "w") as fh:
            fh.write(header)
            fh.write("6 100 m1 A T 10 0.05 0 0.2 0.1 1 1 0.045 1000\n")
            fh.write("6 200 m2 G C 10 0.10 0 0.5 0.2 1 1 0.012 1000\n")
            fh.write("6 300 mdup A G 1 0.01 0 0.1 0.1 1 1 0.3 1000\n")
            fh.write("6 301 mdup A G 1 0.01 0 9 9 1 1 0.9 1000\n")
            fh.write("6 400 mzero A C 1 0.01 0 0.2 0.1 1 1 0 1000\n")
            fh.write("6 500 mbad A C 1 0.01 0 NA 0.1 1 1 0.5 1000\n")
        with open(s1_chr10, "w") as fh:
            fh.write(header)
            fh.write("10 50 m3 A G 4 0.02 0 -0.4 0.2 1 1 0.04 1000\n")
        with open(s1_only, "w") as fh:
            fh.write(header)
            fh.write("6 1 monly A T 1 0.01 0 0.1 0.1 1 1 0.3 800\n")
        with gzip.open(s2, "wt") as fh:
            fh.write(
                "CHR POS SNPID Allele1 Allele2 AC_Allele2 AF_Allele2 imputationInfo "
                "BETA SE Tstat var p.value N_case N_ctrl\n"
            )
            fh.write("6 100 m1 A T 8 0.08 1 0.4 0.2 1 1 0.045 200 800\n")
            fh.write("6 200 m2 G C 8 0.08 1 -0.1 0.3 1 1 0.7 200 800\n")
            fh.write("6 900 m4 A T 8 0.08 1 1.0 0.5 1 1 0.045 200 800\n")
        with open(s3, "w") as fh:
            fh.write(header)
            fh.write("6 100 m1 A T 10 0.02 0 0.1 0.1 1 1 0.3 900\n")
        with open(s1_ae, "w") as fh:
            fh.write(header)
            fh.write("6 1 aonly A T 1 0.01 0 0.2 0.1 1 1 0.05 800\n")
        with open(s3_ae, "w") as fh:
            fh.write(header)
            fh.write("6 1 aonly A T 1 0.01 0 0.3 0.2 1 1 0.1 700\n")
        with open(s3_only, "w") as fh:
            fh.write(header)
            fh.write("6 1 eonly A T 1 0.01 0 0.1 0.1 1 1 0.4 600\n")
        with open(unrecognized, "w") as fh:
            fh.write("not a result\n")
        # Prefix itself contains _step2Out_; the phecode is still GI_519.2.
        weird = os.path.join(tmp, "pre_step2Out_x_step2Out_GI_519.2_chr6_singlevar")
        assert parse_saige_basename(os.path.basename(weird)) == ("GI_519.2", "chr6")

        loaded = load_study([s1_chr10, s1_chr6])
        assert loaded["n_used"] == 5, loaded["n_used"]  # m1 m2 mdup mzero m3
        assert loaded["n_dup"] == 1
        assert loaded["dropped"]["beta"] == 1
        assert list(loaded["rows"]) == ["m1", "m2", "mdup", "mzero", "m3"]
        assert loaded["rows"]["m1"]["allele2"] == "T"
        assert loaded["rows"]["m1"]["beta"] == "0.2"
        assert loaded["rows"]["mzero"]["p"] == "1e-300"
        loaded2 = load_study([s2])
        assert loaded2["rows"]["m1"]["n"] == "1000"
        assert loaded2["rows"]["m1"]["allele2"] == "T"
        loaded3 = load_study([s3])
        assert loaded3["rows"]["m1"]["beta"] == "0.1"

        bundle_dir = os.path.join(tmp, "bundles")
        list1 = os.path.join(tmp, "s1.list")
        list2 = os.path.join(tmp, "s2.list")
        list3 = os.path.join(tmp, "s3.list")
        with open(list1, "w") as fh:
            fh.write("\n".join([s1_chr6, s1_chr10, s1_only, s1_ae, unrecognized, s1_chr6]) + "\n")
        with open(list2, "w") as fh:
            fh.write(s2 + "\n")
        with open(list3, "w") as fh:
            fh.write("\n".join([s3, s3_ae, s3_only]) + "\n")
        ns = argparse.Namespace(
            study=["AFR=" + list1, "AMR=" + list2, "EUR=" + list3],
            bundle_dir=bundle_dir,
            pairing=os.path.join(tmp, "pairing.tsv"),
            direction_order=os.path.join(tmp, "direction.tsv"),
            skipped=os.path.join(tmp, "skipped.tsv"),
        )
        cmd_index(ns)
        assert set(os.listdir(bundle_dir)) == {"GI_519.2.tar.gz", "AE.1.tar.gz"}
        with open(ns.pairing) as fh:
            body = list(csv.DictReader(fh, delimiter="\t"))
        status = {row["phecode"]: row["status"] for row in body}
        assert status == {
            "GI_519.2": "AFR+AMR+EUR",
            "ONLY.1": "only_AFR",
            "AE.1": "AFR+EUR",
            "EUR.9": "only_EUR",
        }, status
        with open(ns.direction_order) as fh:
            direction_rows = list(csv.DictReader(fh, delimiter="\t"))
        assert [(row["position"], row["study"]) for row in direction_rows] == [
            ("1", "AFR"),
            ("2", "AMR"),
            ("3", "EUR"),
        ]
        with open(ns.skipped) as fh:
            skipped = list(csv.DictReader(fh, delimiter="\t"))
        assert skipped == [{"path": unrecognized, "reason": "unrecognized_name"}]
        with tarfile.open(os.path.join(bundle_dir, "AE.1.tar.gz")) as tar:
            names = tar.getnames()
        assert "studies.txt" in names and any(name.startswith("AFR/") for name in names)
        assert any(name.startswith("EUR/") for name in names)
        assert not any(name.startswith("AMR/") for name in names)

        script_path = os.path.join(tmp, "metal_script.txt")
        write_metal_script(script_path, ["AFR", "AMR", "EUR"])
        script = open(script_path).read()
        for line in (
            "SCHEME STDERR",
            "GENOMICCONTROL OFF",
            "MINMAXFREQ OFF",
            "ANALYZE HETEROGENEITY",
            "CUSTOMVARIABLE N",
            "ALLELE EFFECT_ALLELE OTHER_ALLELE",
            "PROCESS AFR.metal.tsv",
            "PROCESS AMR.metal.tsv",
            "PROCESS EUR.metal.tsv",
        ):
            assert line in script, line
        assert script.index("PROCESS AFR.metal.tsv") < script.index("PROCESS AMR.metal.tsv")
        assert script.index("PROCESS AMR.metal.tsv") < script.index("PROCESS EUR.metal.tsv")
        assert "SCHEME SAMPLESIZE" not in script
        assert "OVERLAP ON" not in script

        metal_in = os.path.join(tmp, "study1.metal.tsv")
        write_metal_input(loaded, metal_in)
        with open(metal_in) as fh:
            metal_rows = list(csv.DictReader(fh, delimiter="\t"))
        by_marker = {row["MARKER"]: row for row in metal_rows}
        assert by_marker["m1"]["EFFECT_ALLELE"] == "T"
        assert by_marker["m1"]["OTHER_ALLELE"] == "A"
        assert by_marker["m1"]["BETA"] == "0.2"
        assert by_marker["m3"]["MARKER"] == "m3"

        tbl = os.path.join(tmp, "meta1.tbl")
        with open(tbl, "w") as fh:
            fh.write(
                "MarkerName\tAllele1\tAllele2\tFreq1\tFreqSE\tEffect\tStdErr\tP-value\tDirection\t"
                "HetISq\tHetChiSq\tHetDf\tHetPVal\tN\n"
            )
            # Effect is on Allele1=a, so the Allele2 (T) effect is the negation.
            fh.write("m1\ta\tt\t0.9500\t0.0100\t-0.2400\t0.0894\t0.007289\t+-\t0.0\t0.800\t1\t0.3711\t2000\n")
            fh.write("m4\ta\tt\t0.9200\t0.0000\t-1.0000\t0.5000\t0.0455\t?+\t0.0\t0.000\t0\t1\t1000\n")
        out_tsv = os.path.join(tmp, "metal_result.tsv")
        two = {"AFR": loaded, "AMR": loaded2}
        counts = build_result(
            tbl, "GI_519.2", ["AFR", "AMR"], two, ["AFR", "AMR", "EUR"], False, out_tsv
        )
        assert counts["n_written"] == 2, counts
        assert counts["n_complete"] == 1
        assert counts["n_discordant"] == 0, counts
        with open(out_tsv) as fh:
            got = list(csv.DictReader(fh, delimiter="\t"))
        assert got[0]["META_BETA"] == "0.24", got[0]["META_BETA"]
        assert abs(float(got[0]["META_SE"]) - 1.0 / math.sqrt(125.0)) < 1e-9
        assert got[0]["SAIGE_Allele2"] == "T"
        assert got[0]["AFR_BETA"] == "0.2" and got[0]["AMR_BETA"] == "0.4"
        assert got[0]["EUR_BETA"] == ""
        assert got[0]["studies"] == "AFR,AMR"
        assert got[0]["n_studies"] == "2"
        assert got[1]["MarkerName"] == "m4"
        assert got[1]["META_BETA"] == "1", got[1]["META_BETA"]
        assert got[1]["AFR_BETA"] == "" and got[1]["allele_source"] == "AMR"
        kept = build_result(
            tbl, "GI_519.2", ["AFR", "AMR"], two, ["AFR", "AMR", "EUR"], True, out_tsv + ".both"
        )
        assert kept["n_written"] == 1

        effect3 = format(three_beta, ".12g")
        tbl3 = os.path.join(tmp, "meta3.tbl")
        with open(tbl3, "w") as fh:
            fh.write(
                "MarkerName\tAllele1\tAllele2\tFreq1\tFreqSE\tEffect\tStdErr\tP-value\tDirection\t"
                "HetISq\tHetChiSq\tHetDf\tHetPVal\tN\n"
            )
            fh.write("m1\tt\ta\t0.0400\t0.0100\t%s\t0.0667\t0.008\t+++\t0.0\t1.000\t2\t0.6\t2900\n" % effect3)
            fh.write("m4\ta\tt\t0.9200\t0.0000\t-1.0000\t0.5000\t0.0455\t?+?\t0.0\t0.000\t0\t1\t1000\n")
        three = {"AFR": loaded, "AMR": loaded2, "EUR": loaded3}
        counts3 = build_result(
            tbl3, "GI_519.2", ["AFR", "AMR", "EUR"], three, ["AFR", "AMR", "EUR"], False, out_tsv + ".3"
        )
        assert counts3["n_complete"] == 1 and counts3["n_written"] == 2
        with open(out_tsv + ".3") as fh:
            got3 = list(csv.DictReader(fh, delimiter="\t"))
        assert abs(float(got3[0]["META_BETA"]) - three_beta) < 1e-9
        assert abs(float(got3[0]["META_SE"]) - three_se) < 1e-9
        assert got3[0]["EUR_BETA"] == "0.1"
        assert got3[0]["studies"] == "AFR,AMR,EUR"
        assert got3[0]["n_studies"] == "3"
        assert counts3["n_discordant"] == 0
        kept3 = build_result(
            tbl3, "GI_519.2", ["AFR", "AMR", "EUR"], three, ["AFR", "AMR", "EUR"], True, out_tsv + ".3b"
        )
        assert kept3["n_written"] == 1
    print("self-test ok")


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    sub = parser.add_subparsers(dest="command", required=True)

    index = sub.add_parser(
        "index",
        help="pair phecodes and tar each ancestry set",
        allow_abbrev=False,
    )
    index.add_argument(
        "--study",
        action="append",
        help="LABEL=file_list, repeated in Direction order. Alternative to --studyN-list.",
    )
    for i in range(1, 6):
        index.add_argument("--study%d-list" % i, help="step2 file list for study %d" % i)
        index.add_argument("--study%d-label" % i, help="Direction label for study %d" % i)
    index.add_argument("--bundle-dir", required=True)
    index.add_argument("--pairing", required=True)
    index.add_argument("--direction-order", required=True)
    index.add_argument("--skipped", required=True)
    index.set_defaults(func=cmd_index)

    meta = sub.add_parser("meta", help="run METAL on one phecode bundle")
    meta.add_argument("--bundle", required=True)
    meta.add_argument("--both-studies-only", action="store_true")
    meta.add_argument("--out-dir", required=True)
    meta.set_defaults(func=cmd_meta)

    test = sub.add_parser("self-test", help="check parsing and allele orientation")
    test.set_defaults(func=lambda _args: self_test())

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
