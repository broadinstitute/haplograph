#!/usr/bin/env python3
"""Manhattan and HLA heatmaps from SAIGE output, drawn by HATK.

HATK (https://github.com/WansonChoi/HATK) colors a Manhattan plot by marker
class and draws one amino-acid heatmap per HLA gene. Its plotting modules
read a PLINK association table, so this script rewrites SAIGE into that
table and then runs the modules.

SAIGE step 2 ``*_singlevar`` files and ``saige_metal.py`` meta tables are
both accepted. BETA / META_BETA is the log odds ratio of Allele2 /
SAIGE_Allele2; the PLINK OR column is exp of that coefficient. A p-value
of 0 is written as 1e-300 so HATK's log10 does not overflow.

Marker names HATK will color and match to a maptable:

  HLA-DRB1*04:01          ->  HLA_DRB1*04:01
  AA_DRB1_71_32659900_K   ->  kept; 71 is the IMGT amino-acid position
  AA_DRB1_71_K            ->  AA_DRB1_71_<POS>_K   (POS fills the genomic base)
  SNPS_A_30018457, rsIDs  ->  kept

The heatmap needs both ``HLA_{gene}*`` alleles and ``AA_{gene}_`` markers,
plus the maptable IMGT2Seq writes (``HLA_MAPTABLE_{gene}.hg*.imgt*.txt``).
HATK ships hg19 / IMGT 3.32.0 maptables under ``HLA_Heatmap/data/``. The
allele string after ``HLA_`` must equal a maptable row (``DRB1*04:01``,
same field resolution). Amino-acid positions must equal the maptable's
``AA_rel_pos`` header.

HATK draws the Manhattan on chr6 from 29.60 Mb to 33.20 Mb, which covers
the classical HLA genes on hg19 and hg38. Pass the genome build the
positions are on.

    saige_hatk_plot.py --self-test

    saige_hatk_plot.py \\
        --saige EUR_step2Out_phecode_chr6_singlevar.txt.gz \\
        --hatk ~/HATK --python ~/miniconda3/envs/HATK/bin/python \\
        --hg 38 --out plots/phecode

    saige_hatk_plot.py --saige meta.tsv --hatk ~/HATK --hg 38 \\
        --hla A B C DRB1 DQA1 DQB1 DPA1 DPB1 --out plots/meta
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import re
import subprocess
import sys
from collections import Counter


MISSING = {"", "NA", "NaN", "nan", "NAN", ".", "None"}
P_FLOOR = "1e-300"

# Gene symbols HATK uses without the HLA- prefix. A leading HLA- / HLA_ on
# the marker is enough on its own; this set also catches bare A*02:01.
HLA_GENES = frozenset({
    "A", "B", "C", "E", "F", "G",
    "DMA", "DMB", "DOA", "DOB",
    "DPA1", "DPB1", "DQA1", "DQA2", "DQB1", "DQB2",
    "DRA", "DRB1", "DRB3", "DRB4", "DRB5",
    "HFE", "MICA", "MICB", "TAP1", "TAP2", "V",
})
GENE_ORDER = (
    "A", "B", "C", "E", "F", "G",
    "DRA", "DRB1", "DRB3", "DRB4", "DRB5",
    "DQA1", "DQB1", "DPA1", "DPB1",
    "DMA", "DMB", "DOA", "DOB",
    "MICA", "MICB", "TAP1", "TAP2",
)

MARKER_CANDIDATES = ("MarkerName", "MarkerID", "SNPID", "SNP", "ID")
P_CANDIDATES = ("p.value", "P-value", "P", "pval", "PVALUE", "p_value")
# META_BETA is already on SAIGE Allele2. A raw step 2 file has BETA only.
BETA_CANDIDATES = ("META_BETA", "BETA")
SE_CANDIDATES = ("META_SE", "SE")
ALLELE2_CANDIDATES = ("SAIGE_Allele2", "Allele2")
ALLELE1_CANDIDATES = ("SAIGE_Allele1", "Allele1")

# AA_{gene}_{relative position}_{genomic bp}_{tail}
AA_FULL = re.compile(
    r"^AA[-_](?:HLA[-_])?([A-Z][A-Z0-9]*)_(-?\d+)_(\d+)_(.+)$"
)
# AA_{gene}_{relative position}_{residue}  — genomic bp taken from POS
AA_SHORT = re.compile(
    r"^AA[-_](?:HLA[-_])?([A-Z][A-Z0-9]*)_(-?\d+)_([A-Za-z*]+)$"
)
HLA_ALLELE = re.compile(r"^(?:HLA[-_])?([A-Z][A-Z0-9]*)\*(\S+)$")
RESIDUE = re.compile(r"^[A-Za-z*]{1,3}$")
EXON_TOKEN = re.compile(r"^exon\d+$", re.IGNORECASE)

ASSOC_HEADER = ["CHR", "SNP", "BP", "A1", "TEST", "NMISS", "OR", "STAT", "P"]


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


def norm_chrom(token):
    text = str(token).strip().lower()
    if text.startswith("chr"):
        text = text[3:]
    return text


def strip_novel(allele):
    for tag in (":new", "new"):
        if allele.endswith(tag):
            return allele[: -len(tag)]
    return allele


def pick(header, candidates):
    index = {name: i for i, name in enumerate(header)}
    for name in candidates:
        if name in index:
            return index[name]
    return None


def resolve_columns(header, path):
    if len(set(header)) != len(header):
        raise ValueError("duplicate column names in %s" % path)
    cols = {
        "marker": pick(header, MARKER_CANDIDATES),
        "p": pick(header, P_CANDIDATES),
        "beta": pick(header, BETA_CANDIDATES),
        "se": pick(header, SE_CANDIDATES),
        "allele2": pick(header, ALLELE2_CANDIDATES),
        "allele1": pick(header, ALLELE1_CANDIDATES),
        "pos": pick(header, ("POS", "BP")),
        "chr": pick(header, ("CHR", "CHROM", "#CHROM")),
        "n": pick(header, ("N",)),
        "n_case": pick(header, ("N_case", "N_cases", "N.Cases", "N.case")),
        "n_ctrl": pick(header, ("N_ctrl", "N_ctrls", "N.Controls", "N.ctrl")),
        "phecode": pick(header, ("phecode", "Phecode", "phenotype")),
    }
    missing = [
        what for what, key in (
            ("marker id", "marker"),
            ("p-value", "p"),
            ("BETA", "beta"),
            ("position", "pos"),
            ("effect allele", "allele2"),
        ) if cols[key] is None
    ]
    if missing:
        raise ValueError(
            "%s is missing %s; header is %s"
            % (path, ", ".join(missing), " ".join(header))
        )
    return cols


def split_row(line, use_tab):
    line = line.rstrip("\r\n")
    if use_tab:
        return line.split("\t")
    return line.split()


def hatk_snp(marker, pos):
    """Return (HATK SNP id, class). class is hla, aa, snp, or other."""
    raw = marker.strip()
    cleaned = strip_novel(raw)

    aa = AA_FULL.match(cleaned)
    if aa:
        gene, rel, bp, tail = aa.groups()
        return "AA_%s_%s_%s_%s" % (gene, rel, bp, tail), "aa"
    short = AA_SHORT.match(cleaned)
    if short and pos:
        gene, rel, residue = short.groups()
        return "AA_%s_%s_%s_%s" % (gene, rel, pos, residue), "aa"

    if re.match(r"^SNP(S)?_", cleaned) or re.match(r"^rs\d+$", cleaned):
        return cleaned, "snp"

    allele = HLA_ALLELE.match(cleaned)
    if allele:
        gene, fields = allele.group(1), strip_novel(allele.group(2))
        prefixed = raw.startswith("HLA-") or raw.startswith("HLA_")
        if prefixed or gene in HLA_GENES:
            return "HLA_%s*%s" % (gene, fields), "hla"
    return cleaned, "other"


def residue_token(snp):
    """Last field of an AA id when it is a residue rather than an exon label."""
    if not snp.startswith("AA_"):
        return None
    last = snp.rsplit("_", 1)[-1]
    if EXON_TOKEN.match(last) or not RESIDUE.match(last):
        return None
    return last


def effect_allele(snp, allele1, allele2, beta):
    """A1 and the beta that belongs to it.

    HATK signs a bi-allelic amino-acid column from A1, and a multi-allelic
    column from the last field of the marker. Both are set to the residue
    when the name carries one. BETA stays on that residue: flipped when
    SAIGE's Allele1, rather than Allele2, is the residue itself.
    """
    residue = residue_token(snp)
    if residue is None:
        return allele2, beta
    if allele1 is not None and allele1.upper() == residue.upper() and (
        allele2 is None or allele2.upper() != residue.upper()
    ):
        return residue, -beta
    return residue, beta


def format_or(beta):
    clipped = max(min(beta, 50.0), -50.0)
    return "%.6g" % math.exp(clipped)


def format_stat(beta, se):
    if se is None or se <= 0:
        return "NA"
    return "%.6g" % (beta / se)


def format_n(value):
    if value is None:
        return "NA"
    if abs(value - round(value)) < 1e-6:
        return str(int(round(value)))
    return "%.6g" % value


def format_p(token, value):
    if value == 0:
        return P_FLOOR
    return token


def read_saige(paths, chrom, phecode):
    """Rows ready for the PLINK table. First copy of a marker is kept."""
    rows = []
    seen = set()
    dropped = Counter()
    classes = Counter()
    phecodes = set()
    for path in paths:
        with open_text(path) as fh:
            header_line = None
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                header_line = line
                break
            if header_line is None:
                raise ValueError("no header in %s" % path)
            use_tab = "\t" in header_line
            header = split_row(header_line, use_tab)
            cols = resolve_columns(header, path)
            if phecode and cols["phecode"] is None:
                raise ValueError("%s has no phecode column" % path)
            need = max(i for i in cols.values() if i is not None)
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                tok = split_row(line, use_tab)
                if len(tok) <= need:
                    dropped["short_row"] += 1
                    continue
                if cols["phecode"] is not None:
                    this = tok[cols["phecode"]]
                    phecodes.add(this)
                    if phecode and this != phecode:
                        continue
                    if not phecode and len(phecodes) > 1:
                        raise ValueError(
                            "%s contains more than one phecode (%s). "
                            "Pass --phecode."
                            % (path, ", ".join(sorted(phecodes)[:8]))
                        )
                if chrom != "all" and cols["chr"] is not None:
                    if norm_chrom(tok[cols["chr"]]) != chrom:
                        dropped["chrom"] += 1
                        continue
                elif chrom != "all" and cols["chr"] is None:
                    dropped["chrom"] += 1
                    continue
                marker = tok[cols["marker"]]
                if marker in MISSING:
                    dropped["marker"] += 1
                    continue
                pos_text = tok[cols["pos"]]
                pos_value = parse_float(pos_text)
                if pos_value is None or pos_value <= 0:
                    dropped["pos"] += 1
                    continue
                pos = str(int(pos_value))
                beta = parse_float(tok[cols["beta"]])
                if beta is None:
                    dropped["beta"] += 1
                    continue
                p_token = tok[cols["p"]]
                p_value = parse_float(p_token)
                if p_value is None or p_value < 0 or p_value > 1:
                    dropped["p"] += 1
                    continue
                allele2 = tok[cols["allele2"]]
                allele1 = tok[cols["allele1"]] if cols["allele1"] is not None else ""
                if allele2 in MISSING:
                    dropped["allele"] += 1
                    continue
                se = parse_float(tok[cols["se"]]) if cols["se"] is not None else None
                if cols["n"] is not None:
                    n = parse_float(tok[cols["n"]])
                elif cols["n_case"] is not None and cols["n_ctrl"] is not None:
                    n_case = parse_float(tok[cols["n_case"]])
                    n_ctrl = parse_float(tok[cols["n_ctrl"]])
                    n = None if n_case is None or n_ctrl is None else n_case + n_ctrl
                else:
                    n = None
                snp, kind = hatk_snp(marker, pos)
                a1, beta_a1 = effect_allele(snp, allele1, allele2, beta)
                if snp in seen:
                    dropped["duplicate"] += 1
                    continue
                seen.add(snp)
                classes[kind] += 1
                chrom_out = "6"
                if cols["chr"] is not None:
                    chrom_out = norm_chrom(tok[cols["chr"]]).upper()
                rows.append({
                    "CHR": chrom_out,
                    "SNP": snp,
                    "BP": pos,
                    "A1": a1,
                    "TEST": "ADD",
                    "NMISS": format_n(n),
                    "OR": format_or(beta_a1),
                    "STAT": format_stat(beta_a1, se),
                    "P": format_p(p_token, p_value),
                    "kind": kind,
                    "beta": beta_a1,
                })
    rows.sort(key=lambda rec: (int(rec["BP"]), rec["SNP"]))
    return rows, dropped, classes


def write_assoc(rows, path):
    with open(path, "w", newline="") as fh:
        fh.write("\t".join(ASSOC_HEADER) + "\n")
        for rec in rows:
            fh.write("\t".join(rec[name] for name in ASSOC_HEADER) + "\n")


def gene_of(snp):
    if snp.startswith("HLA_"):
        return snp[4:].split("*", 1)[0]
    if snp.startswith("AA_"):
        body = snp[3:]
        return body.split("_", 1)[0]
    return None


def genes_present(rows):
    hla, aa = set(), set()
    for rec in rows:
        gene = gene_of(rec["SNP"])
        if gene is None:
            continue
        if rec["kind"] == "hla":
            hla.add(gene)
        elif rec["kind"] == "aa":
            aa.add(gene)
    return hla, aa


def order_genes(genes):
    rank = {gene: i for i, gene in enumerate(GENE_ORDER)}
    return sorted(genes, key=lambda gene: (rank.get(gene, len(GENE_ORDER)), gene))


def find_maptables(directory, gene, hg, imgt):
    if not directory or not os.path.isdir(directory):
        return []
    prefix = "HLA_MAPTABLE_%s" % gene
    found = []
    for name in sorted(os.listdir(directory)):
        if not name.startswith(prefix) or not name.endswith(".txt"):
            continue
        rest = name[len(prefix):]
        if rest and not rest.startswith("."):
            continue
        found.append(os.path.join(directory, name))
    if hg:
        tagged = [path for path in found if "hg%s" % hg in os.path.basename(path)]
        if tagged:
            found = tagged
    if imgt:
        tagged = [path for path in found if imgt in os.path.basename(path)]
        if tagged:
            found = tagged
    return found


def maptable_alleles(path):
    alleles = set()
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i < 3 or not line.strip():
                continue
            alleles.add(line.split()[0])
    return alleles


def maptable_relpos(path):
    with open(path) as fh:
        header = fh.readline().split()
    return set(header[1:])


def assoc_alleles(rows, gene):
    prefix = "HLA_%s*" % gene
    return {rec["SNP"][len("HLA_"):] for rec in rows if rec["SNP"].startswith(prefix)}


def assoc_relpos(rows, gene):
    positions = set()
    head = "AA_%s_" % gene
    for rec in rows:
        if not rec["SNP"].startswith(head):
            continue
        parts = rec["SNP"].split("_")
        # AA, gene, rel, bp, ...
        if len(parts) >= 4:
            positions.add(parts[2])
    return positions


def heatmap_jobs(rows, genes, maptable_dir, hg, imgt):
    """One (gene, maptable) pair per gene that can actually be drawn."""
    hla, aa = genes_present(rows)
    jobs = []
    notes = []
    for gene in order_genes(genes):
        if gene not in hla or gene not in aa:
            missing = []
            if gene not in hla:
                missing.append("HLA_%s* alleles" % gene)
            if gene not in aa:
                missing.append("AA_%s_ markers" % gene)
            notes.append("%s: no heatmap (%s)" % (gene, " or ".join(missing)))
            continue
        tables = find_maptables(maptable_dir, gene, hg, imgt)
        if not tables:
            notes.append("%s: no HLA_MAPTABLE_%s file in %s" % (gene, gene, maptable_dir))
            continue
        table = os.path.abspath(tables[0])
        if len(tables) > 1:
            notes.append("%s: using %s" % (gene, os.path.basename(table)))
        alleles = assoc_alleles(rows, gene)
        known = maptable_alleles(table)
        overlap = alleles & known
        if not overlap:
            notes.append(
                "%s: no allele in common with %s "
                "(assoc has %s; maptable rows look like %s)"
                % (
                    gene,
                    os.path.basename(table),
                    ", ".join(sorted(alleles)[:3]),
                    ", ".join(sorted(known)[:3]),
                )
            )
            continue
        rel = assoc_relpos(rows, gene)
        known_rel = maptable_relpos(table)
        if known_rel and not (rel & known_rel):
            notes.append(
                "%s: amino-acid positions %s are not in the maptable header"
                % (gene, ", ".join(sorted(rel)[:6]))
            )
            continue
        notes.append(
            "%s: %d/%d alleles match %s"
            % (gene, len(overlap), len(alleles), os.path.basename(table))
        )
        jobs.append((gene, table))
    return jobs, notes


def hatk_commands(python, assoc, out, hg, jobs, point_color, top_color, point_size, yaxis_unit):
    commands = []
    manhattan = [
        python, "-m", "HLA_Manhattan",
        "--assoc-result", assoc,
        "--out", out + ".manhattan",
        "--hg", hg,
    ]
    if point_color:
        manhattan.extend(["--point-color", point_color])
    if top_color:
        manhattan.extend(["--top-color", top_color])
    if point_size:
        manhattan.extend(["--point-size", point_size])
    if yaxis_unit:
        manhattan.extend(["--yaxis-unit", yaxis_unit])
    commands.append(manhattan)
    for gene, table in jobs:
        commands.append([
            python, "-m", "HLA_Heatmap",
            "--HLA", gene,
            "--maptable", table,
            "--assoc-result", assoc,
            "--out", "%s.heatmap.%s" % (out, gene),
        ])
    return commands


def run_commands(commands, hatk_dir):
    failed = 0
    for command in commands:
        print(" ".join(command), file=sys.stderr)
        try:
            subprocess.run(command, cwd=hatk_dir, check=True)
        except subprocess.CalledProcessError as exc:
            failed += 1
            print("HATK exited %s" % exc.returncode, file=sys.stderr)
    return failed


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--self-test", action="store_true", help="run the built-in checks and exit")
    parser.add_argument("--saige", nargs="+", help="SAIGE step 2 or saige_metal.py table(s), optionally gzipped")
    parser.add_argument("--out", help="output prefix; writes {out}.assoc.logistic and the PDFs beside it")
    parser.add_argument("--hatk", help="clone of https://github.com/WansonChoi/HATK")
    parser.add_argument("--python", default="python3",
                        help="interpreter that can import HATK (conda env HATK)")
    parser.add_argument("--hg", default="38", choices=["18", "19", "38"],
                        help="genome build of the POS column (default 38)")
    parser.add_argument("--chrom", default="6",
                        help="chromosome to keep (default 6). 'all' keeps every chromosome")
    parser.add_argument("--phecode", help="meta table row to plot, when the file has several")
    parser.add_argument("--hla", nargs="+",
                        help="heatmap genes, as HATK spells them (A DRB1). Default: every gene with both alleles and amino acids")
    parser.add_argument("--maptable-dir",
                        help="directory of HLA_MAPTABLE_{gene}.txt files. Default: {hatk}/HLA_Heatmap/data")
    parser.add_argument("--imgt", help="prefer a maptable whose file name contains this version")
    parser.add_argument("--no-manhattan", action="store_true")
    parser.add_argument("--no-heatmap", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="write the association table and print the HATK commands")
    parser.add_argument("--point-color", help="Manhattan point color, e.g. #778899")
    parser.add_argument("--top-color", help="Manhattan top-hit color, e.g. #FF0000")
    parser.add_argument("--point-size", help="Manhattan point size passed to HATK (default 15)")
    parser.add_argument("--yaxis-unit", help="Manhattan -log10 tick spacing passed to HATK (default 10)")
    return parser


def plot(args):
    if not args.saige or not args.out:
        raise SystemExit("--saige and --out are required")
    if args.no_manhattan and args.no_heatmap:
        raise SystemExit("both plots were disabled")
    chrom = norm_chrom(args.chrom) if args.chrom.lower() != "all" else "all"
    rows, dropped, classes = read_saige(args.saige, chrom, args.phecode)
    if not rows:
        raise SystemExit("no usable markers (%s)" % dict(dropped))
    if args.hatk:
        args.hatk = os.path.abspath(args.hatk)
    if os.sep in args.python:
        args.python = os.path.abspath(args.python)
    out = os.path.abspath(args.out)
    out_dir = os.path.dirname(out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    assoc = out + ".assoc.logistic"
    write_assoc(rows, assoc)
    print(
        "wrote %d markers to %s (%s); dropped %s"
        % (len(rows), assoc, ", ".join("%s=%d" % item for item in sorted(classes.items())), dict(dropped)),
        file=sys.stderr,
    )

    hla, aa = genes_present(rows)
    requested = order_genes(hla & aa) if not args.hla else order_genes(
        norm_gene(gene) for gene in args.hla
    )
    maptable_dir = args.maptable_dir
    if maptable_dir is None and args.hatk:
        maptable_dir = os.path.join(args.hatk, "HLA_Heatmap", "data")
    jobs, notes = ([], [])
    if not args.no_heatmap:
        jobs, notes = heatmap_jobs(rows, requested, maptable_dir, args.hg, args.imgt)
        for note in notes:
            print(note, file=sys.stderr)

    if not args.dry_run:
        if not args.hatk or not os.path.isdir(os.path.join(args.hatk, "HLA_Manhattan")):
            raise SystemExit("--hatk must be a HATK checkout (it contains HLA_Manhattan/)")
    commands = hatk_commands(
        args.python, assoc, out, args.hg,
        [] if args.no_heatmap else jobs,
        args.point_color, args.top_color, args.point_size, args.yaxis_unit,
    )
    if args.no_manhattan:
        commands = [cmd for cmd in commands if "HLA_Manhattan" not in cmd]
    if args.dry_run or not commands:
        for command in commands:
            print(" ".join(command))
        return 0
    return run_commands(commands, os.path.abspath(args.hatk))


def norm_gene(gene):
    text = gene.strip()
    for prefix in ("HLA-", "HLA_"):
        if text.startswith(prefix):
            text = text[len(prefix):]
    return text


def self_test():
    import tempfile

    assert hatk_snp("HLA-DRB1*04:01:new", "1") == ("HLA_DRB1*04:01", "hla")
    assert hatk_snp("HLA_A*02:01", "1") == ("HLA_A*02:01", "hla")
    assert hatk_snp("A*02:01", "1") == ("HLA_A*02:01", "hla")
    assert hatk_snp("KIR2DL1*001", "1") == ("KIR2DL1*001", "other")
    assert hatk_snp("AA_HLA-DRB1_71_32659900_K", "9") == ("AA_DRB1_71_32659900_K", "aa")
    assert hatk_snp("AA_DRB1_71_K", "32659900") == ("AA_DRB1_71_32659900_K", "aa")
    assert hatk_snp("AA_A_-15_29910359_exon1", "1") == ("AA_A_-15_29910359_exon1", "aa")
    assert hatk_snp("AA_A_9_29910558_exon2_F", "1") == ("AA_A_9_29910558_exon2_F", "aa")
    assert hatk_snp("SNPS_A_30018457", "1") == ("SNPS_A_30018457", "snp")
    assert hatk_snp("rs41557221", "1") == ("rs41557221", "snp")
    assert residue_token("AA_A_9_29910558_exon2_F") == "F"
    assert residue_token("AA_A_-15_29910359_exon1") is None

    a1, beta = effect_allele("AA_DRB1_71_32659900_K", "A", "G", 0.5)
    assert a1 == "K" and beta == 0.5
    a1, beta = effect_allele("AA_DRB1_71_32659900_K", "K", "G", 0.5)
    assert a1 == "K" and beta == -0.5
    assert abs(float(format_or(math.log(2))) - 2) < 1e-6

    with tempfile.TemporaryDirectory() as tmp:
        saige = os.path.join(tmp, "step2_singlevar")
        with open(saige, "w") as fh:
            fh.write(
                "CHR POS MarkerID Allele1 Allele2 AC_Allele2 AF_Allele2 "
                "MissingRate BETA SE Tstat var p.value N\n"
            )
            fh.write("6 29910558 AA_A_9_29910558_exon2_F A G 10 0.05 0 0.693147 0.1 1 1 1e-8 1000\n")
            fh.write("6 29910247 HLA-DRB1*04:01 A G 10 0.05 0 -0.693147 0.2 1 1 0.01 1000\n")
            fh.write("6 29910247 HLA-A*02:01 A G 10 0.05 0 0.1 0.2 1 1 0 1000\n")
            fh.write("6 29910248 HLA-A*02:01 A G 10 0.05 0 0.9 0.2 1 1 0.2 1000\n")
            fh.write("1 100 rs1 A T 10 0.05 0 0.2 0.1 1 1 0.04 1000\n")
            fh.write("6 29911000 AA_DRB1_71_R R G 10 0.05 0 0.4 0.1 1 1 0.03 800\n")
            fh.write("6 50 bad A G 10 0.05 0 NA 0.1 1 1 0.5 1000\n")
        rows, dropped, classes = read_saige([saige], "6", None)
        by_snp = {rec["SNP"]: rec for rec in rows}
        assert set(by_snp) == {
            "AA_A_9_29910558_exon2_F",
            "HLA_DRB1*04:01",
            "HLA_A*02:01",
            "AA_DRB1_71_29911000_R",
        }, set(by_snp)
        assert by_snp["AA_A_9_29910558_exon2_F"]["A1"] == "F"
        assert abs(float(by_snp["AA_A_9_29910558_exon2_F"]["OR"]) - 2) < 1e-4
        assert by_snp["HLA_A*02:01"]["P"] == P_FLOOR
        assert by_snp["AA_DRB1_71_29911000_R"]["BP"] == "29911000"
        assert by_snp["AA_DRB1_71_29911000_R"]["A1"] == "R"
        assert abs(float(by_snp["AA_DRB1_71_29911000_R"]["beta"]) - (-0.4)) < 1e-9
        assert dropped["chrom"] == 1 and dropped["duplicate"] == 1 and dropped["beta"] == 1
        assert classes["hla"] == 2 and classes["aa"] == 2

        meta = os.path.join(tmp, "meta.tsv")
        with open(meta, "w") as fh:
            fh.write(
                "phecode\tCHR\tPOS\tMarkerName\tSAIGE_Allele1\tSAIGE_Allele2\t"
                "META_BETA\tMETA_SE\tP-value\tN\tAllele1\tBETA\n"
            )
            fh.write("P1\t6\t29910247\tHLA-B*07:02\tA\tG\t0.2\t0.1\t0.04\t1000\tp\t9\n")
            fh.write("P1\t6\t29910558\tAA_B_9_29910558_exon2_Y\tA\tG\t\t0.1\t0.5\t1000\tp\t9\n")
            fh.write("P2\t6\t1\tHLA-B*08:01\tA\tG\t0.2\t0.1\t0.04\t1000\tp\t9\n")
        try:
            read_saige([meta], "6", None)
        except ValueError as exc:
            assert "phecode" in str(exc)
        else:
            raise AssertionError("two phecodes should be rejected")
        meta_rows, _, _ = read_saige([meta], "6", "P1")
        meta_by = {rec["SNP"]: rec for rec in meta_rows}
        assert set(meta_by) == {"HLA_B*07:02"}
        assert abs(float(meta_by["HLA_B*07:02"]["OR"]) - math.exp(0.2)) < 1e-4

        table_dir = os.path.join(tmp, "maps")
        os.makedirs(table_dir)
        for name, alleles, rels in (
            ("HLA_MAPTABLE_A.hg19.imgt3320.txt", ["A*02:01"], ["9"]),
            ("HLA_MAPTABLE_A.hg38.imgt3700.txt", ["A*02:01", "A*01:01"], ["9", "-15"]),
            ("HLA_MAPTABLE_DMA.hg38.imgt3700.txt", ["DMA*01:01"], ["1"]),
            ("HLA_MAPTABLE_DRB1.hg38.imgt3700.txt", ["DRB1*01:01"], ["1"]),
        ):
            with open(os.path.join(table_dir, name), "w") as fh:
                fh.write("AA_rel_pos\t" + "\t".join(rels) + "\n")
                fh.write("genomic_pos\t" + "\t".join("1" for _ in rels) + "\n")
                fh.write("exon\t" + "\t".join("exon2" for _ in rels) + "\n")
                for allele in alleles:
                    fh.write(allele + "\t" + "\t".join("K" for _ in rels) + "\n")
        chosen = find_maptables(table_dir, "A", "38", None)
        assert [os.path.basename(path) for path in chosen] == ["HLA_MAPTABLE_A.hg38.imgt3700.txt"]
        jobs, notes = heatmap_jobs(rows, ["A", "DRB1", "C"], table_dir, "38", None)
        assert [gene for gene, _path in jobs] == ["A"]
        assert any(note.startswith("DRB1:") for note in notes)
        assert any(note.startswith("C:") for note in notes)

        assoc = os.path.join(tmp, "out.assoc.logistic")
        write_assoc(rows, assoc)
        with open(assoc) as fh:
            header = fh.readline().rstrip("\n").split("\t")
        assert header == ASSOC_HEADER
        commands = hatk_commands(
            "python3", assoc, os.path.join(tmp, "out"), "38", jobs, None, None, None, None)
        assert commands[0][2] == "HLA_Manhattan"
        assert commands[0][commands[0].index("--assoc-result") + 1].endswith(".assoc.logistic")
        assert commands[1][2] == "HLA_Heatmap"
        assert "A" in commands[1]

    print("self-test ok")


def main(argv):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.self_test:
        self_test()
        return 0
    return plot(args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]) or 0)
