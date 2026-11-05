"""
Publication-style plot of gene order / copy-number structure across haplotypes
(MHC or KIR), coloured by functional class from a gene table.

Input: ONE annotation table for all assemblies, columns like
    Contig  Tool  Category  Spos  Epos  Info  Strand  Phase  Attribute  Gene  Sample
(by default each contig is one haplotype, so no haplotype id is needed;
 use bounds=("HLA-F", "HLA-A") so only contigs carrying the whole region are counted).

One call:
    plot_gene_structure(annot_df, "HLA_classI", hap_number=21, genes=CLASS_I_GENES,
                        bounds=("HLA-F", "HLA-A"))
or step by step:
    feat, haps = collect_haplotype_features(annot_df, genes=CLASS_I_GENES, hap_number=21)
    plot_haplotype_structures(feat, "/home/jupyter/Gene_Structure/HLA_classI")

Several structures in one figure (subplots)
    plot_haplotype_structures draws into an axes you pass, so panels can be composed:
        feats = [collect_haplotype_features(df, genes=G, hap_number=8)[0] for G in gene_sets]
        fig, axes = make_panel_figure([f["rank"].nunique() for f in feats])
        for a, f, t in zip(axes, feats, ["a  MHC class II", "b  class I", "c  KIR"]):
            plot_haplotype_structures(f, ax=a, title=t, freq_xmax=10000)
        fig.savefig("panels.pdf", facecolor="white", dpi=600)
    make_panel_figure sizes the panels in mm so a haplotype row is the same height in
    every panel; any axes works, though (plt.subplots, GridSpec), in which case the rows
    share whatever height that axes has. Pass show_freq=False for a panel without the
    count bars, legend=False to draw the legend only once, and freq_xmax= to put every
    panel's bars on the same count scale. Don't use tight_layout / constrained_layout:
    the frequency panel is placed in mm inside the axes box.

Ancestry stratification (optional)
    Pass a DataFrame with columns "research_id" (matched to the annotation table's Sample
    column) and "ancestry_pred" and the count bars are split by ancestry:
        plot_gene_structure(annot_df, "HLA_classI", hap_number=21, genes=CLASS_I_GENES,
                            bounds=("HLA-F", "HLA-A"), ancestry_df=ancestry_df)
    Each haplotype is assigned the ancestry of its sample (so both haplotypes of one
    person count once each). Samples without a label are shown as "Unknown".

Colouring
    functional genes                      -> palette colours (default Set2), drawn larger
    pseudogene / fragment / not_expressed -> grey, drawn smaller
    genes missing from the table          -> light grey "Unclassified" (a warning lists them)
    If there are more functional genes than palette colours (e.g. KIR), genes are
    coloured by the table's `group` column instead (inhibitory / activating / ...).

Output: <prefix>.pdf (vector, editable text), <prefix>.png (600 dpi),
        <prefix>_source_data.tsv (gene coordinates, class and colour per structure),
        <prefix>_haplotype_assignments.tsv (every haplotype's structure, rank and n),
        <prefix>_ancestry_counts.tsv (structure rank x ancestry counts; only with ancestry).
"""
from __future__ import annotations

import os
import re

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_hex, to_rgb
from matplotlib.patches import Patch, Polygon
from matplotlib.ticker import FixedLocator, MultipleLocator, NullFormatter

MM = 1 / 25.4  # mm -> inch

GENE_TABLE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hla_kir_gene_classes.tsv")
NONFUNCTIONAL = {"pseudogene", "fragment", "not_expressed"}

PSEUDO_COLOR = "#6f6d68"     # pseudogenes / fragments
UNKNOWN_COLOR = "#cfcdc7"    # genes not in the table

INK = "#1f1f1d"
INK_SECONDARY = "#5c5b56"
BACKBONE = "#c2c0ba"
GRID = "#e4e3df"
ZEBRA = "#f4f4f1"
BAR = "#6b6964"

# Ancestry colours (Okabe-Ito, colour-blind safe), keyed by lower-case label. The order
# below is also the stacking / legend order; unlisted labels follow, then other, Unknown.
ANCESTRY_ORDER = ["afr", "amr", "eas", "eur", "mid", "sas"]
ANCESTRY_COLORS = {
    "afr": "#E69F00", "amr": "#CC79A7", "eas": "#009E73",
    "eur": "#0072B2", "mid": "#56B4E9", "sas": "#D55E00",
    "other": "#8c8a85", "unknown": UNKNOWN_COLOR,
}
ANCESTRY_FALLBACK = ["#7570b3", "#a6761d", "#66a61e", "#e7298a", "#1b9e77", "#666666"]

PUB_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.linewidth": 0.6,
    "axes.edgecolor": INK_SECONDARY,
    "axes.labelcolor": INK,
    "xtick.color": INK_SECONDARY,
    "xtick.labelcolor": INK,
    "xtick.labelsize": 6.5,
    "xtick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.labelsize": 7,
    "ytick.major.size": 0,
    "pdf.fonttype": 42,   # TrueType: text stays editable in Illustrator
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.dpi": 600,
    "mathtext.fontset": "custom",
    "mathtext.rm": "sans",
    "mathtext.it": "sans:italic",
    "mathtext.bf": "sans:bold",
    "mathtext.cal": "sans",
}


# --------------------------------------------------------------------------- #
# 0. Gene table                                                                #
# --------------------------------------------------------------------------- #
def load_gene_table(path=GENE_TABLE):
    """Read the gene class TSV (columns: gene, region, class, group, note)."""
    t = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = {"gene", "class", "group"} - set(t.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    t["class"] = t["class"].str.strip().str.lower()
    t.index = t["gene"].str.strip().str.upper()
    return t


def _lookup(gene, table):
    """Find a gene in the table, tolerating 'A' vs 'HLA-A' and KIR2DL5 vs KIR2DL5A/B."""
    g = str(gene).strip().upper()
    for cand in (g, f"HLA-{g}", re.sub(r"[AB]$", "", g), f"{g}A"):
        if cand in table.index:
            return table.loc[cand]
    return None


def _palette_colors(palette):
    """Colours from a matplotlib colormap name or a list; neutral greys are dropped
    (e.g. Set2's last slot) so no gene colour can be mistaken for a pseudogene."""
    if isinstance(palette, str):
        cmap = matplotlib.colormaps[palette]
        cols = getattr(cmap, "colors", None)
        if cols is None or len(cols) > 20:  # continuous map: sample it
            cols = cmap(np.linspace(0.05, 0.95, 8))
    else:
        cols = palette
    out = []
    for c in cols:
        r, g, b = to_rgb(c)
        if max(r, g, b) - min(r, g, b) > 0.08:
            out.append(to_hex(c))
    return out


def classify_and_color(feat, gene_table=GENE_TABLE, palette="Set2", color_by="auto",
                       colors=None, ref_rank=1):
    """Add class / group / color columns to `feat`; return (feat, legend_entries)."""
    table = gene_table if isinstance(gene_table, pd.DataFrame) else load_gene_table(gene_table)
    feat = feat.copy()

    info = {}
    for g in feat["gene"].unique():
        row = _lookup(g, table)
        info[g] = ("unknown", "", "") if row is None else (row["class"], row["group"], row.get("note", ""))
    feat["class"] = feat["gene"].map(lambda g: info[g][0])
    feat["group"] = feat["gene"].map(lambda g: info[g][1])
    unknown = sorted(g for g, v in info.items() if v[0] == "unknown")
    if unknown:
        print(f"Not in gene table (drawn as 'Unclassified'): {unknown}")

    # colour order = order along the reference haplotype, then any others
    order = pd.concat([feat[feat["rank"] == ref_rank].sort_values("start"),
                       feat.sort_values(["rank", "start"])])
    func_genes = [g for g in dict.fromkeys(order["gene"]) if info[g][0] == "functional"]
    func_groups = list(dict.fromkeys(info[g][1] for g in func_genes))
    pal = _palette_colors(palette)

    if color_by == "auto":
        color_by = "gene" if len(func_genes) <= len(pal) else "group"
        if color_by == "group":
            print(f"{len(func_genes)} functional genes but {len(pal)} colours in palette "
                  f"'{palette}': colouring by group ({', '.join(func_groups)}).")
    keys = func_genes if color_by == "gene" else func_groups
    if len(keys) > len(pal):
        raise ValueError(f"{len(keys)} {color_by}s to colour but palette has {len(pal)} colours; "
                         "pass a larger qualitative palette (e.g. 'tab10') or a list of colours.")
    key_color = dict(zip(keys, pal))
    if colors:  # manual overrides, keyed by gene (or group) name
        key_color.update({k: to_hex(v) for k, v in colors.items()})

    def color_of(g):
        cls, grp, _ = info[g]
        if cls == "functional":
            return key_color.get(g if color_by == "gene" else grp, PSEUDO_COLOR)
        return UNKNOWN_COLOR if cls == "unknown" else PSEUDO_COLOR

    feat["color"] = feat["gene"].map(color_of)

    legend = [(k if color_by == "gene" else k[:1].upper() + k[1:], key_color[k], color_by == "gene")
              for k in keys]
    classes = set(feat["class"])
    if classes & NONFUNCTIONAL:
        lab = "Pseudogene / fragment" if "fragment" in classes else "Pseudogene"
        legend.append((lab, PSEUDO_COLOR, False))
    if unknown:
        legend.append(("Unclassified", UNKNOWN_COLOR, False))
    return feat, legend


# --------------------------------------------------------------------------- #
# 1. Data: one annotation table -> haplotype structures -> tidy gene rows      #
# --------------------------------------------------------------------------- #
ANNOT_COLS = ("Contig", "Spos", "Epos", "Strand", "Gene")


def _strip(name, prefixes):
    for p in prefixes:
        if name.startswith(p):
            return name[len(p):]
    return name


def _hap_from_contig(contigs):
    """Haplotype number from contig names: hifiasm 'h1tg000350l', PanSN 'S#1#ctg',
    or names containing 'hap1' / 'haplotype2'."""
    c = contigs.astype(str)
    hap = c.str.extract(r"^h(\d+)tg", expand=False)
    for pat in (r"#(\d+)#", r"(?i)hap(?:lotype)?[_\-.]?(\d+)"):
        miss = hap.isna()
        if miss.any():
            hap[miss] = c[miss].str.extract(pat, expand=False)
    if hap.isna().any():
        bad = c[hap.isna()].unique()[:5].tolist()
        raise ValueError(f"Can't tell the haplotype from contig names such as {bad}; "
                         "use hap_col=None (each contig = one haplotype) or a column name.")
    return hap


def find_haplotype_structures(annot, genes=None, region=None, *,
                              feature=("gene", "pseudogene"),
                              gene_table=GENE_TABLE, sample_col="Sample", hap_col=None,
                              bounds=None, crop_to_bounds=False,
                              orient_by=None, orient_strand="+", fragmented="drop"):
    """Read every haplotype's gene order from one annotation table.

    annot      : one table for all assemblies with columns Contig, Category, Spos, Epos,
                 Strand, Gene, Sample (e.g. concatenated Immuannot output).
    genes      : genes that define a structure, e.g. the class I segment HLA-F ... HLA-J.
    region     : instead of `genes`: "MHC" or "KIR" = every gene of that region in the
                 gene table. With neither, every gene in `annot` is used.
    feature    : Category value(s) to keep (default "gene" and "pseudogene"); None keeps
                 all rows. Exon / mRNA / CDS rows must be excluded or genes count twice.
    hap_col    : how contigs are grouped into haplotypes.
                 None (default) = no haplotype id: each contig is its own unit
                 (unit id "<Sample>_<Contig>"). With bounds= this counts each contig that
                 carries the complete region once, i.e. once per haplotype.
                 "from_contig" = parse hifiasm-style names (h1tg..., h2tg...).
                 any other string = name of a column holding the haplotype.
    bounds     : (first_gene, last_gene), e.g. ("HLA-F", "HLA-A"): keep only haplotypes with
                 both genes on the same contig (complete annotation); others are left out.
                 Each kept haplotype is oriented so first_gene is on the left.
    crop_to_bounds : also drop genes outside the first_gene..last_gene interval.
    orient_by  : (when no bounds) gene used to orient each haplotype: contigs are flipped
                 so this gene reads on `orient_strand`. Default = the gene present in the
                 most haplotypes.
    fragmented : haplotypes whose genes sit on more than one contig:
                 "drop" (default; their gene order can't be read reliably),
                 "largest" (keep only the contig carrying most of the genes), or
                 "keep" (not recommended: positions from different contigs get mixed).

    Returns (rows, haps): `rows` = one oriented row per gene per haplotype;
    `haps` = one row per unit (hap_id = "<Sample>_<Contig>", or "<Sample>_<hap>" when
    hap_col is given) with its structure,
    the comma-joined gene order read left to right.
    """
    missing = [c for c in ANNOT_COLS + (sample_col,) if c not in annot.columns]
    if missing:
        raise ValueError(f"annotation table is missing columns: {missing}")
    d_all = annot
    d = annot
    if feature is not None and "Category" in d.columns:
        feats = {feature} if isinstance(feature, str) else set(feature)
        d = d[d["Category"].isin(feats)]
    d = d.copy()
    d["gene"] = d["Gene"].astype(str).str.split("*").str[0]

    if genes is not None:
        gene_set = set(genes)
    else:
        table = gene_table if isinstance(gene_table, pd.DataFrame) else load_gene_table(gene_table)
        if region is not None:
            gene_set = set(table.loc[table["region"].str.upper() == region.upper(), "gene"])
        else:
            gene_set = set(d["gene"])
            regions = {table.loc[g.upper(), "region"] for g in gene_set if g.upper() in table.index}
            if len(regions) > 1:
                raise ValueError(f"The table has genes from {sorted(regions)}; "
                                 "pass region='MHC' / 'KIR' or a genes= list.")
    if bounds is not None:
        if len(bounds) != 2:
            raise ValueError("bounds must be (first_gene, last_gene)")
        gene_set |= set(bounds)
    # report requested genes that never appear (often a different Category label)
    if genes is not None or bounds is not None:
        absent = sorted(gene_set - set(d["gene"]))
        if absent:
            gall = d_all["Gene"].astype(str).str.split("*").str[0]
            other = d_all.loc[gall.isin(absent)]
            msg = f"Requested genes not found in the table: {absent}"
            if len(other) and "Category" in other.columns:
                msg += (f" -- they do appear with Category {sorted(other['Category'].unique())};"
                        " pass feature=[...] to include them")
            print(msg)
    d = d[d["gene"].isin(gene_set)]
    if d.empty:
        raise ValueError("None of the requested genes are in the annotation table.")

    d["sample"] = d[sample_col].astype(str)
    if hap_col is None:
        d["hap"] = d["Contig"].astype(str)
        unit = "contigs"
    elif hap_col == "from_contig":
        d["hap"] = _hap_from_contig(d["Contig"])
        unit = "haplotypes"
    else:
        d["hap"] = d[hap_col].astype(str)
        unit = "haplotypes"
    d["hap_id"] = d["sample"] + "_" + d["hap"]
    d["Spos"] = d["Spos"].astype(int)
    d["Epos"] = d["Epos"].astype(int)

    # complete annotation: both bound genes on one contig
    bound_flip = None
    if bounds is not None:
        b0, b1 = bounds
        key = ["hap_id", "Contig"]
        g0, g1 = d[d["gene"] == b0].groupby(key), d[d["gene"] == b1].groupby(key)
        st = pd.concat({"s0": g0["Spos"].min(), "e0": g0["Epos"].max(),
                        "s1": g1["Spos"].min(), "e1": g1["Epos"].max()}, axis=1).dropna()
        st = st[~st.index.get_level_values("hap_id").duplicated()].reset_index()
        fwd = st["s0"] <= st["s1"]
        st["lo"] = np.where(fwd, st["s0"], st["s1"])
        st["hi"] = np.where(fwd, st["e1"], st["e0"])
        bound_flip = pd.Series(~fwd.to_numpy(), index=st["hap_id"])
        n_all = d["hap_id"].nunique()
        where = "" if unit == "contigs" else " on one contig"
        print(f"{n_all - len(st)} of {n_all} {unit} don't have both {b0} and {b1}{where} "
              f"(incomplete; left out).")
        if crop_to_bounds:
            d = d.merge(st[key + ["lo", "hi"]], on=key)
            d = d[(d["Spos"] >= d["lo"]) & (d["Epos"] <= d["hi"])].drop(columns=["lo", "hi"])
        else:
            d = d[d["hap_id"].isin(st["hap_id"])]

    # haplotypes split over several contigs
    n_contig = d.groupby("hap_id")["Contig"].nunique()
    frag = n_contig.index[n_contig > 1]
    if len(frag):
        if fragmented == "drop":
            d = d[~d["hap_id"].isin(frag)]
            what = "left out"
        elif fragmented == "largest":
            main = (d.groupby(["hap_id", "Contig"]).size().reset_index(name="k")
                    .sort_values("k").drop_duplicates("hap_id", keep="last"))
            d = d.merge(main[["hap_id", "Contig"]], on=["hap_id", "Contig"])
            what = "kept the contig with the most genes"
        else:
            what = "kept as is; positions mix contigs"
        print(f"{len(frag)} haplotypes have these genes on more than one contig ({what}).")

    # without bounds: flag haplotypes missing genes that nearly every haplotype has
    if bounds is None:
        n_h = d["hap_id"].nunique()
        frac = d.groupby("gene")["hap_id"].nunique() / n_h
        core = sorted(frac.index[frac >= 0.9])
        n_full = (d[d["gene"].isin(core)].groupby("hap_id")["gene"].nunique() == len(core)).sum()
        if core and n_h - n_full:
            print(f"Note: {n_h - n_full} of {n_h} {unit} lack at least one of {core} "
                  "(each found in >=90% of haplotypes), so they are probably incomplete "
                  "annotations. To keep only complete ones, pass bounds=(first_gene, last_gene), "
                  "e.g. bounds=('HLA-F', 'HLA-A').")

    # orientation: bounds[0] on the left, else flip so the anchor gene reads on orient_strand
    if bound_flip is not None:
        orient_by = bounds[0]
        d["flipped"] = d["hap_id"].map(bound_flip)
    else:
        if orient_by is None:
            orient_by = d.groupby("gene")["hap_id"].nunique().idxmax()
        anchor = d[d["gene"] == orient_by].groupby("hap_id")["Strand"].first()
        majority = d.groupby(["hap_id", "Strand"]).size().unstack(fill_value=0).idxmax(axis=1)
        hap_strand = anchor.reindex(majority.index).fillna(majority)
        d["flipped"] = d["hap_id"].map(hap_strand != orient_strand)

    lo = d.groupby("hap_id")["Spos"].transform("min")
    hi = d.groupby("hap_id")["Epos"].transform("max")
    f = d["flipped"].to_numpy(bool)
    d["start"] = np.where(f, hi - d["Epos"], d["Spos"] - lo)
    d["end"] = np.where(f, hi - d["Spos"], d["Epos"] - lo)
    d["strand"] = np.where(f, d["Strand"].map({"+": "-", "-": "+"}).fillna(d["Strand"]), d["Strand"])
    d["hap_length"] = (hi - lo).astype(int)
    d = d.sort_values(["hap_id", "start", "end"]).reset_index(drop=True)

    haps = d.groupby("hap_id", sort=False).agg(
        sample=("sample", "first"), hap=("hap", "first"), contig=("Contig", "first"),
        flipped=("flipped", "first"), n_genes=("gene", "size"), structure=("gene", ",".join))
    haps.attrs["orient_by"] = orient_by
    return d, haps


def _select_ranks(n, hap_number=None, ranks=None):
    """Structure ranks (1 = most common) to plot. `ranks` wins over `hap_number`.
    hap_number: N -> top N; -N -> the N rarest; None -> all.
    ranks: explicit list; negative values count from the end (-1 = rarest)."""
    if ranks is not None:
        out = [n + 1 + r if r < 0 else r for r in ranks]
        return [r for r in dict.fromkeys(out) if 1 <= r <= n]
    if not hap_number:
        return list(range(1, n + 1))
    if hap_number < 0:
        return list(range(max(1, n + 1 + hap_number), n + 1))
    return list(range(1, min(hap_number, n) + 1))


def _order_ancestry(labels):
    """Stable display order: known labels, any others (alphabetical), other, Unknown."""
    def key(lab):
        k = str(lab).lower()
        if k in ANCESTRY_ORDER:
            return (0, ANCESTRY_ORDER.index(k), "")
        return {"other": (2, 0, ""), "unknown": (3, 0, "")}.get(k, (1, 0, str(lab)))
    return sorted(labels, key=key)


def _ancestry_of(samples, ancestry_df, id_col, anc_col):
    """Ancestry label for each entry of `samples` (matched on id_col, as strings);
    anything unmatched or blank becomes "Unknown"."""
    missing = [c for c in (id_col, anc_col) if c not in ancestry_df.columns]
    if missing:
        raise ValueError(f"ancestry_df is missing columns: {missing}")
    a = ancestry_df[[id_col, anc_col]].dropna().astype(str)
    a[id_col] = a[id_col].str.strip()
    a[anc_col] = a[anc_col].str.strip()
    a = a[a[anc_col] != ""]
    n_lab = a.groupby(id_col)[anc_col].nunique()
    if (n_lab > 1).any():
        print(f"{int((n_lab > 1).sum())} {id_col} values have more than one {anc_col}; "
              "using the first.")
    mapping = a.drop_duplicates(id_col).set_index(id_col)[anc_col]
    out = samples.astype(str).str.strip().map(mapping).astype(object).fillna("Unknown")
    n_miss = samples[out.eq("Unknown")].nunique()
    if n_miss:
        print(f"{n_miss} of {samples.nunique()} samples have no ancestry label "
              "(shown as 'Unknown').")
    return out


def _ancestry_counts_by_rank(counts, anc_of, ranks):
    """Table (index = structure rank, columns = ancestry) of haplotype counts."""
    rows = {}
    for r in ranks:
        members = counts[r - 1][1]
        rows[r] = pd.Series([anc_of.get(m, "Unknown") for m in members]).value_counts()
    tab = pd.DataFrame(rows).T.fillna(0).astype(int)
    tab.index.name = "rank"
    return tab[_order_ancestry(list(tab.columns))]


def collect_haplotype_features(annot, genes=None, region=None, hap_number=None, ranks=None,
                               exclude_ranks=(), Genes_dict=None, label_col="Gene",
                               strip_prefix=("HLA-", "KIR"), ancestry_df=None,
                               ancestry_id_col="research_id", ancestry_col="ancestry_pred",
                               **find_kws):
    """Group haplotypes by gene order, rank structures by count, and return the gene
    coordinates of one representative assembly per structure.

    annot, genes, region, **find_kws : see find_haplotype_structures.
    hap_number    : top N structures; a negative number = the N rarest (None = all).
    ranks         : explicit structure ranks to plot instead, e.g. [1, -5, -4, -3, -2, -1]
                    (negative = counted from the rarest end). Rank numbers are kept.
    exclude_ranks : structure ranks to leave out of the figure (numbering is kept).
    Genes_dict    : optional precomputed groups {structure key: [hap_id, ...]} with ids in
                    the same form as find_haplotype_structures makes them;
                    if given, its groups and counts are used and `annot` only supplies
                    coordinates.
    label_col     : column used for the glyph labels (prefixes in strip_prefix removed).
    ancestry_df   : optional DataFrame with columns ancestry_id_col ("research_id", matched
                    to the annotation table's Sample column) and ancestry_col
                    ("ancestry_pred"). Adds an `ancestry` column to `haps` and a rank x
                    ancestry count dict in feat.attrs["ancestry_counts"].

    Returns (feat, haps): `feat` = one row per gene per plotted structure (input to
    plot_haplotype_structures); `haps` = every haplotype with its structure, rank and n.
    """
    if isinstance(strip_prefix, str):
        strip_prefix = (strip_prefix,)
    rows, haps = find_haplotype_structures(annot, genes, region, **find_kws)
    anc_of = None
    if ancestry_df is not None:
        haps["ancestry"] = _ancestry_of(haps["sample"], ancestry_df, ancestry_id_col, ancestry_col)
        anc_of = haps["ancestry"]  # indexed by hap_id

    if Genes_dict is None:
        groups = {k: list(v.index) for k, v in haps.groupby("structure")}
    else:
        groups = {k: [str(m) for m in v] for k, v in Genes_dict.items()}
    counts = sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    total = sum(len(v) for _, v in counts)
    selected = _select_ranks(len(counts), hap_number, ranks)

    by_hap = rows.set_index("hap_id")
    available = set(haps.index)
    records, no_rep = [], []
    for rank in selected:
        key, members = counts[rank - 1]
        if rank in exclude_ranks:
            continue
        rep = next((m for m in members if m in available), None)
        if rep is None:
            no_rep.append(rank)
            continue
        n = len(members)
        for _, g in by_hap.loc[[rep]].iterrows():
            records.append(dict(
                rank=rank, n=n, pct=100 * n / total,
                gene=g["gene"], label=_strip(str(g[label_col]), strip_prefix),
                start=int(g["start"]), end=int(g["end"]), strand=g["strand"],
                hap_length=int(g["hap_length"]), sample=g["sample"], hap=g["hap"],
                contig=g["Contig"],
            ))
    if no_rep:
        print(f"No single-contig assembly found for structure ranks {no_rep}; skipped.")

    rank_of = {m: i for i, (_, v) in enumerate(counts, 1) for m in v}
    n_of = {m: len(v) for _, v in counts for m in v}
    haps = haps.reset_index()
    haps["rank"] = haps["hap_id"].map(rank_of).astype("Int64")
    haps["n"] = haps["hap_id"].map(n_of).astype("Int64")
    haps = haps.sort_values(["rank", "hap_id"]).reset_index(drop=True)

    feat = pd.DataFrame.from_records(records)
    feat.attrs["n_haplotypes"] = len(counts)
    feat.attrs["total"] = total
    if anc_of is not None and len(feat):
        # plain nested dict {rank: {ancestry: n}}: a DataFrame in .attrs breaks pd.concat
        feat.attrs["ancestry_counts"] = _ancestry_counts_by_rank(
            counts, anc_of, sorted(feat["rank"].unique())).to_dict(orient="index")
    return feat, haps


# --------------------------------------------------------------------------- #
# 2. Figure                                                                    #
# --------------------------------------------------------------------------- #
def _repel_1d(x, w, pad):
    """Spread label centres x (widths w) so they don't overlap; minimal movement."""
    order = np.argsort(x)
    xs = np.asarray(x, float)[order]
    ws = np.asarray(w, float)[order]

    def layout(idx):
        offs = [0.0]
        for a, b in zip(idx[:-1], idx[1:]):
            offs.append(offs[-1] + (ws[a] + ws[b]) / 2 + pad)
        offs = np.array(offs)
        return offs + np.mean(xs[idx] - offs)

    blocks = []
    for i in range(len(xs)):
        blocks.append([i])
        while len(blocks) > 1:
            prev, cur = blocks[-2], blocks[-1]
            if layout(cur)[0] - layout(prev)[-1] < (ws[cur[0]] + ws[prev[-1]]) / 2 + pad:
                blocks[-2] = prev + cur
                blocks.pop()
            else:
                break
    pos = np.empty_like(xs)
    for b in blocks:
        pos[b] = layout(b)
    out = np.empty_like(pos)
    out[order] = pos
    return out


def _gene_glyph(x0, x1, y, half_h, head, strand):
    """Pentagon arrow (box if strand unknown)."""
    head = min(head, 0.6 * (x1 - x0))
    if strand == "+":
        pts = [(x0, y - half_h), (x1 - head, y - half_h), (x1, y),
               (x1 - head, y + half_h), (x0, y + half_h)]
    elif strand == "-":
        pts = [(x1, y - half_h), (x0 + head, y - half_h), (x0, y),
               (x0 + head, y + half_h), (x1, y + half_h)]
    else:
        pts = [(x0, y - half_h), (x1, y - half_h), (x1, y + half_h), (x0, y + half_h)]
    return Polygon(pts, closed=True)


def _ancestry_color_map(labels, overrides=None):
    """{label: hex}; known labels use ANCESTRY_COLORS, the rest cycle ANCESTRY_FALLBACK."""
    fb = iter(ANCESTRY_FALLBACK)
    out = {}
    for lab in labels:
        out[lab] = ANCESTRY_COLORS.get(str(lab).lower()) or next(fb, PSEUDO_COLOR)
    if overrides:
        out.update({k: to_hex(v) for k, v in overrides.items()})
    return out


def _stacked_ancestry_bars(ax, ys, ranks, anc, color_map, height, mode, baseline=0.8):
    """Stacked horizontal bars on the (log) count axis. Each bar spans baseline ..
    baseline + n, like the plain bars.
    mode "count": segment edges sit at cumulative counts (read against the axis; small
                  groups stacked on large ones become thin on a log scale).
    mode "share": the bar's log-length is divided in proportion to each ancestry's share,
                  so composition is comparable across bars; only the bar end reads as n."""
    for y, r in zip(ys, ranks):
        row = anc.loc[r]
        total = row.sum()
        cum = 0.0
        for lab, c in row.items():
            if c <= 0:
                continue
            if mode == "count":
                x0, x1 = baseline + cum, baseline + cum + c
            else:
                span = (baseline + total) / baseline
                x0, x1 = baseline * span ** (cum / total), baseline * span ** ((cum + c) / total)
            ax.barh(y, x1 - x0, left=x0, height=height, color=color_map[lab],
                    edgecolor="white", linewidth=0.3, zorder=2)
            cum += c


def _style_axes(a):
    """Apply the publication tick / spine style explicitly.

    Axes made by the caller (ax=) were created outside plt.rc_context(PUB_RC), so they
    don't pick those rcParams up; setting them here keeps panels identical either way."""
    a.tick_params(axis="x", labelsize=6.5, color=INK_SECONDARY, labelcolor=INK,
                  width=0.6, length=2.5)
    a.tick_params(axis="y", labelsize=7, length=0, labelcolor=INK)
    for s in a.spines.values():
        s.set_linewidth(0.6)
        s.set_color(INK_SECONDARY)
    for lbl in (a.xaxis.label, a.yaxis.label):
        lbl.set_size(7)
        lbl.set_color(INK)


def make_panel_figure(row_counts, width_mm=180, row_mm=6.2, left_mm=15, right_mm=2,
                      bottom_mm=11, top_mm=12, panel_gap_mm=18, dpi=200):
    """Figure with one structure panel per entry of `row_counts` (how many haplotype
    rows each panel will draw), stacked top to bottom and sized in mm so that one row
    is `row_mm` tall in every panel regardless of how many rows it has.

    Returns (fig, axes); pass axes[i] as `ax=` to plot_haplotype_structures.
    `panel_gap_mm` must leave room for a panel's x-axis label plus the next panel's
    legend (~18 mm with the default two legend rows + ancestry row).

    Don't call tight_layout() / constrained_layout on the result: the panels are
    placed in absolute mm, which is what keeps the geometry reproducible.

        fig, axes = make_panel_figure([len(f["rank"].unique()) for f in feats])
        for a, f, t in zip(axes, feats, titles):
            plot_haplotype_structures(f, ax=a, title=t)
        fig.savefig("panels.pdf", facecolor="white", dpi=600)
    """
    counts = [int(n) for n in row_counts]
    panel_h = [n * row_mm for n in counts]
    height_mm = top_mm + bottom_mm + sum(panel_h) + panel_gap_mm * (len(counts) - 1)
    with plt.rc_context(PUB_RC):
        fig = plt.figure(figsize=(width_mm * MM, height_mm * MM), dpi=dpi)
        axes, y = [], height_mm - top_mm
        for h in panel_h:
            y -= h
            axes.append(fig.add_axes([left_mm / width_mm, y / height_mm,
                                      (width_mm - left_mm - right_mm) / width_mm,
                                      h / height_mm]))
            y -= panel_gap_mm
    return fig, axes


def plot_haplotype_structures(feat, out_prefix=None, *, gene_table=GENE_TABLE, palette="Set2",
                              color_by="auto", colors=None, ref_rank=1, align="first_shared",
                              label_genes=True, show_strand=True,
                              min_glyph_bp=1000, x_tick_kb=None,
                              width_mm=180, row_mm=6.2,
                              highlight_genes=None, gene_height=0.36, highlight_height=0.56,
                              highlight_min_mm=2.0, highlight_edge=True,
                              formats=("pdf", "png"), dpi=600, display_dpi=200,
                              save_source_data=True, ancestry_counts=None,
                              ancestry_colors=None, ancestry_stack="share",
                              ax=None, fig=None, show_freq=True, freq_xmax=None,
                              legend="auto", title=None, return_axes=False):
    """Draw haplotype gene structures with a frequency panel on the right.

    Colour:
    gene_table : path to (or DataFrame of) the gene class TSV.
    palette    : matplotlib qualitative colormap name ('Set2', 'Dark2', 'tab10', ...)
                 or a list of colours, used for functional genes.
    color_by   : "gene", "group" (table's group column), or "auto" (gene if the palette
                 has enough colours, else group).
    colors     : optional {gene or group: colour} overrides.

    align : "first_shared" -> each haplotype is shifted so its leftmost gene that also
            occurs in the reference haplotype (ref_rank) sits at that gene's reference
            position, so a missing first gene shows as a gap instead of shifting the row.
            "start" -> every row starts at 0.

    Glyph sizes (heights are fractions of one row; widths in mm on the page):
    highlight_genes  : genes drawn larger; default = functional genes in the table.
    gene_height      : glyph height for all other genes (pseudogenes).
    highlight_height : glyph height for highlighted genes.
    highlight_min_mm : highlighted genes narrower than this on the page are widened to it,
                       centred on their true midpoint (so position is preserved). 0 = to scale.
    highlight_edge   : thin white outline so enlarged genes stay distinct from neighbours.

    Ancestry (count bars stratified by ancestry):
    ancestry_counts : rank x ancestry count table. Normally not passed: it is taken from
                      feat.attrs["ancestry_counts"] (a {rank: {ancestry: n}} dict), which
                      collect_haplotype_features / plot_gene_structure fill when given
                      ancestry_df.
    ancestry_colors : optional {ancestry label: colour} overrides.
    ancestry_stack  : "share" (default) splits each bar's log-length by ancestry share, so
                      composition is comparable across structures; "count" stacks true
                      counts on the log axis (segment edges read against the axis, but
                      minor ancestries on top of big ones get very thin).

    Drawing into an existing figure (subplots):
    ax          : None (default) -> a figure of its own is made, sized width_mm x
                  (rows * row_mm) as before.
                  An Axes -> that axes becomes the structure panel. It is narrowed in
                  place to make room for the frequency bars and the n / % columns, which
                  are added as a second axes inside its original box, so the panel takes
                  exactly the space the caller gave it. Its height is used as is, i.e.
                  row_mm is ignored and rows share whatever height the axes has.
                  A (ax_structure, ax_frequency) pair -> both are used as given and
                  nothing is repositioned.
                  Repositioning is incompatible with tight_layout / constrained_layout;
                  the layout engine is switched off (with a warning) when one is active.
                  See make_panel_figure() for a ready-made stack of panels.
    fig         : figure to draw into when `ax` is None; the usual mm layout is applied
                  to it (its size then stands in for width_mm / the row-derived height).
    show_freq   : False draws only the structures (no bars, no n / % columns) and gives
                  the whole panel width to the gene track.
    freq_xmax   : upper limit of the count axis. Set it to the same value in every panel
                  so bar lengths are comparable across subplots (default: per panel).
    legend      : "figure" (one legend at the top of the figure), "axes" (legend above
                  this panel, for subplots), False (none), or "auto" (default: "axes"
                  when `ax` was given, else "figure").
    title       : optional panel title, placed above the legend.
    return_axes : True returns (fig, ax_structure, ax_frequency) instead of just fig
                  (ax_frequency is None when show_freq=False).

    Output:
    formats     : any of "pdf", "svg" (vector, sharp at any zoom), "png", "tiff" (raster).
    dpi         : resolution of raster files (600 = print standard; 1200 for line-art
                  requirements). Vector files are unaffected.
    display_dpi : resolution of the preview shown in the notebook.
    out_prefix  : with `ax` / `fig` given this saves the whole shared figure, so it is
                  usually left out and the composed figure saved once by the caller.
    """
    if isinstance(formats, str):  # formats="png" -> ("png",); also accepts "pdf,png"
        formats = formats.replace(",", " ").split()
    formats = tuple(f.lower().lstrip(".") for f in formats)
    if ancestry_stack not in ("share", "count"):
        raise ValueError("ancestry_stack must be 'share' or 'count'")
    if legend not in ("auto", "figure", "axes", True, False, None):
        raise ValueError("legend must be 'auto', 'figure', 'axes' or False")
    host_ax, freq_ax = None, None
    if ax is not None:
        if isinstance(ax, (tuple, list)):
            host_ax = ax[0]
            freq_ax = ax[1] if len(ax) > 1 else None
            if freq_ax is not None:
                show_freq = True
        else:
            host_ax = ax
        fig = host_ax.figure
    external = host_ax is not None
    if legend in ("auto", True):
        legend = "axes" if external else "figure"
    anc = ancestry_counts if ancestry_counts is not None else feat.attrs.get("ancestry_counts")
    present = sorted(feat["rank"].unique())
    if anc is not None:
        if not isinstance(anc, pd.DataFrame):  # dict {rank: {ancestry: n}} from feat.attrs
            anc = pd.DataFrame.from_dict(anc, orient="index")
        anc = anc.reindex(present).fillna(0)
        anc = anc.loc[:, anc.sum() > 0]
        anc = anc[_order_ancestry(list(anc.columns))]
        anc_colors = _ancestry_color_map(list(anc.columns), ancestry_colors)
    if ref_rank not in present:  # e.g. only rare structures plotted: align to the commonest shown
        ref_rank = present[0]
    feat, legend_entries = classify_and_color(feat, gene_table, palette, color_by, colors, ref_rank)
    if highlight_genes is None:
        highlight = set(feat.loc[feat["class"] == "functional", "gene"])
    else:
        highlight = set(highlight_genes)
    ranks = list(dict.fromkeys(feat.sort_values("rank")["rank"]))
    n_rows = len(ranks)

    # ---- alignment -------------------------------------------------------
    ref = feat[feat["rank"] == ref_rank].sort_values("start")
    ref_pos = ref.drop_duplicates("gene").set_index("gene")["start"]
    feat["offset"] = 0
    if align == "first_shared":
        for r in ranks:
            sub = feat[feat["rank"] == r].sort_values("start")
            shared = sub[sub["gene"].isin(ref_pos.index)]
            if len(shared):
                feat.loc[sub.index, "offset"] = ref_pos[shared["gene"].iloc[0]] - shared["start"].iloc[0]
    feat["x0"] = (feat["start"] + feat["offset"]) / 1e3
    feat["x1"] = (np.maximum(feat["end"], feat["start"] + min_glyph_bp) + feat["offset"]) / 1e3
    anchor = ref["gene"].iloc[0]

    # ---- layout (explicit mm margins => identical geometry every time) -----
    leg_cols = min(len(legend_entries), 8)
    leg_rows = int(np.ceil(len(legend_entries) / max(leg_cols, 1)))
    left_mm, right_mm, bottom_mm = 15, 2, 11
    top_mm = 9 + 3.2 * (leg_rows - 1) + (4.5 if anc is not None else 0) + (4 if title else 0)
    freq_mm, cols_mm, gap_mm = (18, 21, 3) if show_freq else (0, 0, 0)

    if external or fig is not None:  # caller owns the figure: take its size in mm
        W, H = np.asarray(fig.get_size_inches()) / MM
    if external:
        pos = host_ax.get_position()
        panel_mm, panel_h_mm = pos.width * W, pos.height * H
        if freq_ax is not None:  # explicit pair: each panel keeps its own box
            struct_mm = panel_mm
            freq_mm = freq_ax.get_position().width * W
        else:
            struct_mm = panel_mm - freq_mm - cols_mm - gap_mm
        if struct_mm <= 0:
            raise ValueError(f"the axes is only {panel_mm:.0f} mm wide, too narrow for the "
                             f"frequency panel ({freq_mm + cols_mm + gap_mm:.0f} mm); widen it "
                             "or pass show_freq=False")
    else:
        if fig is not None:
            width_mm, height_mm = W, H
        else:
            height_mm = top_mm + bottom_mm + n_rows * row_mm
            W, H = width_mm, height_mm
        panel_h_mm = n_rows * row_mm
        struct_mm = width_mm - left_mm - right_mm - freq_mm - cols_mm - gap_mm

    feat["highlight"] = feat["gene"].isin(highlight)
    if highlight_min_mm:
        min_kb = highlight_min_mm * (feat["x1"].max() - feat["x0"].min()) * 1.02 / struct_mm
        mid = (feat["x0"] + feat["x1"]) / 2
        wide = np.maximum(feat["x1"] - feat["x0"], min_kb)
        m = feat["highlight"]
        feat.loc[m, "x0"] = mid[m] - wide[m] / 2
        feat.loc[m, "x1"] = mid[m] + wide[m] / 2

    xmin = feat["x0"].min()
    xmax = feat["x1"].max()
    xpad = 0.01 * (xmax - xmin)
    if x_tick_kb is None:  # ~6-9 ticks across the region
        span = xmax - xmin
        x_tick_kb = next(s for s in (5, 10, 20, 25, 50, 100, 200, 250, 500, 1e9) if span / s <= 9)

    with plt.rc_context(PUB_RC):
        if external:
            ax = host_ax
            axf = freq_ax
            if freq_ax is None and show_freq:
                if fig.get_layout_engine() is not None:
                    print("Note: the figure's automatic layout (tight / constrained) was "
                          "switched off; the panels are placed in mm.")
                    fig.set_layout_engine("none")
                if not getattr(ax, "_hapstruct_split", False):  # only narrow it once
                    ax.set_position([pos.x0, pos.y0, struct_mm / W, pos.height])
                    ax._hapstruct_split = True
                axf = fig.add_axes([pos.x0 + (struct_mm + gap_mm) / W, pos.y0,
                                    freq_mm / W, pos.height], sharey=ax)
        else:
            if fig is None:
                fig = plt.figure(figsize=(width_mm * MM, height_mm * MM), dpi=display_dpi)
            b, h = bottom_mm / H, panel_h_mm / H
            ax = fig.add_axes([left_mm / W, b, struct_mm / W, h])
            axf = (fig.add_axes([(left_mm + struct_mm + gap_mm) / W, b, freq_mm / W, h],
                                sharey=ax) if show_freq else None)
        panels = [a for a in (ax, axf) if a is not None]
        for a in panels:
            _style_axes(a)
        if freq_ax is not None:  # a caller-made pair isn't linked yet: rows must line up
            try:
                freq_ax.sharey(ax)
            except ValueError:
                freq_ax.set_ylim(ax.get_ylim())

        ax.set_xlim(xmin - xpad, xmax + xpad)
        ax.set_ylim(n_rows - 0.5, -0.5)  # row 0 on top
        kb_per_mm = (xmax - xmin + 2 * xpad) / struct_mm
        kb_per_pt = kb_per_mm * 25.4 / 72

        # Within a row: label band on top, gene track below.
        half_h, half_hl = gene_height / 2, highlight_height / 2
        tallest = max(half_h, half_hl if highlight else 0)
        y_off_track = min(0.5 - tallest - 0.01, 0.2)  # sit the track as low as fits
        y_off_label = y_off_track - tallest - 0.03
        head = 0.5 * kb_per_mm  # 0.5 mm arrow head
        label_fs = 6
        char_w_pt = 0.62 * label_fs

        for i, r in enumerate(ranks):
            sub = feat[feat["rank"] == r].sort_values("x0")
            yt = i + y_off_track
            if i % 2 == 1:
                for a in panels:
                    a.axhspan(i - 0.5, i + 0.5, color=ZEBRA, lw=0, zorder=0)

            ax.plot([sub["x0"].min(), sub["x1"].max()], [yt, yt], color=BACKBONE,
                    lw=0.8, solid_capstyle="butt", zorder=2)
            for _, g in sub.iterrows():
                hl = g["highlight"]
                p = _gene_glyph(g["x0"], g["x1"], yt, half_hl if hl else half_h,
                                head * (1.4 if hl else 1), g["strand"] if show_strand else None)
                p.set_facecolor(g["color"])
                p.set_edgecolor("white")
                p.set_linewidth(0.4 if (hl and highlight_edge) else 0)
                p.set_zorder(4 if hl else 3)  # enlarged genes drawn on top of neighbours
                ax.add_patch(p)

            if label_genes:
                labels = sub["label"].tolist()
                centers = ((sub["x0"] + sub["x1"]) / 2).to_numpy()
                widths = np.array([len(s) * char_w_pt * kb_per_pt for s in labels])
                xs = _repel_1d(centers, widths, pad=1.0 * kb_per_pt)
                for x, s in zip(xs, labels):
                    ax.text(x, i + y_off_label, s, ha="center", va="bottom",
                            fontsize=label_fs, fontstyle="italic", color=INK, zorder=4)

        # ---- structure axis --------------------------------------------
        ax.set_yticks([i + y_off_track for i in range(n_rows)])
        ax.set_yticklabels([f"Hap {r}" for r in ranks])
        ax.tick_params(axis="y", pad=3)
        ax.xaxis.set_major_locator(MultipleLocator(x_tick_kb))
        ax.grid(axis="x", color=GRID, lw=0.4, zorder=1)
        ax.set_axisbelow(True)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        gene_tex = "-".join("$\\it{" + p + "}$" for p in str(anchor).split("-") if p)
        xlab = f"Position relative to {gene_tex} (kb)" if align == "first_shared" else "Position (kb)"
        ax.set_xlabel(xlab, labelpad=2)

        # ---- frequency panel ---------------------------------------------
        first = feat.drop_duplicates("rank").set_index("rank").loc[ranks]
        ys = [i + y_off_track for i in range(n_rows)]
        if axf is not None:
            if anc is None:
                axf.barh(ys, first["n"], height=2 * half_h, color=BAR, lw=0, zorder=2, left=0.8)
            else:
                if not np.allclose(anc.sum(axis=1).to_numpy(), first["n"].to_numpy()):
                    print("Warning: ancestry counts don't add up to n for every structure.")
                _stacked_ancestry_bars(axf, ys, ranks, anc, anc_colors, 2 * half_h, ancestry_stack)
            axf.set_xscale("log")
            top = freq_xmax if freq_xmax else 10 ** max(1, np.ceil(np.log10(first["n"].max())))
            axf.set_xlim(0.8, top)
            decades = [10 ** k for k in range(int(round(np.log10(axf.get_xlim()[1]))) + 1)]
            major = decades if len(decades) <= 4 else decades[(len(decades) - 1) % 2::2]  # keep the end labelled
            axf.xaxis.set_major_locator(FixedLocator(major))
            axf.set_xticklabels([f"{d:,}" if d < 1000 else f"{d // 1000:,}k" for d in major])
            axf.xaxis.set_minor_locator(FixedLocator([d for d in decades if d not in major]))
            axf.xaxis.set_minor_formatter(NullFormatter())
            axf.tick_params(axis="x", which="minor", length=1.5, width=0.5)
            axf.grid(axis="x", color=GRID, lw=0.4, zorder=1)
            axf.set_axisbelow(True)
            axf.tick_params(axis="y", left=False, labelleft=False)
            for s in ("top", "right", "left"):
                axf.spines[s].set_visible(False)
            axf.set_xlabel("Count (log)", labelpad=2)

            # n / % columns, right of the bars
            tr = axf.get_yaxis_transform()
            x_n, x_p = 1 + 10 / freq_mm, 1 + (cols_mm - 0.5) / freq_mm
            for y, (n, pct) in zip(ys, first[["n", "pct"]].to_numpy()):
                axf.text(x_n, y, f"{int(n):,}", transform=tr, ha="right", va="center",
                         fontsize=6.5, color=INK)
                axf.text(x_p, y, f"{pct:.2f}%", transform=tr, ha="right", va="center",
                         fontsize=6.5, color=INK_SECONDARY)
            for x, s in ((x_n, "n"), (x_p, "%")):
                axf.text(x, -0.62, s, transform=tr, ha="right", va="bottom",
                         fontsize=6.5, color=INK_SECONDARY, fontweight="bold")

        # ---- legend ------------------------------------------------------
        # One figure-wide legend, or (for subplots) a legend stacked above this panel:
        # gene rows on top, the ancestry row just above the axes, title above both.
        handles = [Patch(fc=c, ec="none", label=lab) for lab, c, _ in legend_entries]
        leg_kw = dict(frameon=False, fontsize=6.5, handlelength=1.4, handleheight=0.8,
                      handletextpad=0.5, columnspacing=1.4, borderaxespad=0, labelspacing=0.4)
        y_top = 1.0  # top of what has been drawn so far, in axes fractions of `ax`
        if legend == "figure":
            leg = fig.legend(handles=handles, loc="upper left", ncol=leg_cols,
                             bbox_to_anchor=(left_mm / W, 1 - (1.2 + (4 if title else 0)) / H),
                             **leg_kw)
        elif legend == "axes":
            y_anc = 1 + 1.5 / panel_h_mm
            if anc is not None:
                ax.text(0, y_anc - 0.3 / panel_h_mm, "Ancestry", transform=ax.transAxes,
                        ha="left", va="bottom", fontsize=6.5, fontweight="bold", color=INK)
                anc_leg = ax.legend(handles=[Patch(fc=anc_colors[c], ec="none", label=str(c))
                                             for c in anc.columns],
                                    loc="lower left", ncol=len(anc.columns),
                                    bbox_to_anchor=(12 / struct_mm, y_anc), **leg_kw)
                ax.add_artist(anc_leg)  # keep it: the gene legend below replaces ax.legend_
            y_gene = y_anc + ((3.6 if anc is not None else 0)) / panel_h_mm
            leg = ax.legend(handles=handles, loc="lower left", ncol=leg_cols,
                            bbox_to_anchor=(0, y_gene), **leg_kw)
            y_top = y_gene + (3.2 * leg_rows + 1.0) / panel_h_mm
        else:
            leg = None
        if leg is not None:
            for t, (_, _, italic) in zip(leg.get_texts(), legend_entries):
                if italic:
                    t.set_fontstyle("italic")

        if legend == "figure" and anc is not None:  # second legend line, under the gene legend
            y_leg = 1 - (1.2 + 3.2 * leg_rows + 0.8 + (4 if title else 0)) / H
            fig.text(left_mm / W, y_leg - 1.1 / H, "Ancestry", ha="left", va="top",
                     fontsize=6.5, fontweight="bold", color=INK)
            fig.legend(handles=[Patch(fc=anc_colors[c], ec="none", label=str(c)) for c in anc.columns],
                       loc="upper left", ncol=len(anc.columns),
                       bbox_to_anchor=((left_mm + 12) / W, y_leg), **leg_kw)

        if title:
            if legend == "figure":
                fig.text(left_mm / W, 1 - 1.2 / H, title, ha="left", va="top",
                         fontsize=8, fontweight="bold", color=INK)
            else:
                ax.text(0, y_top, title, transform=ax.transAxes, ha="left", va="bottom",
                        fontsize=8, fontweight="bold", color=INK)

        if out_prefix:
            os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
            for ext in formats:
                kw = dict(facecolor="white", dpi=dpi)
                if ext in ("tif", "tiff"):
                    kw["pil_kwargs"] = {"compression": "tiff_lzw"}
                fig.savefig(f"{out_prefix}.{ext}", **kw)
            if save_source_data:
                cols = ["rank", "n", "pct", "gene", "label", "class", "group", "color",
                        "start", "end", "strand", "hap_length", "offset"]
                cols += [c for c in ("sample", "hap", "contig") if c in feat]
                feat[cols].to_csv(f"{out_prefix}_source_data.tsv", sep="\t", index=False)
                if anc is not None:
                    anc.to_csv(f"{out_prefix}_ancestry_counts.tsv", sep="\t")
            print(f"Saved {out_prefix}.{{{','.join(formats)}}}")
    return (fig, ax, axf) if return_axes else fig


# --------------------------------------------------------------------------- #
# 3. One-call wrapper                                                          #
# --------------------------------------------------------------------------- #
def plot_gene_structure(annot, filename, hap_number=14, genes=None, region=None, t="Gene", ranks=None,
                        r=0, fontsize=None, exclude_ranks=(), Genes_dict=None,
                        sample_col="Sample", hap_col=None, bounds=None, crop_to_bounds=False,
                        orient_by=None, fragmented="drop", out_dir="/home/jupyter/Gene_Structure",
                        ancestry_df=None, ancestry_id_col="research_id",
                        ancestry_col="ancestry_pred", **plot_kws):
    """Annotation table -> figure, in one call.

    annot    : the single annotation DataFrame (Contig, Category, Spos, Epos, Strand,
               Gene, Sample, ...).
    filename : only names the outputs: <out_dir>/<basename without .json>.{pdf,png} plus
               _source_data.tsv and _haplotype_assignments.tsv (every haplotype's structure).
    genes / region : which genes define the structure (see find_haplotype_structures).
    hap_number     : top N structures; -N = the N rarest. ranks= picks explicit ranks.
    bounds         : e.g. ("HLA-F", "HLA-A") keeps only complete annotations (both genes on
                     one contig); crop_to_bounds=True also drops genes outside them.
    Genes_dict     : optional precomputed groups; otherwise groups come from `annot`.
    ancestry_df    : optional DataFrame with columns "research_id" (matched to `sample_col`)
                     and "ancestry_pred"; stratifies the count bars by ancestry. Also
                     accepts ancestry_stack="share"|"count" and ancestry_colors={...}.
    `r` and `fontsize` are accepted but ignored. Extra keywords (gene_table, palette,
    color_by, highlight_height, ...) go to plot_haplotype_structures. For subplots use
    collect_haplotype_features + plot_haplotype_structures(ax=...) directly: this wrapper
    always writes one file set per call.

    Returns (number of distinct structures, minor_haplotype_info) as before:
    minor_haplotype_info = [[(structure, 1), [hap_id]], ...] for singletons in the top N.
    """
    feat, haps = collect_haplotype_features(
        annot, genes, region, hap_number=hap_number, ranks=ranks, exclude_ranks=exclude_ranks,
        Genes_dict=Genes_dict, label_col=t, sample_col=sample_col, hap_col=hap_col,
        bounds=bounds, crop_to_bounds=crop_to_bounds, orient_by=orient_by,
        fragmented=fragmented, ancestry_df=ancestry_df,
        ancestry_id_col=ancestry_id_col, ancestry_col=ancestry_col)
    prefix = os.path.join(out_dir, os.path.basename(filename).replace(".json", ""))
    plot_haplotype_structures(feat, prefix, **plot_kws)
    haps.to_csv(f"{prefix}_haplotype_assignments.tsv", sep="\t", index=False)
    print("All Haplotype Number: %s" % feat.attrs["n_haplotypes"])

    if Genes_dict is not None:
        groups = sorted(((k, list(v)) for k, v in Genes_dict.items()),
                        key=lambda kv: (-len(kv[1]), kv[0]))
    else:
        groups = sorted(((k, list(v)) for k, v in haps.groupby("structure")["hap_id"]),
                        key=lambda kv: (-len(kv[1]), kv[0]))
    sel = [groups[r - 1] for r in _select_ranks(len(groups), hap_number, ranks)]
    minor = [[(k, 1), m] for k, m in sel if len(m) == 1]
    return feat.attrs["n_haplotypes"], minor
