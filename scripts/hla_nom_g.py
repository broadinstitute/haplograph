#!/usr/bin/env python3
"""Map IPD-IMGT/HLA G groups to 2-field (4-digit) alleles.

The WMDA file ``hla_nom_g.txt`` is semicolon-separated::

    A*;01:01:01:01/01:01:01:02N/01:32/01:37:01:01;01:01:01G
    A*;01:01:02;

Column 1 is the locus, column 2 the slash-separated alleles that share the
same antigen-recognition-domain exons, and column 3 the G code. A blank G
code means the allele is not assigned to a G group and is skipped.

Each member is rolled up to two fields, so ``01:01:01:02N`` and ``01:32``
become ``01:01`` and ``01:32``. A G group can contain several 4-digit alleles,
and the same 4-digit allele can sit in more than one G group when those
groups differ by a synonymous change in exon 2 or 3.

Notebook use
------------
    import sys; sys.path.append("scripts")
    from hla_nom_g import g_to_four_digit

    g_groups = g_to_four_digit()
    g_groups["HLA-A*01:01:01G"]   # ["HLA-A*01:01", "HLA-A*01:04", ...]
"""

from __future__ import annotations

import re
import urllib.request
from pathlib import Path

__all__ = ["HLA_NOM_G_URL", "g_to_four_digit"]

HLA_NOM_G_URL = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_g.txt"

# Expression status is written on the last field: N, L, S, C, A, Q.
_EXPR_RE = re.compile(r"[NLSCAQ]$", re.IGNORECASE)


def _four_digit(allele: str) -> str | None:
    """``01:01:01:02N`` -> ``01:01``. None when the name has fewer than two fields."""
    fields = [part for part in allele.split(":") if part]
    if not fields:
        return None
    fields[-1] = _EXPR_RE.sub("", fields[-1])
    fields = [part for part in fields if part]
    if len(fields) < 2:
        return None
    return ":".join(fields[:2])


def _read_nomenclature(source: str | Path | None) -> str:
    if source is None:
        source = HLA_NOM_G_URL
    text = str(source)
    if text.startswith(("http://", "https://")):
        request = urllib.request.Request(text, headers={"User-Agent": "haplograph-hla-nom-g"})
        with urllib.request.urlopen(request, timeout=60) as response:
            return response.read().decode("utf-8")
    path = Path(source).expanduser()
    if not path.is_file():
        raise FileNotFoundError(f"hla_nom_g file not found: {path}")
    return path.read_text(encoding="utf-8")


def g_to_four_digit(
    source: str | Path | None = None,
    *,
    prefix: str = "HLA-",
) -> dict[str, list[str]]:
    """Map each G allele to the distinct 4-digit alleles in that group.

    Keys are G alleles (``HLA-A*01:01:01G``) and values are 2-field alleles in
    the order they first appear in the file. ``source`` is a path or URL; the
    default is the current IPD-IMGT/HLA release on GitHub. ``prefix`` is
    prepended to the WMDA locus (``A`` becomes ``HLA-A``). Pass ``prefix=""``
    to keep the file's ``A*01:01:01G`` names.

    Returns an empty-free dictionary. Raises ``ValueError`` when the file
    contains no G groups.
    """
    mapping: dict[str, list[str]] = {}
    seen: dict[str, set[str]] = {}
    for raw in _read_nomenclature(source).splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        locus, _, rest = line.partition(";")
        alleles, _, g_code = rest.partition(";")
        g_code = g_code.strip()
        if not g_code.endswith("G"):
            continue
        locus_name = locus.strip()
        if locus_name.endswith("*"):
            locus_name = locus_name[:-1]
        if not locus_name:
            continue
        gene = locus_name if not prefix or locus_name.startswith(prefix) else f"{prefix}{locus_name}"
        g_allele = f"{gene}*{g_code}"
        bucket = mapping.setdefault(g_allele, [])
        known = seen.setdefault(g_allele, set())
        for allele in alleles.split("/"):
            four = _four_digit(allele.strip())
            if four is None or four in known:
                continue
            known.add(four)
            bucket.append(f"{gene}*{four}")
    if not mapping:
        raise ValueError(f"no G groups in {source or HLA_NOM_G_URL}")
    return mapping
