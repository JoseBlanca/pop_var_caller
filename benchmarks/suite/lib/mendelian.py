"""Mendelian consistency of a trio's joint calls (plan §4.2).

A site is checked when all three members have a called genotype and at least one carries a
non-reference allele. It is an error when the child's two alleles cannot be one from the mother and
one from the father. The count is split into SNPs and indels by the record's alleles.
"""

import subprocess

from .results import row


def _alleles(gt):
    gt = gt.replace("|", "/")
    parts = gt.split("/")
    if len(parts) != 2 or "." in parts:
        return None
    return parts


def _variant_class(ref, alts):
    return "SNP" if all(len(a) == 1 == len(ref) for a in alts.split(",") if a != "*") else "indel"


def count(vcf, bed, mother, father, child, min_qual=None):
    """Return {class: (sites checked, errors)} for the trio's calls inside `bed`."""
    qual = f"-i 'QUAL>={min_qual}'" if min_qual is not None else ""
    cmd = (f"bcftools view --no-version -f PASS,. -T {bed} -s {mother},{father},{child} {qual} {vcf} "
           f"| bcftools query -f '%REF\\t%ALT[\\t%GT]\\n'")
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    counts = {"SNP": [0, 0], "indel": [0, 0]}
    for line in out.splitlines():
        ref, alts, gm, gf, gc = line.split("\t")
        m, f, c = _alleles(gm), _alleles(gf), _alleles(gc)
        if m is None or f is None or c is None:
            continue
        if all(a == "0" for a in m + f + c):
            continue
        cls = _variant_class(ref, alts)
        counts[cls][0] += 1
        ok = (c[0] in m and c[1] in f) or (c[1] in m and c[0] in f)
        if not ok:
            counts[cls][1] += 1
    return counts


def rows(base, counts):
    out = []
    for cls, (checked, errors) in counts.items():
        out.append(row(base, variant_class=cls, stage="Mendelian", metric="sites_checked",
                       value=checked))
        out.append(row(base, variant_class=cls, stage="Mendelian", metric="mendelian_errors",
                       value=errors))
        out.append(row(base, variant_class=cls, stage="Mendelian", metric="mendelian_error_rate",
                       value=round(errors / checked, 6) if checked else ""))
    return out
