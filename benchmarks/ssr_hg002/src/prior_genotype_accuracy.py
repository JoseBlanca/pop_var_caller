#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# ///
"""PRIOR COMPARISON (derived from genotype_accuracy.py). Genotype-LENGTH accuracy (not just detection): among truth-positive loci a caller
detects, does it get HG002's EXACT allele lengths right? This is where an STR-aware
caller should earn its keep — freebayes can flag "a length changed" without nailing
the repeat count.

Common currency = the genotype as a sorted pair of per-allele bp deltas relative to
the reference genome:
  truth : GIAB is phased; per haplotype, sum (len(ALT_on_that_hap) - len(REF)) over
          the overlapping records -> (hap0_bp, hap1_bp).
  ours  : sorted(len(ALT_i) - len(REF_tract)) over GT  (REF tract == reference).
  HipSTR: the GB field IS bp-diff-from-reference per allele.
  freebayes: sorted bp deltas from the best-overlapping indel record's GT (APPROX —
          freebayes is usually unphased and can split a tract across records).

Reported per coverage, over loci that are truth length-variant (period-aware) AND
have a well-defined phased truth genotype:
  detected      = caller called the locus a length variant (the detection recall)
  gt_correct    = detected AND exact genotype match
  acc | detected = gt_correct / detected  (given you detected it, did you nail it?)

Usage: genotype_accuracy.py   (from ssr_hg002/)
"""
import bisect, gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CATALOG = ROOT / "catalog/HG002_Tier_restricted.cat"
TRUTH = ROOT / "truth/HG002_GRCh38_TandemRepeats_v1.0.1_50000.vcf.gz"
ARMS = {"plug-in": lambda c: ROOT / f"results/ours_prior_comparison/plugin/HG002_{c}x.ssr.vcf",
        "marginal": lambda c: ROOT / f"results/ours_prior_comparison/marginalized/HG002_{c}x.ssr.vcf"}
HIP = lambda c: ROOT / f"results/hipstr/HG002_{c}x.str.vcf.gz"
FB = lambda c: ROOT / f"results/freebayes/HG002_{c}x.fb.vcf.gz"
FB_MINQUAL = 20.0
COVERAGES = [50, 30, 20, 15]
MAX_BACK = 300


def _open(p):
    p = str(p)
    return gzip.open(p, "rt") if p.endswith((".gz", ".bgz", ".cat")) else open(p)


def phased(gt):
    if "|" in gt:
        parts = gt.split("|")
    else:
        parts = gt.split("/")
    if any(x in (".", "") for x in parts):
        return None
    return [int(x) for x in parts], ("|" in gt)


# ---- catalog loci -----------------------------------------------------------
loci = {}
with _open(CATALOG) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        loci[(f[0], int(f[1]) + 1)] = {"start": int(f[1]), "end": int(f[2]), "period": len(f[3])}


# ---- generic overlap index of (s0, e0, record) ------------------------------
def build_index(path, min_qual=None):
    iv = {}
    if not Path(path).exists():
        return iv
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            if min_qual is not None:
                try:
                    if float(f[5]) < min_qual:
                        continue
                except ValueError:
                    pass
            ref, alt = f[3], f[4]
            alleles = [ref] + ([] if alt == "." else alt.split(","))
            s0 = int(f[1]) - 1
            span = max([len(ref)] + [len(a) for a in alleles[1:]])
            iv.setdefault(f[0], []).append((s0, s0 + span, ref, alleles, f[8], f[9]))
    return {c: (sorted(v), [x[0] for x in sorted(v)]) for c, v in iv.items()}


def overlapping(index, chrom, s, e):
    entry = index.get(chrom)
    if not entry:
        return []
    v, st = entry
    j = bisect.bisect_right(st, e) - 1
    out = []
    while j >= 0 and st[j] >= s - MAX_BACK:
        if v[j][0] < e and v[j][1] > s:
            out.append(v[j])
        j -= 1
    return out


def gfield(fmt, sample, key):
    parts = fmt.split(":")
    if key not in parts:
        return None
    return sample.split(":")[parts.index(key)]


# ---- truth phased genotype per locus (bp deltas) ----------------------------
truth_idx = build_index(TRUTH)


def truth_gt(chrom, s, e, period):
    """sorted (hap0_bp, hap1_bp) or None; also whether it's a period-aware variant."""
    hap = [0, 0]
    seen = False
    ok = True
    for s0, e0, ref, alleles, fmt, sample in overlapping(truth_idx, chrom, s, e):
        g = phased(gfield(fmt, sample, "GT"))
        if g is None:
            continue
        idx, is_phased = g
        if not is_phased or len(idx) != 2:
            ok = False
            continue
        for h in (0, 1):
            hap[h] += len(alleles[idx[h]]) - len(ref)
        seen = True
    if not seen or not ok:
        return None
    gt = tuple(sorted(hap))
    is_var = any(d != 0 and (not period or d % period == 0) for d in gt)
    return gt if is_var else None


# ---- caller genotype per locus (bp deltas) ----------------------------------
def ours_hip_gt(path, is_hipstr):
    """(chrom,pos) -> sorted tuple of bp deltas, or None if no-call."""
    out = {}
    if not Path(path).exists():
        return out
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or f[6] not in (".", "PASS", ""):
                continue
            g = phased(gfield(f[8], f[9], "GT"))
            if g is None:
                continue
            idx = g[0]
            if is_hipstr:
                gb = gfield(f[8], f[9], "GB")
                if gb and all(x.lstrip("-").isdigit() for x in gb.replace("|", "/").split("/")):
                    out[(f[0], int(f[1]))] = tuple(sorted(int(x) for x in gb.replace("|", "/").split("/")))
                    continue
            alleles = [f[3]] + ([] if f[4] == "." else f[4].split(","))
            out[(f[0], int(f[1]))] = tuple(sorted(len(alleles[i]) - len(f[3]) for i in idx))
    return out


def fb_gt(index, chrom, s, e):
    """APPROX freebayes genotype: bp deltas from the overlapping indel record with the
    largest length change that is in the GT."""
    best = None
    for s0, e0, ref, alleles, fmt, sample in overlapping(index, chrom, s, e):
        g = phased(gfield(fmt, sample, "GT"))
        if g is None:
            continue
        idx = g[0]
        deltas = tuple(sorted(len(alleles[i]) - len(ref) for i in idx))
        if all(d == 0 for d in deltas):
            continue
        mag = max(abs(d) for d in deltas)
        if best is None or mag > best[0]:
            best = (mag, deltas)
    return best[1] if best else None


# ---- compare ----------------------------------------------------------------
def is_het(gt):
    return gt[0] != gt[1]


print(f"{'cov':>4} {'prior':>9} {'truthpos':>8} {'detected':>8} {'correct':>8} "
      f"{'recall':>7} {'acc|det':>8} | {'homvar>het':>10} {'het>homvar':>10} {'other':>6}")
for cov in COVERAGES:
    calls = {name: ours_hip_gt(path(cov), False) for name, path in ARMS.items()}
    tot = 0
    stat = {n: [0, 0, 0, 0, 0] for n in ARMS}  # det, correct, homvar->het, het->homvar, other
    both_det = 0
    disagree = 0
    for (chrom, pos), m in loci.items():
        tg = truth_gt(chrom, m["start"], m["end"], m["period"])
        if tg is None:
            continue
        tot += 1
        got = {}
        for name in ARMS:
            g = calls[name].get((chrom, pos))
            got[name] = g
            if g is None or all(d == 0 for d in g):
                continue
            s = stat[name]
            s[0] += 1
            if g == tg:
                s[1] += 1
            elif not is_het(tg) and is_het(g):
                s[2] += 1
            elif is_het(tg) and not is_het(g):
                s[3] += 1
            else:
                s[4] += 1
        a, b = (got[n] for n in ARMS)
        if a is not None and b is not None:
            both_det += 1
            if a != b:
                disagree += 1
    for name in ARMS:
        d, c, hh, hg, o = stat[name]
        print(f"{cov:>4} {name:>9} {tot:>8} {d:>8} {c:>8} "
              f"{d/tot if tot else 0:>7.3f} {c/d if d else 0:>8.3f} | {hh:>10} {hg:>10} {o:>6}")
    print(f"{'':>4} {'(both detected ' + str(both_det) + ', genotypes differ ' + str(disagree) + ')':>60}")
