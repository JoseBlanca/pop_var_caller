"""Genotype against genotype, sample by sample, for a simulation whose every genotype is known.

vcfeval scores one sample per run and needs about a minute a run; a cohort of 50 over five depths and
five callers would be 1,250 runs. A simulation that plants only biallelic SNPs needs no haplotype
matching, so each caller's genotype is looked up at every planted SNP for every sample instead:

- the calls are kept as vcfeval keeps them (PASS or no filter; freebayes at QUAL >= 20), split into
  one record per ALT and into single-base records (`bcftools norm -m -any -a`), so a caller that
  writes two nearby SNPs as one record is read at each;
- at a planted SNP, a sample's called genotype is the number of copies of the planted ALT it
  carries; a sample with no record there is homozygous reference, one with a record whose genotype
  is missing is `no_call`, and one carrying only another ALT is `other_allele`;
- a sample carrying an ALT at a position where nothing was planted is a false variant genotype.

The counts are kept per group of samples (`group_of(sample)`), and over all samples.
"""

import subprocess
from collections import defaultdict

from .results import row

TRUE = ["hom_ref", "het", "hom_alt"]
CALLED = ["hom_ref", "het", "hom_alt", "no_call", "other_allele"]


def _alt_copies(gt):
    """Copies of allele 1 in a genotype, or None when any allele is missing."""
    alleles = gt.replace("|", "/").split("/")
    if any(a == "." for a in alleles):
        return None
    return sum(a == "1" for a in alleles)


def _query(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True,
                          executable="/bin/bash").stdout.splitlines()


def truth_genotypes(truth, samples):
    """{(pos, ref, alt): [alt copies per sample]} from a biallelic truth VCF."""
    out = {}
    for line in _query(f"bcftools query -s {','.join(samples)} -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n' "
                       f"{truth}"):
        pos, ref, alt, *gts = line.split("\t")
        out[(int(pos), ref, alt)] = [_alt_copies(g) for g in gts]
    return out


def call_genotypes(calls, reference, samples, min_qual=None):
    """{pos: {alt: [alt copies or None per sample]}}, SNP records only, after the split."""
    qual = f"-i 'QUAL>={min_qual}'" if min_qual is not None else ""
    cmd = (f"bcftools view --no-version -f PASS,. {qual} -s {','.join(samples)} {calls} "
           f"| bcftools norm --no-version -f {reference} -m -any -a --atom-overlaps . "
           f"| bcftools view --no-version -v snps "
           f"| bcftools query -f '%POS\\t%ALT[\\t%GT]\\n'")
    out = defaultdict(dict)
    for line in _query(cmd):
        pos, alt, *gts = line.split("\t")
        copies = [_alt_copies(g) for g in gts]
        seen = out[int(pos)].get(alt)
        if seen is not None:  # two records gave this ALT here: keep the larger count
            copies = [max(a or 0, b or 0) if (a is not None or b is not None) else None
                      for a, b in zip(seen, copies)]
        out[int(pos)][alt] = copies
    return out


def compare(truth, calls, samples):
    """Per sample: {(true, called): count} and the false variant genotypes."""
    table = {s: defaultdict(int) for s in samples}
    false_variants = {s: 0 for s in samples}
    planted = {pos for pos, _, _ in truth}
    for (pos, _, alt), true_copies in truth.items():
        here = calls.get(pos, {})
        mine = here.get(alt)
        for i, s in enumerate(samples):
            t = true_copies[i]
            if t is None:
                continue
            if mine is not None and mine[i] is None:
                called = "no_call"
            elif mine is not None and mine[i] > 0:
                called = TRUE[mine[i]]
            elif any(v[i] for a, v in here.items() if a != alt and v[i] is not None):
                called = "other_allele"
            else:
                called = "hom_ref"
            table[s][(TRUE[t], called)] += 1
    for pos, by_alt in calls.items():
        if pos in planted:
            continue
        for i, s in enumerate(samples):
            if any(v[i] for v in by_alt.values() if v[i] is not None):
                false_variants[s] += 1
    return table, false_variants


def _derived(counts, false_variants):
    """The figures a wrong genotype prior moves, from one group's counts."""
    def c(t, k):
        return counts.get((t, k), 0)

    het = sum(c("het", k) for k in CALLED)
    hom_alt = sum(c("hom_alt", k) for k in CALLED)
    nonref_true = het + hom_alt
    right_nonref = c("het", "het") + c("hom_alt", "hom_alt")
    called_nonref_on_ref = c("hom_ref", "het") + c("hom_ref", "hom_alt")
    nonref_union = nonref_true + called_nonref_on_ref + false_variants
    return {
        "het_called_het": round(c("het", "het") / het, 5) if het else "",
        "hom_alt_called_hom_alt": round(c("hom_alt", "hom_alt") / hom_alt, 5) if hom_alt else "",
        "het_called_homozygous": c("het", "hom_ref") + c("het", "hom_alt"),
        "homozygote_called_het": c("hom_ref", "het") + c("hom_alt", "het"),
        "nonref_concordance": round(right_nonref / nonref_union, 5) if nonref_union else "",
    }


def score(truth_vcf, calls_vcf, reference, samples, group_of, base, out_tsv, min_qual=None):
    """Score one call set; write the per-sample table to `out_tsv` and return the rows."""
    truth = truth_genotypes(truth_vcf, samples)
    calls = call_genotypes(calls_vcf, reference, samples, min_qual=min_qual)
    table, false_variants = compare(truth, calls, samples)
    keys = [(t, k) for t in TRUE for k in CALLED]
    with open(out_tsv, "w") as fh:
        fh.write("sample\tgroup\t" + "\t".join(f"{t}>{k}" for t, k in keys)
                 + "\tfalse_variant_genotypes\n")
        for s in samples:
            fh.write(f"{s}\t{group_of(s)}\t" + "\t".join(str(table[s][k]) for k in keys)
                     + f"\t{false_variants[s]}\n")
    groups = defaultdict(list)
    for s in samples:
        groups[group_of(s)].append(s)
    rows = []
    for label, members in [(f"F = {g:.2f}", groups[g]) for g in sorted(groups)] + [("all", samples)]:
        counts = defaultdict(int)
        for s in members:
            for k, v in table[s].items():
                counts[k] += v
        fv = sum(false_variants[s] for s in members)
        cls = dict(base, variant_class=label, stage="genotype")
        for t, k in keys:
            rows.append(row(cls, metric=f"true {t}, called {k}", value=counts[(t, k)]))
        rows.append(row(cls, metric="false_variant_genotypes", value=fv))
        for metric, value in _derived(counts, fv).items():
            rows.append(row(cls, metric=metric, value=value))
    return rows


def realised_inbreeding(truth_vcf, snvs_tsv, samples, group_of):
    """{group: (F realised, heterozygotes, sites)}: 1 - observed / expected heterozygosity, the
    expectation 2p(1-p) from each SNP's planted frequency."""
    af = {}
    with open(snvs_tsv) as fh:
        next(fh)
        for line in fh:
            _, pos, _, _, _, p = line.rstrip("\n").split("\t")
            af[int(pos)] = float(p)
    truth = truth_genotypes(truth_vcf, samples)
    observed, expected, sites = defaultdict(int), defaultdict(float), defaultdict(int)
    for (pos, _, _), copies in truth.items():
        p = af.get(pos)
        if p is None:
            continue
        for i, s in enumerate(samples):
            if copies[i] is None:
                continue
            g = group_of(s)
            observed[g] += copies[i] == 1
            expected[g] += 2 * p * (1 - p)
            sites[g] += 1
    return {g: (round(1 - observed[g] / expected[g], 4), observed[g], sites[g])
            for g in expected if expected[g]}


def fitted_inbreeding(parameters_toml):
    """{sample: (F, warrant)} from a pvc parameters file's [inbreeding] table."""
    import re

    out = {}
    pattern = re.compile(r'\{ sample = "([^"]+)", inbreeding_coefficient = \{ value = ([^,]+),'
                         r'.*?warrant = "([^"]+)"')
    with open(parameters_toml) as fh:
        for line in fh:
            m = pattern.search(line)
            if m:
                out[m.group(1)] = (float(m.group(2)), m.group(3))
    return out

