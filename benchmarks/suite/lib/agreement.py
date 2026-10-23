"""Genotype agreement between two callers' cohort VCFs, with no truth (plan §4.7).

Both VCFs are normalised the same way — left-aligned against the reference, multi-allelic records
split, complex records decomposed (`bcftools norm -m -any -a`) — so a variant both callers found
has one spelling. The unit counted is **one sample's genotype at one variant**: every (sample,
variant) where at least one of the two callers gives that sample a copy of the allele. Each is one
of:

| outcome         | meaning                                                    |
|-----------------|------------------------------------------------------------|
| `same_genotype` | both give the sample the same number of copies (1 or 2)    |
| `other_genotype`| both give it the allele, in different numbers of copies    |
| `only_a`        | caller A gives it the allele, B gives none (or no record)  |
| `only_b`        | caller B gives it the allele, A gives none (or no record)  |
| `no_call_a` / `no_call_b` | one caller gives it the allele, the other leaves the sample uncalled |

A site one caller does not write at all counts as reference there, which is what a joint caller's
silence means. Split into SNPs and indels by the allele lengths.
"""

import subprocess

from .results import row

OUTCOMES = ["same_genotype", "other_genotype", "only_a", "only_b", "no_call_a", "no_call_b"]


def _genotypes(vcf, reference, regions, samples, min_qual=None):
    """{(chrom, pos, ref, alt): {sample: copies or None}} after normalising."""
    qual = f"-i 'QUAL>={min_qual}'" if min_qual is not None else ""
    cmd = (f"bcftools view --no-version -f PASS,. -T {regions} -s {','.join(samples)} {qual} {vcf} "
           f"| bcftools norm --no-version -f {reference} -m -any -a "
           f"| bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n'")
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    calls = {}
    for line in out.splitlines():
        f = line.split("\t")
        key = (f[0], int(f[1]), f[2], f[3])
        if f[3] == "*":
            continue
        per_sample = calls.setdefault(key, {})
        for sample, gt in zip(samples, f[4:]):
            alleles = gt.replace("|", "/").split("/")
            if "." in alleles:
                copies = None
            else:
                copies = sum(1 for a in alleles if a == "1")
            # A split record can appear twice for one key (two input records); keep the larger.
            old = per_sample.get(sample, 0)
            per_sample[sample] = copies if old is None or copies is None else max(old, copies)
    return calls


def compare(vcf_a, vcf_b, reference, regions, samples, min_qual_a=None, min_qual_b=None):
    """Return {variant class: {outcome: count}}."""
    a = _genotypes(vcf_a, reference, regions, samples, min_qual_a)
    b = _genotypes(vcf_b, reference, regions, samples, min_qual_b)
    counts = {"SNP": dict.fromkeys(OUTCOMES, 0), "indel": dict.fromkeys(OUTCOMES, 0)}
    for key in set(a) | set(b):
        cls = "SNP" if len(key[2]) == 1 and len(key[3]) == 1 else "indel"
        for s in samples:
            ca = a.get(key, {}).get(s, 0)
            cb = b.get(key, {}).get(s, 0)
            if not ca and not cb:
                continue
            if ca and cb:
                counts[cls]["same_genotype" if ca == cb else "other_genotype"] += 1
            elif ca and cb is None:
                counts[cls]["no_call_b"] += 1
            elif cb and ca is None:
                counts[cls]["no_call_a"] += 1
            elif ca:
                counts[cls]["only_a"] += 1
            else:
                counts[cls]["only_b"] += 1
    return counts


def rows(base, counts, label_a, label_b):
    out = []
    names = {"only_a": f"only_{label_a}", "only_b": f"only_{label_b}",
             "no_call_a": f"no_call_{label_a}", "no_call_b": f"no_call_{label_b}"}
    for cls, c in counts.items():
        total = sum(c.values())
        for outcome in OUTCOMES:
            out.append(row(base, variant_class=cls, stage=f"agreement with {label_b}",
                           metric=names.get(outcome, outcome), value=c[outcome]))
        out.append(row(base, variant_class=cls, stage=f"agreement with {label_b}",
                       metric="sample_variants_compared", value=total))
        out.append(row(base, variant_class=cls, stage=f"agreement with {label_b}",
                       metric="share_same_genotype",
                       value=round(c["same_genotype"] / total, 5) if total else ""))
    return out
