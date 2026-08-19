# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Does setting the STR prior's total concentration to the panel's gene diversity
make the prior assert the diversity that was measured?

calling_priors.md §5.1 sets `Σα` = the cohort's STR gene diversity D. D is a
probability (the chance two random copies differ); `Σα` is a count of chromosomes.
This checks what the prior then implies.

For a Dirichlet with concentration α, total A = Σα and Simpson index
c = Σ (α_a/A)², the chance two copies drawn from it differ is

    D_implied = A·(1 − c) / (A + 1)

so `A = D` gives D·(1 − c)/(D + 1), and the A that would reproduce D is
D/(1 − c − D), which has no solution once D ≥ 1 − c.

Per locus: allele frequencies come from the called genotypes (REPCN), the shape
comes from G₀'s geometric decay in repeat-unit offset from the modal length
(src/ssr/cohort/allele_freq_prior.rs).
"""

import sys
from collections import Counter

VCF = "benchmarks/ssr_tomato1/results_ssr15k/ours/cohort/cohort.ssr.vcf"
G0_FLOOR = 1e-12
DECAYS = [0.3, 0.5, 0.7]  # 0.5 is DEFAULT_G0_FALLBACK_P
MIN_CALLED = 10  # samples with a genotype, else the frequencies are noise


def loci(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[6] != "PASS":
                continue
            fmt = f[8].split(":")
            if "REPCN" not in fmt:
                continue
            i = fmt.index("REPCN")
            counts = Counter()
            called = 0
            for s in f[9:]:
                parts = s.split(":")
                if len(parts) <= i or parts[i] == "." or parts[i] == "":
                    continue
                try:
                    cns = [int(x) for x in parts[i].split(",")]
                except ValueError:
                    continue
                called += 1
                counts.update(cns)
            if called >= MIN_CALLED and len(counts) >= 1:
                yield f[0], f[1], counts


def simpson(weights):
    total = sum(weights)
    return sum((w / total) ** 2 for w in weights)


def main():
    rows = []
    for contig, pos, counts in loci(VCF):
        n = sum(counts.values())
        freqs = {a: c / n for a, c in counts.items()}
        d_measured = 1.0 - sum(f * f for f in freqs.values())
        modal = counts.most_common(1)[0][0]
        row = {"locus": f"{contig}:{pos}", "n_alleles": len(counts), "D": d_measured}
        for p in DECAYS:
            w = [max(p ** abs(a - modal), G0_FLOOR) for a in counts]
            c = simpson(w)
            row[f"c@{p}"] = c
            row[f"impl@{p}"] = d_measured * (1 - c) / (d_measured + 1)
            row[f"need@{p}"] = (
                d_measured / (1 - c - d_measured) if (1 - c - d_measured) > 0 else None
            )
        rows.append(row)

    if not rows:
        sys.exit("no loci passed the filters")

    poly = [r for r in rows if r["n_alleles"] > 1]
    print(f"PASS loci with >= {MIN_CALLED} called samples : {len(rows)}")
    print(f"  of those, more than one allele length      : {len(poly)}")
    print()

    def q(vals, f):
        v = sorted(vals)
        return v[min(len(v) - 1, int(f * len(v)))]

    ds = [r["D"] for r in poly]
    print("measured gene diversity D over polymorphic loci")
    print(f"  median {q(ds, 0.5):.3f}   p90 {q(ds, 0.9):.3f}   max {max(ds):.3f}")
    print()

    for p in DECAYS:
        impl = [r[f"impl@{p}"] for r in poly]
        ratio = [r[f"impl@{p}"] / r["D"] for r in poly if r["D"] > 0]
        need = [r[f"need@{p}"] for r in poly if r[f"need@{p}"] is not None]
        nosol = sum(1 for r in poly if r[f"need@{p}"] is None)
        cs = [r[f"c@{p}"] for r in poly]
        print(f"decay p = {p}   (G0 shape: median Simpson c = {q(cs, 0.5):.3f})")
        print(
            f"  Sigma-alpha = D  implies diversity   median {q(impl, 0.5):.3f}"
            f"   p90 {q(impl, 0.9):.3f}"
        )
        print(
            f"  implied / measured                   median {q(ratio, 0.5):.3f}"
            f"   p10 {q(ratio, 0.1):.3f}   p90 {q(ratio, 0.9):.3f}"
        )
        if need:
            print(
                f"  Sigma-alpha that would reproduce D   median {q(need, 0.5):.2f}"
                f"   p90 {q(need, 0.9):.2f}"
            )
        print(f"  loci where no Sigma-alpha reproduces D (D >= 1 - c): {nosol}")
        print()


if __name__ == "__main__":
    main()


def paired():
    """The paired need/D factor, which is 1/(1 - c - D) and not a ratio of medians."""
    for p in DECAYS:
        vals = []
        for contig, pos, counts in loci(VCF):
            n = sum(counts.values())
            d = 1.0 - sum((c / n) ** 2 for c in counts.values())
            if len(counts) < 2:
                continue
            modal = counts.most_common(1)[0][0]
            w = [max(p ** abs(a - modal), G0_FLOOR) for a in counts]
            c = simpson(w)
            if 1 - c - d > 0:
                vals.append(1 / (1 - c - d))
        vals.sort()
        print(
            f"decay {p}: need/D per locus  median {vals[len(vals)//2]:.2f}"
            f"  p90 {vals[int(0.9*len(vals))]:.2f}"
        )
