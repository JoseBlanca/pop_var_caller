"""How much of HG002's tandem-repeat truth falls in STR bundles, and what bundles look like.

    python3 scripts/ng_str_bundles/bundle_census.py <typed.tsv> <truth.vcf.gz>

Typed regions are 1-based inclusive (start, end) as the dump writes them. A truth record is placed in
the typed region its first changed base falls in (POS+1, after the padding base). Reported:
truth variant records by region kind, and bundles by member count and length.

<typed.tsv> is `ng_typed_region_dump <reference> <regions.bed> calling` over the benchmark's regions
(`scripts/ng_str_bundles/type_hg002.sh` makes it). Measured on main at a87ce7b6, 2026-10-11: 2,422
truth records start in 2,473 bundles, against 8,534 in 20,183 STR segments.
"""
import bisect
import collections
import gzip
import sys

typed_path, truth_path = sys.argv[1], sys.argv[2]
regions = collections.defaultdict(list)
for line in open(typed_path):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    regions[f[0]].append((int(f[1]), int(f[2]), f[3], f[6]))
for chrom in regions:
    regions[chrom].sort()
starts = {c: [r[0] for r in rs] for c, rs in regions.items()}


def kind_at(chrom, pos):
    rs = regions.get(chrom)
    if not rs:
        return None
    i = bisect.bisect_right(starts[chrom], pos) - 1
    if i >= 0 and rs[i][0] <= pos <= rs[i][1]:
        return rs[i][2]
    return None


by_kind = collections.Counter()
het_hom = collections.Counter()
with gzip.open(truth_path, "rt") as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        f = line.split("\t", 10)
        kind = kind_at(f[0], int(f[1]) + 1) or "outside_typed"
        by_kind[kind] += 1
print("truth variant records by the typed region their first changed base is in:")
for k, v in by_kind.most_common():
    print(f"  {k}: {v}")

members = collections.Counter()
lengths = []
for rs in regions.values():
    for start, end, kind, copies in rs:
        if kind == "ssr_bundle":
            members[int(copies) if copies.isdigit() else copies] += 1
            lengths.append(end - start + 1)
lengths.sort()
n = len(lengths)
print(f"bundles: {n}; members per bundle: {dict(sorted(members.items()))}")
print(f"bundle length: median {lengths[n // 2]}, 90th percentile {lengths[int(n * 0.9)]}, "
      f"max {lengths[-1]}")
