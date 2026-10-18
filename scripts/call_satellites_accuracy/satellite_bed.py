"""Satellite ground of a typed-regions table, clipped to a BED of benchmark regions.

Usage: satellite_bed.py <typed.tsv> <regions.bed> <out.bed>

Both inputs are 0-based half-open. Writes each satellite's overlap with each region, sorted,
and prints how many spans and bases that is.
"""
import sys
from collections import defaultdict

typed, regions_path, out = sys.argv[1:4]

regions = defaultdict(list)
with open(regions_path) as handle:
    for line in handle:
        if line.startswith(("#", "track", "browser")) or not line.strip():
            continue
        chrom, start, end = line.split("\t")[:3]
        regions[chrom].append((int(start), int(end)))
for spans in regions.values():
    spans.sort()

pieces = []
with open(typed) as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        chrom, start, end, kind = line.split("\t")[:4]
        if kind != "satellite":
            continue
        start, end = int(start), int(end)
        for region_start, region_end in regions.get(chrom, []):
            if region_start >= end:
                break
            low, high = max(start, region_start), min(end, region_end)
            if low < high:
                pieces.append((chrom, low, high))

order = {}
pieces = sorted(set(pieces), key=lambda piece: (piece[0], piece[1]))
with open(out, "w") as handle:
    for chrom, low, high in pieces:
        handle.write(f"{chrom}\t{low}\t{high}\n")
print(f"{len(pieces)} spans, {sum(high - low for _, low, high in pieces)} bases")
