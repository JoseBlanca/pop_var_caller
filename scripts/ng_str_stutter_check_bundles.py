"""Check that the region typing calls every simulated STR bundle one bundle, with the drawn members.

    python3 scripts/ng_str_stutter_check_bundles.py <run-dir>

Reads <run-dir>/bundles.bed (the simulator's bundles) and <run-dir>/typing_check.tsv
(`ng_typed_region_dump <reference> <confident.bed> calling`, which scripts/ng_str_stutter_benchmark.sh
writes before calling this). Exits 1, with the first few failures, unless every bundle passes and
the typing makes no bundle the simulator did not lay.

**A simulated bundle passes when exactly one typed `ssr_bundle` overlaps its hull, and that bundle's
members are the drawn members**: the same number of them, in order, each with the drawn motif
length and with a start and an end each within TOLERANCE bases of the drawn member's. The typed
members are the repeat detector's intervals, which can differ from a drawn run by part of a unit
when a gap base happens to continue the motif; the simulator draws gap bases so that none does
(`draw_gap` in examples/ng_str_stutter_benchmark.rs), so the tolerance is 0 by default and a
non-zero one is for measuring how far off a failing run is.

Env: BUNDLE_MEMBER_TOLERANCE (bases, default 0).
"""
import bisect
import collections
import json
import os
import sys

TOLERANCE = int(os.environ.get("BUNDLE_MEMBER_TOLERANCE", "0"))


def read_simulated(path):
    """Per contig, (start, end, name, members) with members [(start, end, period)], zero-based half-open."""
    bundles = collections.defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, start, end, name, periods, _motifs, _repeats, starts, ends = line.rstrip("\n").split("\t")
            members = [(int(s), int(e), int(p))
                       for s, e, p in zip(starts.split(","), ends.split(","), periods.split(","))]
            bundles[chrom].append((int(start), int(end), name, members))
    return bundles


def read_typed(path):
    """Per contig, typed bundles as (start, end, members), zero-based half-open, in start order."""
    typed = collections.defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[3] != "ssr_bundle":
                continue
            if len(fields) < 9:
                sys.exit(f"{path}: an ssr_bundle row has no members column; the typing dump was "
                         "built before it wrote one (examples/ng_typed_region_dump.rs)")
            members = sorted((one["start"], one["end"], one["period"]) for one in json.loads(fields[8]))
            typed[fields[0]].append((int(fields[1]) - 1, int(fields[2]), members))
    for chrom in typed:
        typed[chrom].sort()
    return typed


def main():
    run = sys.argv[1]
    simulated = read_simulated(os.path.join(run, "bundles.bed"))
    typed = read_typed(os.path.join(run, "typing_check.tsv"))
    starts = {chrom: [one[0] for one in rows] for chrom, rows in typed.items()}
    outcomes = collections.Counter()
    failures = []
    largest_offset = 0
    claimed = set()
    for chrom, rows in simulated.items():
        for start, end, name, drawn in rows:
            candidates = typed.get(chrom, [])
            # Typed bundles start-sorted; one overlapping this hull starts less than 1,000 bases
            # before it (a repeat of 1,000 bases is a satellite, not a bundle member).
            overlapping = []
            for index in range(bisect.bisect_left(starts.get(chrom, []), start - 1000), len(candidates)):
                one = candidates[index]
                if one[0] >= end:
                    break
                if one[1] > start:
                    overlapping.append((index, one))
            claimed.update((chrom, index) for index, _ in overlapping)
            if len(overlapping) != 1:
                outcome = "no_typed_bundle" if not overlapping else "several_typed_bundles"
            else:
                members = overlapping[0][1][2]
                if len(members) != len(drawn) or any(m[2] != d[2] for m, d in zip(members, drawn)):
                    outcome = "other_members"
                else:
                    offset = max(max(abs(m[0] - d[0]), abs(m[1] - d[1])) for m, d in zip(members, drawn))
                    largest_offset = max(largest_offset, offset)
                    outcome = "matched" if offset <= TOLERANCE else "members_off_by_more_than_tolerance"
            outcomes[outcome] += 1
            if outcome != "matched":
                failures.append(f"{name} {chrom}:{start}-{end} drawn {drawn} -> {outcome} "
                                f"{[one for _, one in overlapping]}")
    invented = [(chrom, one) for chrom, rows in typed.items()
                for index, one in enumerate(rows) if (chrom, index) not in claimed]
    total = sum(outcomes.values())
    print(f"       typing: {outcomes['matched']} of {total} STR bundles typed as one bundle with the "
          f"drawn members (tolerance {TOLERANCE} bases; largest member-end offset seen "
          f"{largest_offset}); {len(invented)} typed bundles the simulator did not lay")
    for outcome, count in sorted(outcomes.items()):
        if outcome != "matched":
            print(f"       {outcome}: {count}")
    if failures or invented:
        for failure in failures[:10]:
            print("!! " + failure, file=sys.stderr)
        for chrom, one in invented[:10]:
            print(f"!! typed bundle not laid: {chrom}:{one[0]}-{one[1]} {one[2]}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
