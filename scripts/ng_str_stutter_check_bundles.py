"""Check that the region typing calls every simulated STR bundle one bundle over the laid hull, and
that the caller places stutter where the simulator did.

    python3 scripts/ng_str_stutter_check_bundles.py <run-dir>

Reads, in <run-dir>: bundles.bed and bundle_runs.tsv (the simulator's bundles, and each one's shape
and reference stutter sites); typing_check.tsv (`ng_typed_region_dump <reference> <confident.bed>
calling`); and bundle_runs_check.tsv (`ng_str_bundle_runs <reference> <confident.bed>`, the caller's
reference runs of every typed bundle). scripts/ng_str_stutter_benchmark.sh writes the last two
before calling this. Exits 1, with the first few failures, unless every bundle passes and the typing
makes no bundle the simulator did not lay.

**A simulated bundle passes when:**

1. **exactly one typed `ssr_bundle` overlaps its hull, and starts and ends where the hull does**;
2. **in the separated shape** (pure members a few bases apart, the shape laid before hard shapes),
   **that bundle's members are the drawn members**: the same number, in order, each with the drawn
   motif length, with a start and an end each within TOLERANCE bases of the drawn member's, and —
   when the typing writes one — the drawn motif itself. A hard shape's members are not compared:
   the typing reports the repeats its scanner sees, which nest, overlap and join pieces, as on
   HG002, and the simulator's members are its pure pieces, which no typing lists;
3. **the caller's reference runs over that hull are the simulator's stutter sites on the reference**
   (bundle_runs.tsv), each the same motif, offset and length. The caller prices stutter in those
   runs (doc/devel/ng/spec/str_bundle_loci.md §6.1); the simulator stutters a separated bundle in
   each member and a hard one in the runs its own copy of that rule finds — so this says the two
   agree on where the reference slips. Skipped, saying so, when a run predates either file.

The typed members are the repeat detector's intervals cut back to whole copies of their motif, which
can differ from a drawn run by part of a unit when a gap base happens to continue the motif; the
simulator draws gap bases so that none does (`draw_gap` in examples/ng_str_stutter_benchmark.rs), so
the tolerance is 0 by default and a non-zero one is for measuring how far off a failing run is.

Env: BUNDLE_MEMBER_TOLERANCE (bases, default 0).
"""
import bisect
import collections
import json
import os
import sys

TOLERANCE = int(os.environ.get("BUNDLE_MEMBER_TOLERANCE", "0"))


def read_simulated(path):
    """Per contig, (start, end, name, members) with members [(start, end, period, motif)], zero-based
    half-open."""
    bundles = collections.defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, start, end, name, periods, motifs, _repeats, starts, ends = line.rstrip("\n").split("\t")
            members = [(int(s), int(e), int(p), m)
                       for s, e, p, m in zip(starts.split(","), ends.split(","), periods.split(","),
                                             motifs.split(","))]
            bundles[chrom].append((int(start), int(end), name, members))
    return bundles


def read_typed(path):
    """Per contig, typed bundles as (start, end, members), zero-based half-open, in start order. A
    member's motif is None when the typing predates writing one."""
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
            members = sorted((one["start"], one["end"], one["period"], one.get("motif"))
                             for one in json.loads(fields[8]))
            typed[fields[0]].append((int(fields[1]) - 1, int(fields[2]), members))
    for chrom in typed:
        typed[chrom].sort()
    return typed


def read_runs(path, key_columns):
    """{(chrom, zero-based hull start): the row's fields} from a runs table; None when it is absent."""
    if not os.path.exists(path):
        return None
    runs = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("bundle\t"):
                continue
            fields = line.rstrip("\n").split("\t")
            runs[key_columns(fields)] = fields
    return runs


def main():
    run = sys.argv[1]
    simulated = read_simulated(os.path.join(run, "bundles.bed"))
    typed = read_typed(os.path.join(run, "typing_check.tsv"))
    # The simulator's: bundle chrom hull_start hull_end shape runs (zero-based hull start).
    laid_runs = read_runs(os.path.join(run, "bundle_runs.tsv"), lambda f: (f[1], int(f[2])))
    # The caller's: chrom start end members runs (one-based start).
    caller_runs = read_runs(os.path.join(run, "bundle_runs_check.tsv"), lambda f: (f[0], int(f[1]) - 1))
    compare_runs = laid_runs is not None and caller_runs is not None
    shapes = collections.Counter()
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
            laid = laid_runs.get((chrom, start)) if laid_runs is not None else None
            if laid_runs is not None and laid is None:
                sys.exit(f"{name} {chrom}:{start}-{end} is in bundles.bed but not in bundle_runs.tsv")
            shape = laid[4] if laid is not None else "separated"
            # The caller's row: chrom start end members runs; None when the tool reported no bundle here.
            caller = caller_runs.get((chrom, start)) if caller_runs is not None else None
            shapes[shape] += 1
            if len(overlapping) != 1:
                outcome = "no_typed_bundle" if not overlapping else "several_typed_bundles"
            elif overlapping[0][1][:2] != (start, end):
                outcome = "other_hull"
            elif compare_runs and caller is None:
                outcome = "no_caller_runs_row"
            elif compare_runs and caller[4] != laid[5]:
                outcome = "other_runs"
            elif shape != "separated":
                outcome = "matched"
            else:
                members = overlapping[0][1][2]
                if len(members) != len(drawn) or any(m[2] != d[2] for m, d in zip(members, drawn)):
                    outcome = "other_members"
                elif any(m[3] is not None and m[3] != d[3] for m, d in zip(members, drawn)):
                    outcome = "other_motifs"
                else:
                    offset = max(max(abs(m[0] - d[0]), abs(m[1] - d[1])) for m, d in zip(members, drawn))
                    largest_offset = max(largest_offset, offset)
                    outcome = "matched" if offset <= TOLERANCE else "members_off_by_more_than_tolerance"
            outcomes[outcome] += 1
            if outcome != "matched":
                detail = ""
                if outcome == "other_runs":
                    detail = f" laid runs {laid[5]} caller's runs {caller[4]}"
                failures.append(f"{name} {chrom}:{start}-{end} {shape} drawn {drawn} -> {outcome} "
                                f"{[one for _, one in overlapping]}{detail}")
    invented = [(chrom, one) for chrom, rows in typed.items()
                for index, one in enumerate(rows) if (chrom, index) not in claimed]
    total = sum(outcomes.values())
    print(f"       typing: {outcomes['matched']} of {total} STR bundles typed as one bundle over the laid "
          f"hull ({', '.join(f'{count} {shape}' for shape, count in sorted(shapes.items()))}); the "
          f"separated ones with the drawn members (tolerance {TOLERANCE} bases; largest member-end "
          f"offset seen {largest_offset}); "
          + ("matched means the caller's reference runs equal the simulator's stutter sites too"
             if compare_runs else "reference runs NOT compared (bundle_runs.tsv or "
             "bundle_runs_check.tsv missing)")
          + f"; {len(invented)} typed bundles the simulator did not lay")
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
