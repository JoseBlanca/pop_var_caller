#!/usr/bin/env python3
"""Split a dhat heap profile's peak into the parts of the caller that hold it.

`dhat-heap.json` records, for every allocation site, how many bytes that site
held at the instant the whole program's live heap was largest (`gb`, "bytes at
global max"). Summing that column gives the peak; grouping it by which part of
the caller allocated the bytes says where the peak actually is.

Sites are attributed by the **innermost project source file** in the stack, not
by symbol name. Every stack's innermost frames are the allocator shim and the
standard library's `Vec` growth path — the `__rust_alloc` frame carries the
profiled example's own path, so matching the crate name against the whole frame
string attributes every site in the program to whatever example was profiled.
Matching the file path against the project's own module directories is what
separates them.

    uv run tmp/milestone_z/attribute_peak.py <dhat-heap.json> [--top N]
"""

import argparse
import json
import re
import sys
from collections import defaultdict

# A dhat frame reads "0xADDR: symbol (path/to/file.rs:LINE:COL)".
FRAME = re.compile(r"^0x[0-9a-f]+:\s*(?P<symbol>.*?)\s*\((?P<path>[^()]*):\d+(?::\d+)?\)$")

# Which module a project source path belongs to, and what to call it. Order
# matters: the first pattern that matches the path wins.
# Some files hold more than one thing. `pipeline.rs` carries both the per-sample
# summary parse and the reference fetch behind the low-complexity filter, and
# they belong to different groups, so these symbol rules are checked against the
# whole stack before any path rule. Getting this wrong once made a reference-span
# copy look like sample metadata and inverted the comparison between two runs.
SYMBOLS = [
    ("reference / dust", ("sdust_mask_for_span", "dust_mask", "ReferenceFetcher",
                          "fetch_span", "sdust")),
    ("per-sample metadata", ("SampleSummary", "decompress_metadata", "from_toml_bytes")),
]

MODULES = [
    # Per-sample metadata: the .psp's own metadata section, decompressed and
    # parsed once per open file and held for the whole run. It scales with the
    # cohort size and with nothing else — not the block size, not the record
    # encoding, not the depth.
    ("per-sample metadata", ("src/psp/metadata", "src/psp/header")),
    ("block decode", ("src/psp/",)),
    ("per-sample columns", ("src/var_calling/sample_reader", "src/var_calling/from_psp/")),
    ("cohort chunk", ("src/var_calling/cohort_chunk", "src/var_calling/chunk",
                      "src/var_calling/producer")),
    ("merger", ("src/var_calling/per_group_merger", "src/var_calling/variant_caller",
                "src/var_calling/variant_grouping")),
    ("posterior", ("src/var_calling/posterior_engine",)),
    ("dust", ("src/var_calling/dust_filter",)),
    ("pileup record", ("src/pileup_record", "src/pileup/")),
    ("ng", ("src/ng/",)),
]

# Paths that are the project's own source rather than a dependency's. dhat
# reports std as `src/vec/mod.rs`, `alloc/src/...` — indistinguishable from ours
# by prefix alone, so ours are named explicitly.
PROJECT_DIRS = (
    "src/psp/", "src/var_calling/", "src/ng/", "src/pileup/", "src/pileup_record",
    "src/pop_var_caller/", "src/reference", "src/region", "src/vcf",
)


def parse(frame):
    match = FRAME.match(frame)
    if match:
        return match.group("symbol"), match.group("path")
    return frame, ""


def is_project(path):
    return any(d in path for d in PROJECT_DIRS)


def classify(frames):
    """Return (group, representative frame) for one allocation site."""
    parsed = [parse(f) for f in frames]

    # Symbol rules first: they cross file boundaries and settle the files that
    # hold more than one kind of thing.
    for group, needles in SYMBOLS:
        for symbol, path in parsed:
            if any(n in symbol for n in needles):
                inner = next((f"{s}  ({p})" for s, p in parsed if is_project(p)), symbol)
                return group, inner

    for frame in frames:
        symbol, path = parse(frame)
        if not is_project(path):
            continue
        for group, patterns in MODULES:
            if any(p in path for p in patterns):
                return group, f"{symbol}  ({path})"
        return "other project code", f"{symbol}  ({path})"
    # No project frame at all: a dependency allocating on its own behalf, or a
    # stack dhat could not walk.
    for frame in frames:
        symbol, path = parse(frame)
        if path and "dhat-" not in path and "/alloc/" not in path:
            return "dependencies", f"{symbol}  ({path})"
    return "unattributed", frames[0] if frames else "<empty>"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("profile")
    ap.add_argument("--top", type=int, default=15)
    args = ap.parse_args()

    doc = json.load(open(args.profile))
    ftbl, sites = doc["ftbl"], doc["pps"]

    by_group, by_site, total = defaultdict(int), defaultdict(int), 0
    for site in sites:
        held = site.get("gb", 0)
        if not held:
            continue
        total += held
        group, frame = classify([ftbl[i] for i in site["fs"]])
        by_group[group] += held
        by_site[(group, frame)] += held

    if not total:
        sys.exit("no bytes held at peak — is this a dhat heap profile?")

    print(f"# peak live heap {total / 1e6:.1f} MB, "
          f"{sum(1 for s in sites if s.get('gb', 0))} sites holding at peak")
    print()
    print(f"{'group':<22} {'MB at peak':>11} {'share':>7}")
    print("-" * 43)
    for group, held in sorted(by_group.items(), key=lambda kv: -kv[1]):
        print(f"{group:<22} {held / 1e6:>11.1f} {held / total:>6.1%}")
    print("-" * 43)
    print(f"{'total':<22} {total / 1e6:>11.1f} {1.0:>6.1%}")

    print(f"\ntop {args.top} sites")
    print("-" * 43)
    for (group, frame), held in sorted(by_site.items(), key=lambda kv: -kv[1])[: args.top]:
        print(f"{held / 1e6:>8.1f} MB {held / total:>6.1%}  [{group}] {frame}")


if __name__ == "__main__":
    main()
