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

    uv run scripts/attribute_peak.py <dhat-heap.json> [--top N]
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
#
# The rows for production's per-sample file, calling engine and pileup — and the
# symbol rules that split `var_calling/pipeline.rs` between two of them — went
# with that code in promotion Milestone D, and the single `src/ng/` row became
# one row per stage when Milestone E moved ng's modules up to `src/`.
MODULES = [
    ("psp", ("src/psp/",)),
    ("cohort merge and run", ("src/run/",)),
    ("calling", ("src/calling/", "src/genetics")),
    ("parameter fit", ("src/parameter_estimation/",)),
    ("locus generation", ("src/locus_generation/",)),
    ("reads and alignment", ("src/read/", "src/alignment/")),
    ("paralog filter", ("src/paralog/", "src/window_coverage/")),
    ("reference and repeats", ("src/ref_seq", "src/raw_chrom_reader", "src/reference_info",
                               "src/tandem_repeat", "src/repeat_catalog", "src/region_typing",
                               "src/segmentation_inputs")),
    ("vcf", ("src/vcf/",)),
    ("read input", ("src/bam/", "src/fasta/")),
]

# Paths that are the project's own source rather than a dependency's. dhat
# reports std as `src/vec/mod.rs`, `alloc/src/...` — indistinguishable from ours
# by prefix alone, so ours are named explicitly.
PROJECT_DIRS = (
    "src/psp/", "src/run/", "src/calling/", "src/genetics", "src/parameter_estimation/",
    "src/locus_generation/", "src/read/", "src/alignment/", "src/paralog/",
    "src/window_coverage/", "src/ref_seq", "src/raw_chrom_reader", "src/reference_info",
    "src/tandem_repeat", "src/repeat_catalog", "src/region", "src/segmentation_inputs",
    "src/vcf/", "src/types", "src/bam/", "src/fasta/", "src/cli/",
)


def parse(frame):
    match = FRAME.match(frame)
    if match:
        return match.group("symbol"), match.group("path")
    return frame, ""


def is_project(path):
    # The module fragments below also occur inside dependencies' own trees — noodles-sam has a
    # `src/alignment/`, arrow a `src/types.rs` — so a frame from the cargo registry, a git
    # checkout or the standard library is never the project's, whatever its path contains.
    if "/registry/src/" in path or "/git/checkouts/" in path or "/rustc/" in path:
        return False
    return any(d in path for d in PROJECT_DIRS)


def classify(frames):
    """Return (group, representative frame) for one allocation site."""
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
