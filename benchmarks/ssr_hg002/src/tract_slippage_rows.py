#!/usr/bin/env python3
"""Write the `slippage_by_stratum_and_group` rows of a parameters file, at a stated setting.

A repeat tract is scored under three numbers — how often a read reports a length
other than its allele's, which way the misreports go, and how fast two-repeat
misreports fall off against one-repeat ones. With `--defaults` a run has none of
them fitted and takes HipSTR's shipped 0.10 / 0.50 / 0.05, **one pair of numbers
for every class of tract**: a 30-copy mononucleotide run and a 5-copy
tetranucleotide are scored identically.

This writes those rows so a run can be given something else, which is how the
sweep in `doc/devel/reports/ng_tract_genotype_improvement_2026-09-02.md` §2 was
made. To use them: delete the run's own
`slippage_by_stratum_and_group = []` line and append this file to the end of its
`.parameters.toml`. **The end, not in place** — an array-of-tables closes the
table it sits in, so pasting these where the empty array was would make every
key after them a key of the last row.

    tract_slippage_rows.py --share 0.10 --shorter 0.50 --fall-off 0.05 --out rows.toml
    tract_slippage_rows.py --base 0.04 --slope 0.004 --out rows.toml

**Two shapes, and the second is the one that turned out to matter.** `--share`
gives every stratum one number. `--base`/`--slope` make it rise with tract
length, `clamp(base + slope * repeats / period, 0.01, 0.60)` — because slippage
rises steeply as a tract lengthens and as its motif period falls, and because a
flat change is a dial that trades a spurious heterozygote for a collapsed one
without moving the total (that report's §2).

**A row is written for every repeat count a candidate can reach, not just every
reference length.** A tract's scoring parameters are looked up by the
*candidate's* repeat count, so a candidate several repeats either side of the
reference needs a row of its own or it silently falls back to the caller's
default and the run stops being the setting it claims to be.
"""

from __future__ import annotations

import argparse

# The lowest and highest reference repeat count each motif period reaches in the
# HG002 tandem-repeat benchmark's tract ground, from
# `examples/ng_typed_region_dump`'s `ssr_locus` rows over its 50,000 intervals.
GROUND = {1: (8, 43), 2: (6, 50), 3: (6, 32), 4: (6, 24), 5: (5, 17), 6: (4, 16)}

# How far past the ground's own range a row is still written, so no candidate
# selection can offer falls off the end.
MARGIN = 12

ROW = """[[repeat_tracts.slippage_by_stratum_and_group]]
period = {period}
reference_repeats = {repeats}
slippage_group = 0
share_of_reads_that_slip = {share}
shorter_share = {shorter}
fall_off = {fall_off}
share_of_reads_that_slip_origin = {{ smoothing = "this_stratum", expected_slipped_reads = 0.0 }}
shorter_share_and_fall_off_origin = {{ expected_slipped_reads = 0.0, \
shorter_share_smoothing = "this_stratum", fall_off_smoothing = "this_stratum" }}

"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--share", type=float, help="one slip share for every stratum")
    parser.add_argument("--base", type=float, help="with --slope: the gradient's intercept")
    parser.add_argument("--slope", type=float, help="with --base: rise per repeat, per period")
    parser.add_argument("--shorter", type=float, default=0.50)
    parser.add_argument("--fall-off", type=float, default=0.05)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    if (args.share is None) == (args.base is None):
        parser.error("give either --share, or --base with --slope")
    if args.base is not None and args.slope is None:
        parser.error("--base needs --slope")

    def share_at(period: int, repeats: int) -> float:
        if args.share is not None:
            return args.share
        return min(0.60, max(0.01, args.base + args.slope * repeats / period))

    with open(args.out, "w", encoding="utf-8") as handle:
        if args.share is not None:
            handle.write(f"# slip share {args.share} at every stratum,")
        else:
            handle.write(
                f"# slip share clamp({args.base} + {args.slope} * repeats / period, 0.01, 0.60),"
            )
        handle.write(f" shorter {args.shorter}, fall-off {args.fall_off}\n")
        rows = 0
        for period, (low, high) in GROUND.items():
            for repeats in range(max(1, low - MARGIN), high + MARGIN + 1):
                handle.write(
                    ROW.format(
                        period=period,
                        repeats=repeats,
                        share=f"{share_at(period, repeats):.4f}",
                        shorter=args.shorter,
                        fall_off=args.fall_off,
                    )
                )
                rows += 1
    print(f"{rows} rows -> {args.out}")


if __name__ == "__main__":
    main()
