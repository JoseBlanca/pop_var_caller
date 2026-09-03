#!/usr/bin/env python3
"""Which tracts changed verdict between two arms, and which way.

The tract-accuracy program's rule 3 (`doc/devel/ng/research/tract_accuracy_program.md`
§2): every arm reports its verdict FLIPS against baseline beside the headline,
because a headline can rise while correcting nothing — the one proposed
correction judged by its headline alone raised it 4.7 points and corrected not
one verdict. This joins two `--verdicts-out` dumps from
`tract_qual_experiment.py` on (contig, start, end) and prints:

* the crosstab of baseline verdict against arm verdict, per period class,
  off-diagonal cells only — the diagonal is what did not move;
* with `--list`, every flipped tract as one line, so a case can be pulled and
  read whole.

A tract present in one dump and absent from the other is a flip too (the truth
side of the comparison is the same ground, so absence means the truth's
genotype column changed between runs — which should not happen and is worth
seeing loudly), and is printed with the missing side as `absent`.

Usage:

    tract_verdict_flips.py --baseline base_verdicts.tsv --arm arm_verdicts.tsv
        [--baseline-arm LABEL] [--arm-arm LABEL] [--list]

The two files may be the same file when it carries several arms appended;
`--baseline-arm` / `--arm-arm` then say which rows are which.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path


def read_verdicts(path: Path, arm: str | None) -> dict[tuple[str, int, int], tuple[str, str]]:
    """(contig, start, end) -> (period_class, verdict), for one arm's rows."""
    out: dict[tuple[str, int, int], tuple[str, str]] = {}
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        wanted = {name: header.index(name) for name in
                  ("arm", "contig", "start", "end", "period_class", "verdict")}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if arm is not None and fields[wanted["arm"]] != arm:
                continue
            key = (fields[wanted["contig"]], int(fields[wanted["start"]]),
                   int(fields[wanted["end"]]))
            if key in out:
                raise SystemExit(
                    f"{path}: tract {key} appears twice — pass --baseline-arm/"
                    f"--arm-arm to pick one arm of a multi-arm file"
                )
            out[key] = (fields[wanted["period_class"]], fields[wanted["verdict"]])
    if not out:
        raise SystemExit(f"{path}: no rows" + (f" for arm {arm}" if arm else ""))
    return out


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--arm", type=Path, required=True)
    parser.add_argument("--baseline-arm", default=None,
                        help="which arm's rows to read from --baseline")
    parser.add_argument("--arm-arm", default=None,
                        help="which arm's rows to read from --arm")
    parser.add_argument("--list", action="store_true",
                        help="print every flipped tract, one a line")
    args = parser.parse_args()

    base = read_verdicts(args.baseline, args.baseline_arm)
    arm = read_verdicts(args.arm, args.arm_arm)

    flips: Counter[tuple[str, str, str]] = Counter()
    flipped: list[tuple[str, int, int, str, str, str]] = []
    for key in sorted(set(base) | set(arm)):
        period_class, was = base.get(key, (None, "absent"))
        arm_class, now = arm.get(key, (period_class, "absent"))
        if period_class is None:
            period_class = arm_class
        if was != now:
            flips[(period_class, was, now)] += 1
            flipped.append((*key, period_class, was, now))

    unchanged = len(set(base) & set(arm)) - sum(
        1 for one in flipped if one[4] != "absent" and one[5] != "absent"
    )
    print(f"tracts: {len(base)} baseline, {len(arm)} arm, "
          f"{unchanged} unchanged, {len(flipped)} flipped")
    for (period_class, was, now), count in sorted(
        flips.items(), key=lambda one: -one[1]
    ):
        print(f"  {period_class:12} {was:>14} -> {now:<14} {count}")
    if args.list:
        for contig, start, end, period_class, was, now in flipped:
            print(f"{contig}\t{start}\t{end}\t{period_class}\t{was}\t{now}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
