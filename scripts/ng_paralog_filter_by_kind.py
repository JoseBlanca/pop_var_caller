#!/usr/bin/env python3
"""**What the hidden-duplication filter scored and flagged, by kind of record** — plan step D4.

Spec §3.2 scores every record, and every record on both its signals except a repeat tract, which is
scored on coverage alone because slippage smears the allele split there. This counts what that comes
to on a real run: how many records of each kind were scored, how many were flagged, and where each
kind's likelihood ratios sit.

**Deletions are counted apart from repeat tracts** even though both span several bases, because
their depth is biased in opposite directions: a tract's observation depth runs below its read depth,
while a deletion depresses depth across its own span in exactly the samples that carry it. Lumping
them would average one bias against the other and show neither (spec §3.2, plan step D4).

**Equal-length multi-base substitutions are their own kind.** They are neither an insertion nor a
deletion and classing them as a deletion — which a `len(ALT) <= len(REF)` test does — inflates the
deletion count and hides a kind whose depth is not biased at all.

Usage:  ng_paralog_filter_by_kind.py <tagged.vcf> [label]

Reads the **tagging** run, since a dropped record leaves no ratio in the file.
"""

import re
import sys
from collections import defaultdict

FILTER_ID = "hiddenParalog"
CHROM, POS, REF, ALT, FILTER_COLUMN, INFO_COLUMN = 0, 1, 3, 4, 6, 7

# The order kinds are reported in: the four that carry both signals, then the one scored on
# coverage alone, so a reader can see the coverage-only arm against the rest.
KINDS = [
    "biallelic SNP",
    "multiallelic",
    "insertion",
    "deletion",
    "equal-length substitution",
    "repeat tract",
]


def kind_of(reference, alternatives, info):
    """Which of the six a record is. **The tract flag wins**, because that is what the filter
    itself keys on — `is_repeat_tract` decides which signals the record's samples carry."""
    if "STR" in info.split(";"):
        return "repeat tract"
    alts = alternatives.split(",")
    if len(alts) > 1:
        return "multiallelic"
    alt = alts[0]
    if len(reference) == 1 and len(alt) == 1:
        return "biallelic SNP"
    if len(alt) > len(reference):
        return "insertion"
    if len(alt) < len(reference):
        return "deletion"
    return "equal-length substitution"


def percentile(sorted_values, share):
    """The value at sorted index `int(share * n)` — a value some record really has.

    **On fewer than a hundred records the ninety-ninth percentile is the maximum**, because
    `int(0.99 * n)` is the last index for every `n <= 100`. The table marks such cells rather than
    printing two columns that are one measurement.
    """
    return sorted_values[min(len(sorted_values) - 1, int(share * len(sorted_values)))]


def main(path, label):
    ratios = defaultdict(list)
    flagged = defaultdict(int)
    unscored = defaultdict(int)
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            kind = kind_of(columns[REF], columns[ALT], columns[INFO_COLUMN])
            found = re.search(r"PARALOG_LR=([^;\t]*)", columns[INFO_COLUMN])
            if found is None:
                unscored[kind] += 1
            else:
                ratios[kind].append(float(found.group(1)))
            if FILTER_ID in columns[FILTER_COLUMN].split(";"):
                flagged[kind] += 1

    written = {kind: len(ratios[kind]) + unscored[kind] for kind in KINDS}
    total = sum(written.values())
    print(f"=== {label} — {total} records ===")
    print(
        f"{'kind':<27}{'written':>8}{'scored':>8}{'flagged':>8}"
        f"{'lowest':>10}{'median':>10}{'p90':>10}{'p99':>10}{'highest':>10}"
    )
    for kind in KINDS:
        if written[kind] == 0:
            continue
        values = sorted(ratios[kind])
        if values:
            # `int(0.99 * n)` is the last index for n <= 100, so p99 is the maximum there and the
            # cell is marked rather than printed as if it were a separate measurement.
            p99 = f"{percentile(values, .99):>10.2f}" if len(values) > 100 else f"{'= max':>10}"
            spread = (
                f"{values[0]:>10.2f}{percentile(values, .5):>10.2f}"
                f"{percentile(values, .9):>10.2f}{p99}{values[-1]:>10.2f}"
            )
        else:
            spread = f"{'—':>10}" * 5
        print(f"{kind:<27}{written[kind]:>8}{len(values):>8}{flagged[kind]:>8}{spread}")

    # **The coverage-only arm, on its own.** A tract carries no read counts, so its allele term is
    # zero under every hypothesis and the ratio rests on coverage. Where the carrier hypothesis is
    # far enough away in σ₀ units to underflow, the ratio stops depending on the coverage at all
    # and answers the prior — the same number for every tract in a band. **How many share the
    # modal value is the measure of how much of the population the arm cannot see.**
    tracts = sorted(ratios["repeat tract"])
    if tracts:
        # The file prints four decimals, so "the same value" can only ever mean "the same to four
        # decimals"; the count within a tolerance says whether that is a rounding coincidence.
        rounded = [round(v, 4) for v in tracts]
        modal = max(set(rounded), key=rounded.count)
        share = rounded.count(modal)
        near = sum(1 for x in tracts if abs(x - modal) <= 1e-3)
        print(f"\nrepeat tracts: {len(tracts)} scored, highest ratio {tracts[-1]:.4f}")
        if share == 1:
            # **Said, not left to be inferred from a count of one.** Reporting "1 of them share one
            # ratio" reads as a repetition where there is none, which is how a first draft of D4
            # came to claim a coincidence that does not exist.
            print("  no two tracts share a ratio: the arm responds to coverage over this cohort")
        else:
            print(
                f"  {share} of them print the same ratio ({modal}), and {near} are within 0.001 "
                f"of it — the arm returns one answer for that whole band"
            )


if __name__ == "__main__":
    if not 2 <= len(sys.argv) <= 3:
        print(__doc__)
        sys.exit(2)
    main(sys.argv[1], sys.argv[2] if len(sys.argv) == 3 else sys.argv[1])
