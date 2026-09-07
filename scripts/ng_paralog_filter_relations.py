#!/usr/bin/env python3
"""**The three relations between the filter's files, checked on real reads.**

Companion to `ng_paralog_filter_runs.sh`, which calls the same cohort three ways — the filter off,
the filter dropping what it flags, and the filter tagging it instead. Spec §10 states relations
between those files, and until step D they had only ever been checked on a fixture cohort so small
that nothing was ever flagged, which makes them assertions about the empty case
(`doc/devel/ng/spec/hidden_paralog_filter.md` §10; step C5's report).

The three:

* **The filter changes no record it does not flag.** Take the dropping run, remove the two INFO
  keys the filter adds, and what is left must be the filter-off file with exactly the flagged
  records missing.
* **Dropping is tagging minus the tagged.** The records the tagging run did not flag must be, byte
  for byte, the records the dropping run wrote.
* **A flagged record is its own off-run record with the verdict spliced in.** Undo the two INFO
  keys and the `FILTER` splice and a tagged record must be byte for byte what the off run wrote.
  Without this the flagged records — the only place a ratio survives, and where every per-record
  number in the reports is read — are checked by nothing at all.

**What the tool refuses to do is pass by omission.** The two INFO keys and the four header lines
are stripped before comparing, so a filter that simply stopped writing them would satisfy every
relation above; they are therefore *asserted present* first, on the runs that must carry them, and
asserted absent on the run that must not. The same reasoning covers the record counts: a run that
wrote nothing satisfies everything.

Usage:  ng_paralog_filter_relations.py <off.vcf> <drop.vcf> <tag.vcf>

Exits non-zero if any relation or any precondition fails, and prints what differed when it does.
"""

import sys

# The header lines the filter adds, and the two INFO keys it adds to a scored record — spec §3.5.
FILTER_HEADER_PREFIXES = (
    "##paralogFilter=",
    "##INFO=<ID=PARALOG_LR,",
    "##INFO=<ID=PARALOG_POST,",
    "##FILTER=<ID=hiddenParalog,",
)
FILTER_INFO_KEYS = ("PARALOG_LR", "PARALOG_POST")
FILTER_ID = "hiddenParalog"

# What `patch.rs` puts in an empty column, and the two FILTER values it replaces rather than
# joins on (`append_filter_column`, `append_info_column`).
MISSING = "."
PASS = "PASS"

CHROM, POS, REF, ALT, QUAL, FILTER_COLUMN, INFO_COLUMN = 0, 1, 3, 4, 5, 6, 7


class Failed(Exception):
    """A relation or a precondition that did not hold. Carries what to print."""


def read_vcf(path):
    """The header lines and the record lines of a VCF, each kept whole and in order."""
    header, records = [], []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            (header if line.startswith("#") else records).append(line)
    return header, records


def without_the_filters_header_lines(header):
    """The header, less the four lines the filter adds and the one line no two runs share."""
    return [
        line
        for line in header
        if not line.startswith("##commandline=")
        and not line.startswith(FILTER_HEADER_PREFIXES)
    ]


def without_the_filters_info_keys(record):
    """One record line with `PARALOG_LR` and `PARALOG_POST` taken back out of its INFO column.

    A record whose INFO held nothing else goes back to `.`, which is what the writer would have
    put there — the filter replaces a `.` rather than appending to it (`patch.rs`).
    """
    columns = record.split("\t")
    fields = [
        field
        for field in columns[INFO_COLUMN].split(";")
        if field.split("=", 1)[0] not in FILTER_INFO_KEYS
    ]
    columns[INFO_COLUMN] = ";".join(fields) if fields else MISSING
    return "\t".join(columns)


def without_the_filters_id(record):
    """One record line with the `hiddenParalog` id taken back out of its FILTER column.

    The inverse of `patch.rs`'s `append_filter_column`: the id is joined on with `;` unless what
    was there was `PASS` or `.`, in which case it replaced the column outright. So removing the
    id leaves either the run's own filters or an empty column, and an empty one has to become
    `PASS` or `.` — which of the two is not recoverable from the tagged line alone, so this
    returns the column as `None` and lets the caller accept either.
    """
    columns = record.split("\t")
    ids = [one for one in columns[FILTER_COLUMN].split(";") if one != FILTER_ID]
    columns[FILTER_COLUMN] = ";".join(ids) if ids else None
    return columns


def is_flagged(record):
    return FILTER_ID in record.split("\t")[FILTER_COLUMN].split(";")


def carries_a_filter_info_key(record):
    keys = {field.split("=", 1)[0] for field in record.split("\t")[INFO_COLUMN].split(";")}
    return [key for key in FILTER_INFO_KEYS if key in keys]


def check_the_filters_own_lines_are_there(name, header, records, the_filter_ran):
    """**The filter's own marks, asserted rather than assumed.**

    Every relation below strips these before comparing, so a filter that stopped writing them
    would satisfy all three. This is what makes that a failure instead.
    """
    present = [line for line in header if line.startswith(FILTER_HEADER_PREFIXES)]
    if the_filter_ran and len(present) != len(FILTER_HEADER_PREFIXES):
        raise Failed(
            f"the {name} run's header carries {len(present)} of the filter's "
            f"{len(FILTER_HEADER_PREFIXES)} lines, not all of them:\n    "
            + "\n    ".join(present or ["(none)"])
        )
    if not the_filter_ran and present:
        raise Failed(
            f"the {name} run did not run the filter, yet its header carries "
            f"{len(present)} of the filter's lines"
        )

    with_both = with_one = 0
    for record in records:
        keys = carries_a_filter_info_key(record)
        if len(keys) == len(FILTER_INFO_KEYS):
            with_both += 1
        elif keys:
            with_one += 1
    if with_one:
        raise Failed(
            f"{with_one} record(s) of the {name} run carry one of "
            f"{', '.join(FILTER_INFO_KEYS)} and not the other; they go on together (spec §3.5)"
        )
    if the_filter_ran and with_both == 0 and records:
        raise Failed(
            f"the {name} run wrote {len(records)} record(s) and not one carries "
            f"{' or '.join(FILTER_INFO_KEYS)}, so no ratio reached the file"
        )
    if not the_filter_ran and with_both:
        raise Failed(
            f"the {name} run did not run the filter, yet {with_both} of its records carry "
            f"the filter's INFO keys"
        )
    return with_both


def compare(name, expected, actual, limit=5):
    """Two lists of lines, compared in order. Raises `Failed` naming the first few differences."""
    if len(expected) != len(actual):
        raise Failed(
            f"{name}: {len(expected)} line(s) expected, {len(actual)} found — "
            "the two runs did not write the same set of lines"
        )
    differing = [
        (at, one, two) for at, (one, two) in enumerate(zip(expected, actual)) if one != two
    ]
    if differing:
        shown = "\n".join(
            f"    line {at}:\n      expected {one[:200]}\n      found    {two[:200]}"
            for at, one, two in differing[:limit]
        )
        raise Failed(
            f"{name}: {len(differing)} of {len(expected)} line(s) differ\n{shown}"
            + (f"\n    … and {len(differing) - limit} more" if len(differing) > limit else "")
        )
    print(f"  {name}: HOLDS ({len(actual)} line(s) compared)")


def check_a_flagged_record_is_its_off_record_with_the_verdict_on_it(off_records, tag_records):
    """**Relation three** — the only thing that looks at a flagged record's other columns.

    Relations one and two both work by *removing* the flagged records from one side, so between
    them a flagged record's QUAL, its remaining INFO, its FORMAT and every sample column are
    compared against nothing.
    """
    checked = 0
    for off, tagged in zip(off_records, tag_records):
        if not is_flagged(tagged):
            continue
        checked += 1
        undone = without_the_filters_id(without_the_filters_info_keys(tagged))
        off_columns = off.split("\t")
        if undone[FILTER_COLUMN] is None:
            # The id replaced the column, so the off run's value must be one of the two the
            # splice replaces rather than joins on.
            if off_columns[FILTER_COLUMN] not in (PASS, MISSING):
                raise Failed(
                    f"a flagged record at {off_columns[CHROM]}:{off_columns[POS]} carries only "
                    f"{FILTER_ID}, but the off run filtered it as "
                    f"{off_columns[FILTER_COLUMN]!r} — the splice joins, it does not replace"
                )
            undone[FILTER_COLUMN] = off_columns[FILTER_COLUMN]
        if "\t".join(undone) != off:
            raise Failed(
                f"a flagged record is not its off-run record with the verdict spliced in, at "
                f"{off_columns[CHROM]}:{off_columns[POS]}:\n"
                f"      off    {off[:200]}\n      undone {chr(9).join(undone)[:200]}"
            )
    print(
        f"  a flagged record is its off-run record with the verdict spliced in: HOLDS "
        f"({checked} flagged record(s) compared whole)"
    )


def main(off_path, drop_path, tag_path):
    off_header, off_records = read_vcf(off_path)
    drop_header, drop_records = read_vcf(drop_path)
    tag_header, tag_records = read_vcf(tag_path)

    try:
        # **The runs are compared position by position, not by a key.** A key of
        # (CHROM, POS, REF, ALT) is not an identity: spec §3.2 says a repeat tract may share a
        # position with the generic locus owning its anchor base, and two records sharing all
        # four would make a keyed comparison wrong in both directions. All three runs call the
        # same cohort over the same ground and write in the same order, so the off run and the
        # tagging run are the same list of records and can be walked together.
        if len(off_records) != len(tag_records):
            raise Failed(
                f"the off run wrote {len(off_records)} record(s) and the tagging run "
                f"{len(tag_records)} — tagging removes nothing, so these must be equal"
            )

        with_keys = check_the_filters_own_lines_are_there(
            "tagging", tag_header, tag_records, the_filter_ran=True
        )
        check_the_filters_own_lines_are_there(
            "dropping", drop_header, drop_records, the_filter_ran=True
        )
        check_the_filters_own_lines_are_there("off", off_header, off_records, the_filter_ran=False)

        flagged = [record for record in tag_records if is_flagged(record)]
        print(
            f"  {len(off_records)} record(s) off, {len(drop_records)} kept by dropping, "
            f"{len(tag_records)} written by tagging, {len(flagged)} of them flagged; "
            f"{with_keys} record(s) carry a ratio"
        )

        # **Relation one.**
        compare(
            "the filter changed no record it did not flag",
            without_the_filters_header_lines(off_header)
            + [
                off
                for off, tagged in zip(off_records, tag_records)
                if not is_flagged(tagged)
            ],
            without_the_filters_header_lines(drop_header)
            + [without_the_filters_info_keys(record) for record in drop_records],
        )

        # **Relation two.**
        compare(
            "dropping is tagging minus the tagged",
            [record for record in tag_records if not is_flagged(record)],
            drop_records,
        )

        # **Relation three.**
        check_a_flagged_record_is_its_off_record_with_the_verdict_on_it(off_records, tag_records)
    except Failed as why:
        print(f"  FAILED — {why}")
        return 1

    if not flagged:
        print(
            "  NOTE: nothing was flagged, so the three relations are the empty case — they say "
            "the filter left the file alone, not that it removes what it should."
        )
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(*sys.argv[1:]))
