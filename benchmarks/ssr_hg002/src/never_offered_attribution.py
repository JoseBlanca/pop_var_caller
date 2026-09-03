#!/usr/bin/env python3
"""Where each never-offered truth sequence was lost — P2 of the tract-accuracy program.

At 463 of the baseline's 834 wrong tracts, a sequence the truth carries is not
among the ones the caller's record can spell (`never_offered` in
`tract_qual_experiment.py --verdicts-out`) — the genotype was decided before the
model was consulted. This joins those tracts against the caller's own candidate
machinery to say where each was lost, replacing the split produced by
`ng_tract_candidate_recall.py`, whose two recorded window defects
(`tract_genotype_accuracy_2026-09-03.md` §3.5) this derivation avoids by reusing
the corrected instrument's reconstruction on one side and `select_ssr`'s own
answers on the other.

Inputs:

* the same truth/query/ground the instrument scores, from which the missing
  sequences are re-derived exactly as the verdict dump derived them — the tract
  set is asserted equal to the dump's `never_offered` rows;
* `ng_candidate_selection_probe`'s per-tract dump (`NG_TRACT_DUMP=...`): one row
  per (tract, merge-table allele) carrying the allele's tract-span bases,
  whether any sample's reads cleared the support bar for it, and whether
  selection kept it. Its contig column is an index into the reference's contig
  order, resolved here against the reference `.fai`.

Each tract goes to exactly one class, by the furthest-downstream loss among its
missing sequences — the cheapest recovery, which is what sizes a lever:

  merge_refused    no locus was built; selection never ran
  tabled_kept      the sequence survived selection and the record still cannot
                   express it — real disagreements (§3.4c), counted apart so
                   they cannot inflate any lever's ceiling
  top_ploidy_cut   in the table, cleared the support bar, dropped by the
                   per-sample top-ploidy rung cut — the ONLY class a discovery
                   round (L7) can reach
  support_bar      in the table; no sample's reads cleared the bar for it
  never_tabled     no read carried the spelling into the merge — subdivided by
                   whether the truth's LENGTH is in the table (a spelling loss,
                   the realigner's territory) or not (absent from the evidence
                   as the merge saw it)

With `--bam-low` (the alignment the callset was made from), the length-absent
tracts get one further column: how many raw-aligned spanning reads carry the
missing truth length (`spurious_read_provenance.spanning_reads`'s rule). Reads
that carry it while the merge's table does not are an admission or alignment
loss — the realigner's, from the recall side; zero carriers is a limit of the
evidence, not of the caller.

Usage:

    never_offered_attribution.py --reference ref.fa --truth truth.vcf.gz \\
        --query calls.vcf --confident-bed tier_sorted.bed --tract-bed tier.bed \\
        --verdicts verdicts.tsv --verdicts-arm baseline --verdicts-depth 30x \\
        --candidates p2_candidates.tsv --per-tract-out out.tsv
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from tract_qual_experiment import (  # noqa: E402
    TractGround,
    haplotype_pairs,
    offered_sequences,
    period_class_of,
    prepared_vcf,
    read_vcf,
    records_near_each_tract,
    sample_columns_of,
    tract_reference_bases,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from spurious_read_provenance import spanning_reads  # noqa: E402


def read_candidate_dump(
    path: Path, contig_names: list[str]
) -> dict[tuple[str, int, int], dict]:
    """Per tract: whether it was built, and each tabled allele's fate.

    Dump coordinates are the tract's own bases, one-based inclusive — the same
    interval the instrument reads a genotype out over.
    """
    out: dict[tuple[str, int, int], dict] = {}
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        at = {name: header.index(name) for name in
              ("contig", "start", "end", "built", "verdict", "kept",
               "cleared_bar", "bases")}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            key = (contig_names[int(fields[at["contig"]])],
                   int(fields[at["start"]]), int(fields[at["end"]]))
            tract = out.setdefault(key, {"built": False, "alleles": {}})
            if fields[at["built"]] == "0":
                continue
            tract["built"] = True
            bases = fields[at["bases"]].upper()
            tract["alleles"][bases] = (
                fields[at["kept"]] == "1", fields[at["cleared_bar"]] == "1"
            )
    return out


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--query", type=Path, required=True)
    parser.add_argument("--confident-bed", type=Path, required=True)
    parser.add_argument("--tract-bed", type=Path, required=True)
    parser.add_argument("--verdicts", type=Path, required=True)
    parser.add_argument("--verdicts-arm", required=True)
    parser.add_argument("--verdicts-depth", required=True)
    parser.add_argument("--candidates", type=Path, required=True,
                        help="ng_candidate_selection_probe's NG_TRACT_DUMP file")
    parser.add_argument("--genotype-sample", default="HG002")
    parser.add_argument("--per-tract-out", type=Path, required=True)
    parser.add_argument("--bam-low", type=Path, default=None,
                        help="count raw carriers of the missing length at the "
                        "length-absent tracts")
    args = parser.parse_args()

    expected_tracts = set()
    with open(args.verdicts, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        column = {name: header.index(name) for name in
                  ("arm", "depth", "contig", "start", "end", "verdict")}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if (fields[column["arm"]] == args.verdicts_arm
                    and fields[column["depth"]] == args.verdicts_depth
                    and fields[column["verdict"]] == "never_offered"):
                expected_tracts.add((fields[column["contig"]],
                                     int(fields[column["start"]]),
                                     int(fields[column["end"]])))
    if not expected_tracts:
        raise SystemExit(f"{args.verdicts}: no never_offered rows for "
                         f"{args.verdicts_arm} at {args.verdicts_depth}")

    contig_names = [line.split("\t")[0]
                    for line in open(f"{args.reference}.fai", encoding="utf-8")]
    candidates = read_candidate_dump(args.candidates, contig_names)

    with tempfile.TemporaryDirectory() as scratch:
        work = Path(scratch)
        query_path = prepared_vcf(
            args.query, reference=args.reference, masks=[args.confident_bed],
            split_multiallelics=False, pass_only=False, work=work, name="query",
        )
        truth_path = prepared_vcf(
            args.truth, reference=args.reference, masks=[args.confident_bed],
            split_multiallelics=False, pass_only=True, work=work, name="truth",
        )
        query_records, _ = read_vcf(query_path, require_qual=True)
        truth_records, _ = read_vcf(truth_path, require_qual=False)
        query_column = sample_columns_of(query_path).index(args.genotype_sample)
        truth_column = sample_columns_of(truth_path).index(args.genotype_sample)
        ground = TractGround(args.tract_bed)
        bases_of = tract_reference_bases(args.reference, ground, work)
    truth_near = records_near_each_tract(truth_records, ground)
    query_near = records_near_each_tract(query_records, ground)

    found: dict[tuple[str, int, int], tuple[int, list[str]]] = {}
    for key, truths in truth_near.items():
        contig, start, end, period = key
        queries = query_near.get(key)
        if queries is None:
            continue
        truth_genotypes = [one.genotype_indices(truth_column) for one in truths]
        query_genotypes = [one.genotype_indices(query_column) for one in queries]
        if any(one is None for one in truth_genotypes):
            continue
        if any(one is None for one in query_genotypes):
            continue
        window = bases_of.get((contig, start, end))
        if window is None:
            continue
        first, reference = window
        tract = (start + 1, end)
        expected = haplotype_pairs(
            truths, truth_genotypes, first, reference,
            all("|" in one.samples[truth_column] for one in truths), tract,
        )
        called = haplotype_pairs(
            queries, query_genotypes, first, reference, False, tract
        )
        if expected is None or called is None or expected & called:
            continue
        offered = offered_sequences(queries, first, reference, tract)
        wanted = {one for pair in expected for one in pair}
        missing = sorted(one for one in wanted if one not in offered)
        if missing:
            found[(contig, start, end)] = (period, missing)

    if set(found) != expected_tracts:
        gone = expected_tracts - set(found)
        extra = set(found) - expected_tracts
        raise SystemExit(
            f"the derived tracts do not match the verdict dump: "
            f"{len(gone)} missing, {len(extra)} extra"
        )
    print(f"control: {len(found)} tracts, exactly the dump's never_offered rows",
          file=sys.stderr)

    tally: Counter = Counter()
    rows = []
    unjoined = 0
    for (contig, start, end), (period, missing) in sorted(found.items()):
        period_class = period_class_of(period)
        entry = candidates.get((contig, start + 1, end))
        if entry is None:
            unjoined += 1
            tally[(period_class, "no_dump_entry")] += 1
            rows.append((contig, start, end, period_class, len(missing),
                         "no_dump_entry", ".", "."))
            continue
        if not entry["built"]:
            klass, detail = "merge_refused", "."
        else:
            fates = [entry["alleles"].get(one) for one in missing]
            if any(one is not None and one[0] for one in fates):
                klass, detail = "tabled_kept", "."
            elif any(one is not None and one[1] for one in fates):
                klass, detail = "top_ploidy_cut", "."
            elif any(one is not None for one in fates):
                klass, detail = "support_bar", "."
            else:
                klass = "never_tabled"
                tabled_lengths = {len(one) for one in entry["alleles"]}
                detail = ("length_tabled"
                          if all(len(one) in tabled_lengths for one in missing)
                          else "length_absent")
        tally[(period_class, klass)] += 1
        raw_carriers = "."
        if klass == "never_tabled":
            tally[(period_class, f"  {detail}")] += 1
            if detail == "length_absent" and args.bam_low is not None:
                reads = spanning_reads(args.bam_low, contig, start + 1, end)
                absent = [len(one) for one in missing
                          if len(one) not in {len(two) for two in entry["alleles"]}]
                raw_carriers = max(
                    sum(1 for one in reads if one.tract_length == length)
                    for length in absent
                )
                tally[(period_class, "    raw carriers 0" if raw_carriers == 0
                       else "    raw carriers 1" if raw_carriers == 1
                       else "    raw carriers 2+")] += 1
        rows.append((contig, start, end, period_class, len(missing), klass,
                     detail, raw_carriers))

    with open(args.per_tract_out, "w", encoding="utf-8") as handle:
        handle.write("contig\tstart\tend\tperiod_class\tmissing_sequences\t"
                     "class\tdetail\traw_carriers\n")
        for row in rows:
            handle.write("\t".join(str(one) for one in row) + "\n")

    print(f"\n=== the {len(found)} never-offered tracts, attributed")
    for period_class in ("homopolymer", "period2plus"):
        for klass in ("merge_refused", "tabled_kept", "top_ploidy_cut",
                      "support_bar", "never_tabled", "  length_tabled",
                      "  length_absent", "    raw carriers 0",
                      "    raw carriers 1", "    raw carriers 2+",
                      "no_dump_entry"):
            count = tally[(period_class, klass)]
            if count or not klass.startswith(" "):
                print(f"  {period_class:12s} {klass:18s} {count:5d}")
    if unjoined:
        print(f"\n!! {unjoined} tracts had no dump entry — the join is incomplete",
              file=sys.stderr)
    print(f"\nper-tract table: {args.per_tract_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
