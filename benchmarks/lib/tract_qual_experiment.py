#!/usr/bin/env python3
"""Is a repeat tract's QUAL worth believing, and can it be gated on?

`doc/devel/ng/spec/calling_loop_ssr.md` §3.3 asks two questions about the site
quality a caller writes at a repeat tract, and this program answers both on one
callset:

* **Calibration** — bin the emitted records by QUAL and, in each bin, count how
  many sit at a repeat tract the truth set really does carry a variant in. A
  record written at QUAL 30 claims to be wrong about one time in a thousand; the
  bin says how often it actually is. It does **not** ask whether the alleles or
  the genotype are right, which are different questions with their own
  instruments.

  **The unit is the whole tract, not the record's own span**, because that is
  what the claim is about: a repeat tract's QUAL says the samples here are not
  all homozygous reference *at this tract*, and one length change can be spelled
  at either end of a compound repeat. Scored on the record's own span, a call
  and the truth record describing the same event nine bases away read as a false
  positive. So a record counts as truly variant when any truth record falls
  inside the tract it sits at, padded by one base.

* **Gateability** — sweep a QUAL threshold and report precision and recall at
  each step. This one is at the allele: a true positive is a truth record whose
  chromosome, position, REF and ALT a query record reproduces exactly, after
  both sides are left-aligned and multi-allelic records split. That is
  `benchmarks/giab/src/score_ng_recall.sh`'s rule, which is in turn
  `accuracy_dashboard.py`'s, so a number here is read against the ones already
  recorded for these callsets.

Both are split by **motif period** — a homopolymer against everything with a
repeat unit of two bases or more — because slippage rises steeply as the period
falls, and the risk §3.3 names is a model that under-prices slip products where
they are commonest.

**The ground is a BED of repeat tracts, not the caller's own annotation.** Two
callers put different flags on their records and one of them puts none at all,
so the tract ground is supplied from outside — `examples/ng_typed_region_dump`'s
`ssr_locus` rows — and every arm is scored on exactly the same intervals. The
BED's fourth column is the period. Records outside it are dropped and counted,
which is how a run that routed differently from the BED shows up rather than
silently scoring a different set of positions.

**Why a variant is allowed to sit one base outside its tract.** A left-aligned
insertion at the start of a repeat is anchored on the base *before* the repeat,
so its position falls outside the tract by one. Every overlap test here is
therefore padded by one base on each side. Without the pad an insertion at a
homopolymer's first base is scored as if it were not at the tract at all.

**Left-alignment runs before the region masks, and that is a departure from
`score_ng_recall.sh` worth stating.** That script masks first, and the two
orders do not treat the two sides alike at a confident interval's first base. A
truth record already written on the anchor base sits one base *outside* the
interval and is dropped; the query record describing the same event is not
left-aligned yet, so it is still inside, is kept, and only then moves onto the
anchor — where there is nothing left to match it. On this benchmark's
tandem-repeat ground, `chr1:69,233,430` is such a site: the truth set carries
`TATAATAATA -> T` there, the Tier interval starts at 69,233,431, and ng's
identical call was scored a false positive at QUAL 922. Left-aligning first
makes the two sides agree — at that site both are now dropped, because the
event's own position is outside the confident region and neither side should
speak for it.

**Every calibration row says whether the record is a repeat tract's.** A record
carrying the `STR` flag or a `PERIOD` key is a tract's, and everything else on
tract ground is `other` — a substitution inside a repeat, or an indel the SNP
path opened beside one. They are not the same question and the column keeps them
apart without throwing either away.

**A false positive that no sample was called with is counted apart.** A caller
may list an alternative allele in the ALT column and then give no sample a copy
of it; after `bcftools norm -m -any` that allele is its own record, and by the
scorer's rule it is a false positive. It is a different failure from a sample
genotyped wrong, so the sweep carries `fp_with_no_called_copy` beside `fp` —
records whose `AC` is zero — and the two are read together.

Usage:

    tract_qual_experiment.py --reference <ref.fa> --truth <truth.vcf.gz> \\
        --query <calls.vcf> --confident-bed <sample.bed> --tract-bed <tracts.bed> \\
        --arm ng --ground giab_per_sample --depth 30x --sample HG002 \\
        --calibration-out calib.tsv --sweep-out sweep.tsv

Both output files are appended to when they already carry a header, so a driver
loops over arms and depths and ends with one table of each.

`bcftools` does the left-alignment and the region restriction; everything else is
the standard library, so `uv run --no-project python` runs it with nothing
installed.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

# ---------------------------------------------------------------------------
# The bins and the ladder
# ---------------------------------------------------------------------------

# **Left-closed QUAL bin edges**, chosen so the decade a gate would be set in is
# not one lump. The last bin runs to infinity.
QUAL_BIN_EDGES: tuple[float, ...] = (0.0, 1.0, 3.0, 10.0, 20.0, 30.0, 50.0, 100.0, 200.0)

# The thresholds the sweep reports, as `QUAL >= t`. `0` is the ungated arm and
# is what a run's `.raw.vcf` holds; `30` is the floor the benchmark runners
# apply to every caller, so it is the row that reads against the standing
# numbers.
SWEEP_THRESHOLDS: tuple[float, ...] = (
    0.0, 1.0, 3.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0,
)

# How far outside its tract a record may sit and still be counted at it — see
# the module docstring's paragraph on the anchor base.
ANCHOR_PAD = 1


def error_probability_of(qual: float) -> float:
    """The chance QUAL claims this record is *not* a variant.

    Phred, so QUAL 30 is one in a thousand. A QUAL of zero claims nothing and
    comes back as 1.0; the `.` a caller writes when it has no number never
    reaches here, because such records are dropped when the VCF is read.
    """
    return 10.0 ** (-qual / 10.0)


# ---------------------------------------------------------------------------
# Reading the pieces
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TractInterval:
    """One repeat tract of the scored ground, half-open as BED writes it."""

    start: int
    end: int
    period: int


class TractGround:
    """The tract intervals, by contig, searchable by position.

    Built once and asked once per record, so the intervals are sorted at
    construction and looked up by binary search rather than scanned — a
    20,000-interval ground against 20,000 records is 4 x 10^8 comparisons the
    other way, for an answer that takes a second.
    """

    def __init__(self, path: Path) -> None:
        by_contig: dict[str, list[TractInterval]] = {}
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    raise SystemExit(
                        f"{path}: a tract-ground BED needs four columns "
                        f"(chrom, start, end, period); got {len(fields)}"
                    )
                by_contig.setdefault(fields[0], []).append(
                    TractInterval(int(fields[1]), int(fields[2]), int(fields[3]))
                )
        self._by_contig = {
            contig: sorted(intervals, key=lambda one: one.start)
            for contig, intervals in by_contig.items()
        }
        self._starts = {
            contig: [one.start for one in intervals]
            for contig, intervals in self._by_contig.items()
        }

    def tract_at(self, contig: str, start: int, end: int) -> TractInterval | None:
        """The tract this span touches, or `None` if none does.

        `start` and `end` are one-based inclusive, as a VCF record's own span
        is; the pad of the module docstring is applied here so no caller has to
        remember it. Where a padded span touches two tracts the rightmost of
        those starting before it is taken — a record spanning a tract boundary
        is one record and has to be charged to one of them.
        """
        intervals = self._by_contig.get(contig)
        if not intervals:
            return None
        low = start - 1 - ANCHOR_PAD  # to half-open, then padded
        high = end + ANCHOR_PAD
        index = bisect.bisect_right(self._starts[contig], high)
        for one in reversed(intervals[:index]):
            if one.start < high and one.end > low:
                return one
            if one.end <= low:
                # Sorted by start, so an earlier interval can still reach only
                # across a run of nested spans, which a typed-region dump does
                # not produce. Stop at the first that cannot.
                break
        return None


@dataclass(frozen=True)
class VcfRecord:
    """One VCF data line, reduced to what either question needs."""

    contig: str
    pos: int
    ref: str
    alt: str
    qual: float
    info: str

    @property
    def end(self) -> int:
        """The last reference base this record spans, one-based inclusive."""
        return self.pos + len(self.ref) - 1

    @property
    def key(self) -> tuple[str, int, str, str]:
        """What an allele-level match compares — the standing scorer's rule."""
        return (self.contig, self.pos, self.ref.upper(), self.alt.upper())

    @property
    def is_tract_record(self) -> bool:
        """Whether the caller says this record is a repeat tract's.

        ng writes the `STR` flag; production's repeat-tract caller writes
        `PERIOD`. A caller that annotates neither has every record read as
        `other`, which is the honest answer rather than a guess from the
        alleles' shape.
        """
        keys = {field.split("=", 1)[0] for field in self.info.split(";")}
        return "STR" in keys or "PERIOD" in keys

    @property
    def called_copies(self) -> int | None:
        """`AC` — how many chromosomes of the cohort carry this allele.

        `None` where the caller writes no `AC`, which is what a truth set does;
        the count that reads it treats an absent `AC` as "not known to be zero".
        """
        for field in self.info.split(";"):
            if field.startswith("AC="):
                try:
                    return int(float(field[3:].split(",")[0]))
                except ValueError:
                    return None
        return None


def read_vcf(path: Path, *, require_qual: bool) -> tuple[list[VcfRecord], int]:
    """Every data line of `path`, and how many were dropped for having no QUAL.

    A record whose QUAL is `.` states no claim, so it cannot be binned against
    one. Truth sets are read with `require_qual=False`, where the column is not
    a claim about anything and is not read.
    """
    records: list[VcfRecord] = []
    without_qual = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:  # type: ignore[operator]
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            qual_text = fields[5]
            if qual_text in (".", ""):
                if require_qual:
                    without_qual += 1
                    continue
                qual = 0.0
            else:
                qual = float(qual_text)
            records.append(
                VcfRecord(fields[0], int(fields[1]), fields[3], fields[4], qual, fields[7])
            )
    return records, without_qual


# ---------------------------------------------------------------------------
# bcftools: left-alignment, splitting, and the two region masks
# ---------------------------------------------------------------------------


def prepared_vcf(
    source: Path,
    *,
    reference: Path,
    masks: list[Path],
    split_multiallelics: bool,
    pass_only: bool,
    work: Path,
    name: str,
) -> Path:
    """`source` left-aligned and then restricted to every mask, under `work`.

    One stage at a time so a failure names the stage: drop everything but `PASS`
    where asked (the truth sets are read that way and the queries are not, which
    is what the standing scorer does), left-align against the reference, split
    multi-allelic records where the question is at the allele, and only then
    restrict to each region file in turn.

    **Not splitting is deliberate for the calibration half.** A multi-allelic
    record carries one QUAL, and splitting it would enter that one claim in the
    tally two or three times.
    """
    current = source
    if pass_only:
        nxt = work / f"{name}.pass.vcf.gz"
        run_bcftools(["view", "-f", "PASS", "-Oz", "-o", str(nxt), str(current)])
        index_vcf(nxt)
        current = nxt
    normed = work / f"{name}.norm.vcf.gz"
    arguments = ["norm", "-f", str(reference)]
    if split_multiallelics:
        arguments += ["-m", "-any"]
    arguments += ["-Oz", "-o", str(normed), str(current)]
    run_bcftools(arguments)
    index_vcf(normed)
    current = normed
    for index, mask in enumerate(masks):
        nxt = work / f"{name}.mask{index}.vcf.gz"
        run_bcftools(["view", "-T", str(mask), "-Oz", "-o", str(nxt), str(current)])
        index_vcf(nxt)
        current = nxt
    return current


def run_bcftools(arguments: list[str]) -> None:
    """One `bcftools` call, with its own stderr kept if it fails."""
    completed = subprocess.run(
        ["bcftools", *arguments], capture_output=True, text=True, check=False
    )
    if completed.returncode != 0:
        raise SystemExit(
            f"bcftools {' '.join(arguments)} failed:\n{completed.stderr.strip()}"
        )


def index_vcf(path: Path) -> None:
    run_bcftools(["index", "-f", "-t", str(path)])


# ---------------------------------------------------------------------------
# The two answers
# ---------------------------------------------------------------------------


def period_class_of(period: int) -> str:
    """The split §3.3 asks for: a homopolymer, or a unit of two bases or more.

    **Period 0 means the interval is not a repeat tract at all**, and comes back
    as `generic`. That is how ng's SNP and indel quality is binned on the same
    instrument as its tract quality: hand the ground BED the `generic` rows of a
    typed-region dump with 0 in the period column. §3.3's decision rule asks
    whether the tract quality reaches the standard the corrected SNP quality
    reaches on the same benchmark, and this is what makes the two comparable —
    one scorer, one truth set, one matching rule.
    """
    if period == 0:
        return "generic"
    return "homopolymer" if period == 1 else "period2plus"


class TruthSites:
    """The truth records' spans per contig, asked "is there a variant here?"."""

    def __init__(self, records: list[VcfRecord]) -> None:
        by_contig: dict[str, list[tuple[int, int]]] = {}
        for record in records:
            by_contig.setdefault(record.contig, []).append((record.pos, record.end))
        self._by_contig = {contig: sorted(spans) for contig, spans in by_contig.items()}
        self._starts = {
            contig: [span[0] for span in spans]
            for contig, spans in self._by_contig.items()
        }

    def covers(self, contig: str, start: int, end: int) -> bool:
        """Whether any truth record's span touches `[start, end]`, padded.

        The walk back from the insertion point is bounded, because truth spans
        are short — a GIAB record reaches tens of bases — and an unbounded walk
        over a contig's records would be quadratic. The bound is generous
        against the longest span these truth sets carry.
        """
        spans = self._by_contig.get(contig)
        if not spans:
            return False
        low = start - ANCHOR_PAD
        high = end + ANCHOR_PAD
        index = bisect.bisect_right(self._starts[contig], high)
        for span_start, span_end in reversed(spans[max(0, index - 256) : index]):
            if span_start <= high and span_end >= low:
                return True
        return False


@dataclass
class CalibrationBin:
    records: int = 0
    truly_variant: int = 0
    claimed_error: float = 0.0

    def add(self, is_true: bool, qual: float) -> None:
        self.records += 1
        self.truly_variant += int(is_true)
        self.claimed_error += error_probability_of(qual)


def bin_index_of(qual: float) -> int:
    """Which QUAL bin a record falls in — the last bin absorbs everything above."""
    return max(0, bisect.bisect_right(QUAL_BIN_EDGES, qual) - 1)


def bin_label(index: int) -> str:
    low = QUAL_BIN_EDGES[index]
    if index + 1 == len(QUAL_BIN_EDGES):
        return f"{low:g}+"
    return f"{low:g}-{QUAL_BIN_EDGES[index + 1]:g}"


def calibrate(
    query: list[VcfRecord], truth: TruthSites, ground: TractGround
) -> tuple[dict[tuple[str, str, int], CalibrationBin], int]:
    """Per (record kind, period class, QUAL bin): records, and how many were right.

    Returns the bins and the number of query records that fell on no tract of
    the ground, which the caller reports rather than absorbs.
    """
    bins: dict[tuple[str, str, int], CalibrationBin] = {}
    off_ground = 0
    for record in query:
        tract = ground.tract_at(record.contig, record.pos, record.end)
        if tract is None:
            off_ground += 1
            continue
        kind = "tract" if record.is_tract_record else "other"
        key = (kind, period_class_of(tract.period), bin_index_of(record.qual))
        # The tract's own span is the unit — see the module docstring. BED is
        # half-open and zero-based, so the first base is `start + 1`.
        bins.setdefault(key, CalibrationBin()).add(
            truth.covers(record.contig, tract.start + 1, tract.end), record.qual
        )
    return bins, off_ground


@dataclass
class SweepCell:
    tp: int = 0
    fp: int = 0
    fp_with_no_called_copy: int = 0
    fn: int = 0


def sweep(
    query: list[VcfRecord], truth: list[VcfRecord], ground: TractGround
) -> dict[tuple[str, float], SweepCell]:
    """Precision and recall's counts, per period class and threshold.

    Both sides are restricted to the tract ground first, so the denominator is
    the truth this experiment is about rather than the whole callset. A truth
    record and a query record match when their contig, position, REF and ALT are
    equal — both having been left-aligned and split by then.
    """
    truth_by_class: dict[str, set[tuple[str, int, str, str]]] = {}
    for record in truth:
        tract = ground.tract_at(record.contig, record.pos, record.end)
        if tract is None:
            continue
        truth_by_class.setdefault(period_class_of(tract.period), set()).add(record.key)

    query_by_class: dict[str, list[VcfRecord]] = {}
    for record in query:
        tract = ground.tract_at(record.contig, record.pos, record.end)
        if tract is None:
            continue
        query_by_class.setdefault(period_class_of(tract.period), []).append(record)

    cells: dict[tuple[str, float], SweepCell] = {}
    for period_class in set(truth_by_class) | set(query_by_class):
        truth_keys = truth_by_class.get(period_class, set())
        records = query_by_class.get(period_class, [])
        for threshold in SWEEP_THRESHOLDS:
            matched: set[tuple[str, int, str, str]] = set()
            false_positives = 0
            without_a_called_copy = 0
            for one in records:
                if one.qual < threshold:
                    continue
                if one.key in truth_keys:
                    matched.add(one.key)
                else:
                    false_positives += 1
                    if one.called_copies == 0:
                        without_a_called_copy += 1
            cells[(period_class, threshold)] = SweepCell(
                tp=len(matched),
                fp=false_positives,
                fp_with_no_called_copy=without_a_called_copy,
                fn=len(truth_keys) - len(matched),
            )
    return cells


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

CALIBRATION_HEADER = (
    "arm\tground\tdepth\tsample\trecord_kind\tperiod_class\tqual_bin\trecords\t"
    "truly_variant\tshare_truly_variant\tmean_claimed_error\n"
)

SWEEP_HEADER = (
    "arm\tground\tdepth\tsample\tperiod_class\tmin_qual\ttp\tfp\t"
    "fp_with_no_called_copy\tfn\tprecision\trecall\n"
)


def append_rows(path: Path, header: str, rows: list[str]) -> None:
    """Write `rows` to `path`, writing `header` first if the file is new."""
    fresh = not path.exists() or path.stat().st_size == 0
    with open(path, "a", encoding="utf-8") as handle:
        if fresh:
            handle.write(header)
        handle.writelines(rows)


def ratio(numerator: int, denominator: int) -> str:
    return f"{numerator / denominator:.4f}" if denominator else "."


# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--query", type=Path, required=True)
    parser.add_argument(
        "--confident-bed",
        type=Path,
        required=True,
        help="the sample's own confident regions; both sides are cut to it",
    )
    parser.add_argument(
        "--tract-bed",
        type=Path,
        required=True,
        help="chrom/start/end/period, from examples/ng_typed_region_dump's ssr_locus rows",
    )
    parser.add_argument("--arm", required=True, help="which caller and setting, e.g. ng")
    parser.add_argument("--ground", required=True, help="which dataset, e.g. giab_per_sample")
    parser.add_argument("--depth", required=True, help="e.g. 30x")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--calibration-out", type=Path, required=True)
    parser.add_argument("--sweep-out", type=Path, required=True)
    parser.add_argument(
        "--keep-work",
        type=Path,
        default=None,
        help="where to leave the intermediate VCFs; a temporary directory by default",
    )
    args = parser.parse_args()

    for path in (args.reference, args.truth, args.query, args.confident_bed, args.tract_bed):
        if not path.exists():
            raise SystemExit(f"missing input: {path}")

    with tempfile.TemporaryDirectory() as scratch:
        work = args.keep_work if args.keep_work else Path(scratch)
        work.mkdir(parents=True, exist_ok=True)
        masks = [args.confident_bed]

        def prepare(source: Path, name: str, *, split: bool, pass_only: bool) -> Path:
            return prepared_vcf(
                source,
                reference=args.reference,
                masks=masks,
                split_multiallelics=split,
                pass_only=pass_only,
                work=work,
                name=name,
            )

        truth_sites = prepare(args.truth, "truth_sites", split=False, pass_only=True)
        truth_alleles = prepare(args.truth, "truth_alleles", split=True, pass_only=True)
        query_sites = prepare(args.query, "query_sites", split=False, pass_only=False)
        query_alleles = prepare(args.query, "query_alleles", split=True, pass_only=False)

        ground = TractGround(args.tract_bed)
        truth_site_records, _ = read_vcf(truth_sites, require_qual=False)
        truth_allele_records, _ = read_vcf(truth_alleles, require_qual=False)
        query_site_records, without_qual = read_vcf(query_sites, require_qual=True)
        query_allele_records, _ = read_vcf(query_alleles, require_qual=True)

        bins, off_ground = calibrate(
            query_site_records, TruthSites(truth_site_records), ground
        )
        cells = sweep(query_allele_records, truth_allele_records, ground)

    calibration_rows = []
    for key in sorted(bins):
        kind, period_class, index = key
        bucket = bins[key]
        calibration_rows.append(
            f"{args.arm}\t{args.ground}\t{args.depth}\t{args.sample}\t{kind}\t"
            f"{period_class}\t{bin_label(index)}\t{bucket.records}\t"
            f"{bucket.truly_variant}\t{ratio(bucket.truly_variant, bucket.records)}\t"
            f"{bucket.claimed_error / bucket.records:.6f}\n"
        )
    append_rows(args.calibration_out, CALIBRATION_HEADER, calibration_rows)

    sweep_rows = []
    for key in sorted(cells):
        period_class, threshold = key
        cell = cells[key]
        sweep_rows.append(
            f"{args.arm}\t{args.ground}\t{args.depth}\t{args.sample}\t{period_class}\t"
            f"{threshold:g}\t{cell.tp}\t{cell.fp}\t{cell.fp_with_no_called_copy}\t"
            f"{cell.fn}\t{ratio(cell.tp, cell.tp + cell.fp)}\t"
            f"{ratio(cell.tp, cell.tp + cell.fn)}\n"
        )
    append_rows(args.sweep_out, SWEEP_HEADER, sweep_rows)

    print(
        f"{args.arm} {args.ground} {args.depth} {args.sample}: "
        f"{sum(one.records for one in bins.values())} records on tract ground, "
        f"{off_ground} off it, {without_qual} without a QUAL",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
