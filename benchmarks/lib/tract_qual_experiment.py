#!/usr/bin/env python3
"""Is a repeat tract's QUAL worth believing, can it be gated on, and are the
genotypes right?

`doc/devel/ng/spec/calling_loop_ssr.md` §3.3 asks two questions about the site
quality a caller writes at a repeat tract, and this program answers both on one
callset — plus a third the same files can answer and the site-level two cannot:

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

* **Genotype accuracy** — where the truth set and the caller both call a tract,
  how often the caller's genotype is the truth's. Compared as multisets of
  allele sequences, so two files that order their ALT columns differently still
  agree and phase is ignored. Written only when `--genotype-out` asks for it.

  **This is the question a discovery round moves and the other two barely see.**
  Admitting an allele that was hiding under stutter does not usually add a
  variant site — the site was already variant — it turns a homozygote into the
  heterozygote it really is. A calibration that counts sites cannot see that at
  all, and the threshold sweep sees only half of it: the allele appears in the
  ALT column, and whether the sample was given a copy of it is a different fact.

All three are split by **motif period** — a homopolymer against everything with a
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
        --calibration-out calib.tsv --sweep-out sweep.tsv --genotype-out gt.tsv

Every output file is appended to when it already carries a header, so a driver
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

    def intervals_by_contig(self) -> dict[str, list[TractInterval]]:
        """Every tract of the ground, by contig — what a batched reference read needs."""
        return self._by_contig

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
    samples: tuple[str, ...] = ()

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


    def genotype_indices(self, column: int) -> tuple[int, ...] | None:
        """This record's genotype at one sample column, as allele indices.

        `None` where the sample is not called (`./.`), where the record has no
        such column, or where an index names an allele the record does not
        carry — the last being a malformed record rather than a no-call.
        """
        if column >= len(self.samples):
            return None
        genotype = self.samples[column].split(":", 1)[0]
        if genotype in (".", "./.", ".|."):
            return None
        allele_count = 1 + len(self.alt.split(","))
        called = []
        for index in genotype.replace("|", "/").split("/"):
            if index == ".":
                return None
            try:
                position = int(index)
            except ValueError:
                return None
            if position >= allele_count:
                return None
            called.append(position)
        return tuple(called)

    def allele_bases(self, index: int) -> str:
        """The bases of one of this record's alleles, the reference at index 0."""
        if index == 0:
            return self.ref.upper()
        return self.alt.split(",")[index - 1].upper()


def tract_reference_bases(
    reference: Path, ground: TractGround, work: Path
) -> dict[tuple[str, int, int], tuple[int, str]]:
    """Every tract's reference bases, in one `samtools faidx` call.

    **The reference is opened because nothing else can settle the comparison.**
    Two sides describe one tract with records at different positions and over
    different spans — and the truth set writes a two-allele heterozygote as two
    phased records where the caller writes one multi-allelic record. The only
    representation both can be brought into is the tract's own sequence, and
    building that needs the bases between the records.

    One batched call rather than one a tract: 20,000 regions is 20,000
    subprocesses the other way.

    **The window reaches [`ANCHOR_PAD`] bases either side of the tract**, because
    a left-aligned insertion at a repeat's first base is anchored on the base
    *before* it. Without the pad every such record falls outside the window and
    the tract is scored as incomparable — measured, that was 3,354 of 3,648
    homopolymer tracts.

    Returned with each window's own first position rather than the tract's, so
    a tract at a contig's very start — where the pad is clipped — still indexes
    correctly instead of shifting every allele by one.
    """
    regions = work / "tract_regions.txt"
    keys: list[tuple[tuple[str, int, int], int]] = []
    with open(regions, "w", encoding="utf-8") as handle:
        for contig, intervals in sorted(ground.intervals_by_contig().items()):
            for one in intervals:
                # BED is half-open and zero-based, so `one.start` is already the
                # one-based position of the base *before* the tract.
                first = max(1, one.start + 1 - ANCHOR_PAD)
                keys.append(((contig, one.start, one.end), first))
                handle.write(f"{contig}:{first}-{one.end + ANCHOR_PAD}\n")
    completed = subprocess.run(
        ["samtools", "faidx", str(reference), "-r", str(regions)],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise SystemExit(f"samtools faidx failed:\n{completed.stderr.strip()}")
    bases: dict[tuple[str, int, int], tuple[int, str]] = {}
    current: list[str] = []
    index = -1
    for line in completed.stdout.splitlines():
        if line.startswith(">"):
            if index >= 0:
                key, first = keys[index]
                bases[key] = (first, "".join(current).upper())
            current = []
            index += 1
        else:
            current.append(line)
    if index >= 0:
        key, first = keys[index]
        bases[key] = (first, "".join(current).upper())
    return bases


def haplotypes_over_tract(
    records: list[VcfRecord],
    genotypes: list[tuple[int, ...]],
    tract: tuple[int, int],
    reference: str,
    phased: bool,
) -> set[tuple[str, ...]] | None:
    """Every tract sequence pair this side's records could describe.

    A **set** of pairs, because an unphased side with several records does not
    say which allele sits on which copy: every assignment is returned and the
    two sides agree if any pair is shared. A phased side returns the one pair
    its phase names.

    **Only the records carrying a non-reference allele on a copy are applied to
    it**, and that is what makes the common shape work rather than be refused.
    A truth set writes a two-allele heterozygote as two records at the same
    position — `AGT -> A` phased `0|1` and `AGTGTGT -> A` phased `1|0`. As
    edits they overlap and cannot both be applied; as haplotypes they do not,
    because each copy takes a non-reference allele from exactly one of them.
    Refusing on the records' spans alone threw out 1,412 of this benchmark's
    tracts, which are precisely the two-allele heterozygotes.

    `None` where a copy really does need two overlapping edits at once, or where
    a record reaches outside the window — a refusal rather than a guess.
    """
    start, end = tract
    for record in records:
        if record.pos < start or record.end > end:
            return None
    copies = len(genotypes[0])
    if any(len(one) != copies for one in genotypes):
        return None

    def sequence(assignment: list[int]) -> tuple[str, ...] | None:
        built = []
        for copy in range(copies):
            applied = [
                (record, genotypes[index][(copy + assignment[index]) % copies])
                for index, record in enumerate(records)
            ]
            applied = [(record, allele) for record, allele in applied if allele != 0]
            applied.sort(key=lambda pair: pair[0].pos)
            for (left, _), (right, _) in zip(applied, applied[1:]):
                if left.end >= right.pos:
                    return None
            out = []
            cursor = start
            for record, allele in applied:
                out.append(reference[cursor - start : record.pos - start])
                out.append(record.allele_bases(allele))
                cursor = record.end + 1
            out.append(reference[cursor - start :])
            built.append("".join(out))
        return tuple(sorted(built))

    if phased or len(records) == 1:
        one = sequence([0] * len(records))
        return None if one is None else {one}
    # Unphased with several records: every way of rotating each record's called
    # alleles across the copies, with the first held fixed so one assignment is
    # not counted twice.
    pairs: set[tuple[str, ...]] = set()
    for mask in range(copies ** (len(records) - 1)):
        assignment = [0] * len(records)
        remaining = mask
        for index in range(1, len(records)):
            assignment[index] = remaining % copies
            remaining //= copies
        one = sequence(assignment)
        if one is not None:
            pairs.add(one)
    return pairs or None


def sample_columns_of(path: Path) -> list[str]:
    """The sample names of a VCF, in column order."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:  # type: ignore[operator]
        for line in handle:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
            if not line.startswith("#"):
                break
    return []


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
                VcfRecord(
                    fields[0],
                    int(fields[1]),
                    fields[3],
                    fields[4],
                    qual,
                    fields[7],
                    tuple(fields[9:]),
                )
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


@dataclass
class GenotypeCell:
    """One (period class) cell of the genotype comparison.

    **The three wrong-answer counts are the point of the cell, not decoration.**
    A discovery round's whole purpose is to find an allele hiding under stutter,
    and a sample carrying one is called *homozygous where the truth is
    heterozygous* — so `called_homozygous_truth_heterozygous` is the number that
    says whether the mechanism is aimed at the failure this data actually has.
    Its opposite, `called_heterozygous_truth_homozygous`, is the damage a
    wrongly admitted allele does, and the same round can only make it worse.
    """

    tracts_truth_calls: int = 0
    tracts_both_call: int = 0
    genotype_right: int = 0
    no_call: int = 0
    not_comparable: int = 0
    truth_allele_never_offered: int = 0
    called_homozygous_truth_heterozygous: int = 0
    called_heterozygous_truth_homozygous: int = 0
    wrong_some_other_way: int = 0

    def score_haplotypes(
        self,
        called: set[tuple[str, ...]],
        expected: set[tuple[str, ...]],
        offered: frozenset[str],
    ) -> None:
        """Charge one tract to the right counter.

        Each side offers every tract-sequence pair its records could describe —
        one pair where the phase is known, several where it is not — and they
        agree if any pair is shared.

        **The four wrong-answer counters partition the errors by what would have
        to change to fix them**:

        - `truth_allele_never_offered` — a sequence the truth carries is not
          among the ones the caller's records name, so no genotype over that set
          could have been right. **This is candidate selection's**, and it is
          the ceiling on what a wider set can buy.
        - the three below are errors made over a set that *did* hold the right
          sequences, so they are the genotyper's: the likelihood and the prior
          picked the wrong pair from a set containing the right one.
        """
        self.tracts_both_call += 1
        if called & expected:
            self.genotype_right += 1
            return
        wanted = {sequence for pair in expected for sequence in pair}
        if any(sequence not in offered for sequence in wanted):
            self.truth_allele_never_offered += 1
            return
        homozygous_call = all(len(set(pair)) == 1 for pair in called)
        homozygous_truth = all(len(set(pair)) == 1 for pair in expected)
        if homozygous_call and not homozygous_truth:
            self.called_homozygous_truth_heterozygous += 1
        elif homozygous_truth and not homozygous_call:
            self.called_heterozygous_truth_homozygous += 1
        else:
            self.wrong_some_other_way += 1


def records_by_tract(
    records: list[VcfRecord], ground: TractGround
) -> dict[tuple[str, int, int, int], list[VcfRecord]]:
    """Every record a tract holds, by tract — **all of them, not the best one**.

    A tract can carry several records a side and they are not alternatives: the
    truth set writes a two-allele heterozygote as two phased records where the
    caller writes one multi-allelic record, and a tract can hold a substitution
    at one end and a length change at the other. Keeping one a tract compares
    different events and calls the result a genotype error.
    """
    out: dict[tuple[str, int, int, int], list[VcfRecord]] = {}
    for record in records:
        tract = ground.tract_at(record.contig, record.pos, record.end)
        if tract is None:
            continue
        key = (record.contig, tract.start, tract.end, tract.period)
        out.setdefault(key, []).append(record)
    return out


def compare_genotypes(
    query: list[VcfRecord],
    truth: list[VcfRecord],
    ground: TractGround,
    query_column: int,
    truth_column: int,
    tract_bases: dict[tuple[str, int, int], tuple[int, str]],
) -> dict[str, GenotypeCell]:
    """How often the caller says the tract holds what the truth says it holds.

    **This is the question a discovery round moves, and the site-level ones are
    not.** Admitting an allele that was hiding under stutter does not usually
    add a variant site — the site was already variant — it turns a homozygote
    into the heterozygote it really is. A calibration that counts sites cannot
    see that at all, and the threshold sweep sees only half of it.

    **Compared as the tract's own two sequences, rebuilt from the reference.**
    Nothing shorter works, and three attempts at something shorter each gave a
    wrong answer on this benchmark's own tandem-repeat ground:

    - **Allele strings, as written.** Two records describing one event over
      different spans — `AGT -> A` against `AGTGTGT -> A,AGTGT` — do not share a
      string, and `bcftools norm` does not bring them together because it trims
      each record against its own ALT column. 324 of 6,303 tracts read as
      genotype errors that were only a difference of spelling.
    - **Padding both records to their union span.** Fixes that pair and still
      fails wherever the two sides put their records at different places in the
      tract, because the span between them is in neither record's REF.
    - **One record a tract a side.** Throws away 1,711 tracts of 6,303, and the
      1,412 of them where the truth writes two phased records against the
      caller's one multi-allelic record are exactly the two-allele
      heterozygotes — the class most worth measuring.

    So each side's records are laid on the tract's reference bases and the two
    resulting sequences compared. A side whose records leave the phase open
    offers every assignment, and the two agree if any pair is shared.
    """
    truth_by_tract = records_by_tract(truth, ground)
    query_by_tract = records_by_tract(query, ground)
    cells: dict[str, GenotypeCell] = {}
    for key, truth_records in truth_by_tract.items():
        contig, start, end, period = key
        cell = cells.setdefault(period_class_of(period), GenotypeCell())
        truth_genotypes = [record.genotype_indices(truth_column) for record in truth_records]
        if any(one is None for one in truth_genotypes):
            continue
        cell.tracts_truth_calls += 1
        query_records = query_by_tract.get(key)
        if query_records is None:
            continue
        query_genotypes = [record.genotype_indices(query_column) for record in query_records]
        if any(one is None for one in query_genotypes):
            cell.no_call += 1
            continue
        window_bases = tract_bases.get((contig, start, end))
        if window_bases is None:
            cell.tracts_both_call += 1
            cell.not_comparable += 1
            continue
        first, reference = window_bases

        # **A truth record is phased where its own genotype separator says so**,
        # and this truth set writes `0|1`. An unphased side offers every
        # assignment instead.
        truth_phased = all("|" in record.samples[truth_column] for record in truth_records)
        window = (first, first + len(reference) - 1)
        expected = haplotypes_over_tract(
            truth_records,
            [one for one in truth_genotypes if one is not None],
            window,
            reference,
            truth_phased,
        )
        called = haplotypes_over_tract(
            query_records,
            [one for one in query_genotypes if one is not None],
            window,
            reference,
            False,
        )
        if expected is None or called is None:
            cell.tracts_both_call += 1
            cell.not_comparable += 1
            continue

        offered = frozenset(
            reference[: record.pos - window[0]]
            + record.allele_bases(index)
            + reference[record.end - window[0] + 1 :]
            for record in query_records
            for index in range(1 + len(record.alt.split(",")))
        )
        cell.score_haplotypes(called, expected, offered)
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

GENOTYPE_HEADER = (
    "arm\tground\tdepth\tsample\tperiod_class\ttracts_truth_calls\t"
    "tracts_both_call\tgenotype_right\tno_call\tgenotype_accuracy\t"
    "not_comparable\ttruth_allele_never_offered\tcalled_hom_truth_het\t"
    "called_het_truth_hom\twrong_some_other_way\n"
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
        "--genotype-out",
        type=Path,
        default=None,
        help="where to write the genotype comparison; skipped when not given",
    )
    parser.add_argument(
        "--genotype-sample",
        default=None,
        help="which of the query's sample columns to compare; its first by default",
    )
    parser.add_argument(
        "--genotype-truth-sample",
        default=None,
        help="the truth column holding the same individual, where the two files "
        "name it differently (production writes HG002_30x for the truth's HG002); "
        "the query's own name by default",
    )
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

        genotype_cells: dict[str, GenotypeCell] = {}
        genotype_note = "not asked for"
        if args.genotype_out is not None:
            query_samples = sample_columns_of(query_sites)
            truth_samples = sample_columns_of(truth_sites)
            wanted = args.genotype_sample or (query_samples[0] if query_samples else None)
            # **Both columns are named, never matched by position.** Comparing
            # column 0 with column 0 would silently score one individual's
            # genotypes against another's wherever two files order their samples
            # differently, and every number after that would look plausible. The
            # two names may differ — production writes `HG002_30x` where the
            # truth set writes `HG002` — so saying they are the same individual
            # is the caller's statement to make rather than this program's guess.
            in_truth = args.genotype_truth_sample or wanted
            if wanted is None:
                genotype_note = "the query names no sample"
            elif wanted not in query_samples:
                genotype_note = f"the query has no sample {wanted}"
            elif in_truth not in truth_samples:
                genotype_note = f"the truth has no sample {in_truth}"
            else:
                genotype_cells = compare_genotypes(
                    query_site_records,
                    truth_site_records,
                    ground,
                    query_samples.index(wanted),
                    truth_samples.index(in_truth),
                    tract_reference_bases(args.reference, ground, work),
                )
                genotype_note = (
                    f"query {wanted} against truth {in_truth}"
                    if in_truth != wanted
                    else f"sample {wanted}"
                )

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

    if args.genotype_out is not None:
        genotype_rows = []
        for period_class in sorted(genotype_cells):
            cell = genotype_cells[period_class]
            genotype_rows.append(
                f"{args.arm}\t{args.ground}\t{args.depth}\t{args.sample}\t"
                f"{period_class}\t{cell.tracts_truth_calls}\t{cell.tracts_both_call}\t"
                f"{cell.genotype_right}\t{cell.no_call}\t"
                f"{ratio(cell.genotype_right, cell.tracts_both_call)}\t"
                f"{cell.not_comparable}\t{cell.truth_allele_never_offered}\t"
                f"{cell.called_homozygous_truth_heterozygous}\t"
                f"{cell.called_heterozygous_truth_homozygous}\t"
                f"{cell.wrong_some_other_way}\n"
            )
        append_rows(args.genotype_out, GENOTYPE_HEADER, genotype_rows)
        print(f"  genotypes: {genotype_note}", file=sys.stderr)

    print(
        f"{args.arm} {args.ground} {args.depth} {args.sample}: "
        f"{sum(one.records for one in bins.values())} records on tract ground, "
        f"{off_ground} off it, {without_qual} without a QUAL",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
