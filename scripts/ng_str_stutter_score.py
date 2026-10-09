"""Score one configuration of the stuttering benchmark (doc/devel/ng/spec/str_stutter_benchmark.md §6).

    python3 scripts/ng_str_stutter_score.py <run-dir>

Run it inside the dev container, after scripts/ng_str_stutter_benchmark.sh has filled <run-dir>.

## What is scored

**The unit is one sample at one STR segment.** The truth is the sample's two sequences over the
segment's window (the segment and 30 bases either side, `segment_alleles.tsv` and
`segment_genotypes.tsv`). A caller's answer is rebuilt from its records overlapping that window: each
record's called alleles are applied to the reference, and every way of assigning a heterozygous
record's two alleles to the two copies is tried, because no caller here phases reliably. **Two answers
are the same genotype when some assignment gives the truth's two sequences** — the same as
`tract_dual_path.md` §6's test, edits slid left and counted per copy, since two sets of edits that
spell the same sequence are the same edits once slid left. So an STR record and two generic records
saying the same thing score alike.

Each (caller, sample, segment) cell gets one outcome:

| outcome | when |
|---|---|
| `right` | the called sequences are the truth's |
| `length_right` | not the truth's sequences, but the same two lengths (whatever the bases) |
| `wrong` | a variant in the truth, called with other lengths |
| `called_reference` | a variant in the truth, every record over the window called homozygous reference |
| `false_call` | no variant in the truth, a variant called |
| `no_call` | a record over the window leaves the sample uncalled (`./.`), or the truth has a variant and no record covers the window |
| `ref_right` | no variant in the truth, none called |

A record that reaches outside the window widens the stretch compared, and the truth over the widening
is the reference: the simulator puts no variant within 150 bases of a segment beyond its window.
Overlapping records that would apply two edits to one base of one copy (GATK writes `*` for "deleted
by an upstream record", which is skipped) leave the cell `wrong`, or `false_call` on reference truth.

Every variant call carries its genotype quality, so the confidence of each caller's errors can be read
off: GQ for ng and GATK; for HipSTR, its posterior probability of the genotype written on the Phred
scale, capped at 99 as the others are. Where a cell takes several records, the smallest GQ among the
records whose call is not homozygous reference.

**SNPs and indels**: RTG vcfeval per sample per caller against the truth VCF, once inside the
segments' windows and once outside them (the rest of the genome). HipSTR inside only.

## What it writes into <run-dir>

- `score_cells.tsv` — one row per caller, sample and segment.
- `score_summary.tsv` — counts by caller, motif length, repeat band and outcome.
- `score_vcfeval.tsv` — vcfeval's counts per caller, sample and stretch.
"""

import collections
import gzip
import math
import os
import subprocess
import sys

WINDOW_FLANK_BP = 30
CALLING_FLOOR = [8, 6, 6, 6, 5, 4]
PLATEAU_REPEATS = [16, 16, 10, 10, 10, 10]
CALLERS = [
    ("ng_defaults", "ng_defaults.vcf"),
    ("ng_fitted", "ng_fitted.vcf"),
    ("ng_true_stutter", "ng_true_stutter.vcf"),
    ("gatk", "gatk.vcf.gz"),
    ("hipstr", "hipstr.vcf.gz"),
    ("hipstr_fixed_stutter", "hipstr_fixed_stutter.vcf.gz"),
]
HIPSTR = {"hipstr", "hipstr_fixed_stutter"}
RTG_JAR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tmp", "tools",
                       "rtg-tools-3.13", "RTG.jar")


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    contigs, name, chunks = {}, None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    contigs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        contigs[name] = "".join(chunks)
    return contigs


def repeat_band(period, repeats):
    """Where a reference repeat count sits on the stutter curve: below its middle, rising to the
    plateau, or on it (spec §3.2)."""
    floor, plateau = CALLING_FLOOR[period - 1], PLATEAU_REPEATS[period - 1]
    middle = (floor + plateau) // 2
    if repeats < middle:
        return f"{floor}-{middle - 1}"
    if repeats < plateau:
        return f"{middle}-{plateau - 1}"
    return f"{plateau}+"


class Segment:
    def __init__(self, name, chrom, window_start, window_end, period, motif, reference_repeats):
        self.name, self.chrom = name, chrom
        self.window_start, self.window_end = window_start, window_end
        self.period, self.motif, self.reference_repeats = period, motif, reference_repeats
        self.alleles = {}  # index -> (sequence, change)


def read_truth(run):
    segments = {}
    with open(os.path.join(run, "segment_alleles.tsv")) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            segment = segments.get(row["segment"])
            if segment is None:
                segment = segments[row["segment"]] = Segment(
                    row["segment"], row["chrom"], int(row["window_start"]), int(row["window_end"]),
                    int(row["period"]), row["motif"], int(row["reference_repeats"]))
            segment.alleles[int(row["allele"])] = (row["sequence"], row["change"])
    genotypes = collections.defaultdict(dict)  # segment -> sample -> (first, second)
    with open(os.path.join(run, "segment_genotypes.tsv")) as handle:
        handle.readline()
        for line in handle:
            name, sample, first, second = line.rstrip("\n").split("\t")
            genotypes[name][sample] = (int(first), int(second))
    return segments, genotypes


class Record:
    __slots__ = ("chrom", "start", "end", "alleles", "calls")

    def __init__(self, chrom, start, alleles, calls):
        self.chrom, self.start = chrom, start
        self.end = start + len(alleles[0])
        self.alleles, self.calls = alleles, calls


def genotype_quality(caller, fields):
    if caller in HIPSTR:
        value = fields.get("Q", ".")
        if value in (".", ""):
            return None
        posterior = float(value)
        return 99 if posterior >= 1.0 else min(99, round(-10 * math.log10(1.0 - posterior)))
    value = fields.get("GQ", ".")
    return None if value in (".", "") else int(value)


def read_calls(path, caller, samples):
    """Per contig, the records in position order, each with every sample's (allele indices or None, GQ)."""
    records = collections.defaultdict(list)
    columns = None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                columns = {name: index for index, name in enumerate(fields)}
                continue
            keys = fields[8].split(":")
            calls = {}
            for sample in samples:
                values = dict(zip(keys, fields[columns[sample]].split(":")))
                gt = values.get("GT", ".").replace("|", "/").split("/")
                indices = None if any(one in (".", "") for one in gt) or len(gt) != 2 \
                    else (int(gt[0]), int(gt[1]))
                calls[sample] = (indices, genotype_quality(caller, values))
            alleles = [fields[3]] + fields[4].split(",")
            records[fields[0]].append(Record(fields[0], int(fields[1]) - 1, alleles, calls))
    for chrom in records:
        records[chrom].sort(key=lambda record: record.start)
    return records


class CallIndex:
    """Records of one VCF, indexed for overlap queries by contig."""

    def __init__(self, records):
        import bisect
        self.bisect = bisect
        self.records = records
        self.starts = {chrom: [record.start for record in rows] for chrom, rows in records.items()}
        self.longest = {chrom: max((record.end - record.start for record in rows), default=0)
                        for chrom, rows in records.items()}

    def overlapping(self, chrom, start, end):
        rows = self.records.get(chrom, [])
        if not rows:
            return []
        first = self.bisect.bisect_left(self.starts[chrom], start - self.longest[chrom])
        last = self.bisect.bisect_left(self.starts[chrom], end)
        return [record for record in rows[first:last] if record.end > start and record.start < end]


def spell(reference, start, end, edits):
    """The sequence over [start, end) with `edits` — (record start, record end, allele) — applied,
    or None if two edits touch one base."""
    pieces, cursor = [], start
    for edit_start, edit_end, allele in sorted(edits):
        if edit_start < cursor:
            return None
        pieces.append(reference[cursor:edit_start])
        pieces.append(allele)
        cursor = edit_end
    pieces.append(reference[cursor:end])
    return "".join(pieces)


def called_pairs(records, sample, reference, start, end):
    """Every unordered pair of sequences the sample's calls over [start, end) can spell, trying both
    assignments of each heterozygous record; None in a pair where edits collide."""
    choices = []
    for record in records:
        indices, _ = record.calls[sample]
        first, second = indices
        options = [(first, second)] if first == second else [(first, second), (second, first)]
        choices.append((record, options))
    pairs = set()
    def walk(position, edits_one, edits_two):
        if position == len(choices):
            one = spell(reference, start, end, edits_one)
            two = spell(reference, start, end, edits_two)
            pairs.add(tuple(sorted((one or "\0", two or "\0"))))
            return
        record, options = choices[position]
        for first, second in options:
            def edit(index):
                allele = record.alleles[index]
                return [] if index == 0 or allele == "*" else [(record.start, record.end, allele)]
            walk(position + 1, edits_one + edit(first), edits_two + edit(second))
    walk(0, [], [])
    return pairs


def score_cell(caller, segment, sample, truth_pair, calls, contigs):
    reference = contigs[segment.chrom]
    records = calls.overlapping(segment.chrom, segment.window_start, segment.window_end)
    truth_variant = any(index != 0 for index in truth_pair)
    if not records:
        return ("no_call" if truth_variant else "ref_right"), None
    if any(record.calls[sample][0] is None for record in records):
        return "no_call", None
    start = min([segment.window_start] + [record.start for record in records])
    end = max([segment.window_end] + [record.end for record in records])
    truth = tuple(sorted(
        reference[start:segment.window_start] + segment.alleles[index][0] + reference[segment.window_end:end]
        for index in truth_pair))
    variant_records = [record for record in records if record.calls[sample][0] != (0, 0)]
    qualities = [record.calls[sample][1] for record in variant_records if record.calls[sample][1] is not None]
    quality = min(qualities) if qualities else None
    if not variant_records:
        return ("called_reference" if truth_variant else "ref_right"), None
    pairs = called_pairs(records, sample, reference, start, end)
    if truth in pairs:
        return "right", quality
    if not truth_variant:
        return "false_call", quality
    truth_lengths = sorted(len(one) for one in truth)
    if any("\0" not in pair and sorted(len(one) for one in pair) == truth_lengths for pair in pairs):
        return "length_right", quality
    return "wrong", quality


def write_windows_beds(run, segments, contigs):
    inside = os.path.join(run, "score_windows.bed")
    outside = os.path.join(run, "score_outside_windows.bed")
    by_chrom = collections.defaultdict(list)
    for segment in segments.values():
        by_chrom[segment.chrom].append((segment.window_start, segment.window_end))
    with open(inside, "w") as one, open(outside, "w") as other:
        for chrom in contigs:
            cursor = 0
            for start, end in sorted(by_chrom.get(chrom, [])):
                one.write(f"{chrom}\t{start}\t{end}\n")
                if start > cursor:
                    other.write(f"{chrom}\t{cursor}\t{start}\n")
                cursor = end
            if cursor < len(contigs[chrom]):
                other.write(f"{chrom}\t{cursor}\t{len(contigs[chrom])}\n")
    return inside, outside


def run_vcfeval(run, caller, path, samples, stretches):
    """vcfeval per sample and stretch; returns rows of counts."""
    work = os.path.join(run, "vcfeval")
    os.makedirs(work, exist_ok=True)
    sdf = os.path.join(work, "reference.sdf")
    rtg = ["java", "-Xmx4g", "-jar", RTG_JAR]
    if not os.path.isdir(sdf):
        subprocess.run(rtg + ["format", "-o", sdf, os.path.join(run, "reference.fa")],
                       check=True, stdout=subprocess.DEVNULL)
    truth = os.path.join(work, "truth.vcf.gz")
    if not os.path.exists(truth + ".tbi"):
        subprocess.run(["bcftools", "view", "-Oz", "-o", truth, os.path.join(run, "truth.vcf")], check=True)
        subprocess.run(["bcftools", "index", "-t", truth], check=True)
    if not path.endswith(".gz"):
        packed = os.path.join(work, caller + ".vcf.gz")
        subprocess.run(["bcftools", "view", "-Oz", "-o", packed, path], check=True)
        path = packed
    if not os.path.exists(path + ".tbi"):
        subprocess.run(["bcftools", "index", "-t", "-f", path], check=True)
    rows = []
    for sample in samples:
        for stretch, bed in stretches:
            out = os.path.join(work, f"{caller}.{sample}.{stretch}")
            if not os.path.exists(os.path.join(out, "summary.txt")):
                subprocess.run(["rm", "-rf", out])
                subprocess.run(rtg + ["vcfeval", "-b", truth, "-c", path, "-t", sdf, "-o", out,
                                      "--sample", sample, "--bed-regions", bed, "--all-records"],
                               check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            with open(os.path.join(out, "summary.txt")) as handle:
                lines = [line.split() for line in handle if line.strip() and not line.startswith("-")]
            header, values = lines[0], lines[-1]
            row = dict(zip(header, values[-len(header):]))
            rows.append((caller, sample, stretch, row.get("True-pos-baseline", "0"),
                         row.get("False-pos", "0"), row.get("False-neg", "0")))
    return rows


def main():
    run = sys.argv[1]
    contigs = read_fasta(os.path.join(run, "reference.fa"))
    segments, genotypes = read_truth(run)
    with open(os.path.join(run, "libraries.tsv")) as handle:
        handle.readline()
        samples = [line.split("\t")[0] for line in handle]

    # SCORE_SELF_CHECK=1 scores the truth VCF as a caller: every cell must come out right or
    # ref_right, or the scorer, not a caller, is wrong.
    callers = CALLERS + ([("truth", "truth.vcf")] if os.environ.get("SCORE_SELF_CHECK") == "1" else [])
    summary = collections.Counter()
    with open(os.path.join(run, "score_cells.tsv"), "w") as cells:
        cells.write("caller\tsegment\tsample\tperiod\treference_repeats\trepeat_band\ttruth\t"
                    "truth_changes\toutcome\tgenotype_quality\n")
        for caller, name in callers:
            path = os.path.join(run, name)
            if not os.path.exists(path):
                print(f"[score] no {name}; {caller} not scored", file=sys.stderr)
                continue
            calls = CallIndex(read_calls(path, caller, samples))
            for segment in segments.values():
                band = repeat_band(segment.period, segment.reference_repeats)
                for sample in samples:
                    truth_pair = genotypes[segment.name][sample]
                    changes = ",".join(sorted({segment.alleles[i][1] for i in truth_pair if i != 0})) or "."
                    outcome, quality = score_cell(caller, segment, sample, truth_pair, calls, contigs)
                    truth = "variant" if changes != "." else "reference"
                    cells.write(f"{caller}\t{segment.name}\t{sample}\t{segment.period}\t"
                                f"{segment.reference_repeats}\t{band}\t{truth}\t{changes}\t{outcome}\t"
                                f"{'.' if quality is None else quality}\n")
                    summary[(caller, segment.period, band, truth, outcome)] += 1
    with open(os.path.join(run, "score_summary.tsv"), "w") as out:
        out.write("caller\tperiod\trepeat_band\ttruth\toutcome\tcells\n")
        for key in sorted(summary):
            out.write("\t".join(map(str, key)) + f"\t{summary[key]}\n")

    if os.environ.get("SKIP_VCFEVAL") == "1":
        return
    inside, outside = write_windows_beds(run, segments, contigs)
    rows = []
    for caller, name in CALLERS:
        path = os.path.join(run, name)
        if not os.path.exists(path):
            continue
        stretches = [("segments", inside)] + ([] if caller in HIPSTR else [("outside", outside)])
        rows += run_vcfeval(run, caller, path, samples, stretches)
    with open(os.path.join(run, "score_vcfeval.tsv"), "w") as out:
        out.write("caller\tsample\tstretch\ttrue_positives\tfalse_positives\tfalse_negatives\n")
        for row in rows:
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
