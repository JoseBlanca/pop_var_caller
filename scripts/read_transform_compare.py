"""Compare two pvc VCFs under one of the relations `read_transform_oracle.sh` checks.

    python3 scripts/read_transform_compare.py <relation> <base.vcf> <transformed.vcf> [sample]

Run it inside the dev container (or anywhere with Python 3), on two VCFs written with
`--min-site-quality 0 --hidden-dup-tag`, so that no site is missing because it fell below a cutoff
or was filtered away: a call that got worse shows up as a worse call, not as an absence.

Records are matched by contig, position and REF. Genotypes are compared as the allele sequences
they name, not as indices, because the ALT list can legitimately differ between the two runs: a
deeper run can admit an alternative a shallower one did not. A genotype is unphased here, so
`0/1` and `1/0` are the same.

The relations:

  doubled  every read fed twice. No record's QUAL may fall, and no genotype with GQ >= 20 in the
           base may change. A QUAL that falls is split by whether the evidence fell or the
           artifact penalties grew: the quality before the penalties is QUAL + ABPEN + SPPEN.
  empty    one more sample, with no reads. Every record, with the empty sample's column dropped,
           must be the same as in the base, field by field.
  strand   every read's strand flag flipped. Every genotype and every SPPEN must be the same.
  alone    `sample` called alone against the same sample in the cohort. A measurement: the
           differences are counted and characterised, and the exit status is always 0.

Exit status 1 when a pass/fail relation is broken, 0 otherwise.
"""

import sys
from collections import Counter, namedtuple

MIN_CONFIDENT_GQ = 20
EXAMPLES = 12

# One allele a sample carries, as `alone` reports it: the call, and where it was made.
CarriedAt = namedtuple("CarriedAt", "call qual af where is_str")


class Record:
    """One VCF data line, parsed only as far as the comparisons need."""

    def __init__(self, line, samples):
        fields = line.rstrip("\n").split("\t")
        self.line = line.rstrip("\n")
        self.chrom, self.pos, self.ref = fields[0], int(fields[1]), fields[3]
        self.alts = [] if fields[4] == "." else fields[4].split(",")
        self.qual = None if fields[5] == "." else float(fields[5])
        self.filter = fields[6]
        self.info = {}
        for item in fields[7].split(";"):
            key, _, value = item.partition("=")
            self.info[key] = value if value else True
        self.info_text = fields[7]
        keys = fields[8].split(":") if len(fields) > 8 else []
        self.calls = {}
        self.call_text = {}
        for name, value in zip(samples, fields[9:]):
            self.calls[name] = dict(zip(keys, value.split(":")))
            self.call_text[name] = value
        self.fixed = fields[:8]

    def key(self):
        return (self.chrom, self.pos, self.ref)

    def where(self):
        return f"{self.chrom}:{self.pos} {self.ref}>{','.join(self.alts) or '.'}"

    def is_str(self):
        return "STR" in self.info

    def float_info(self, name):
        value = self.info.get(name)
        return 0.0 if value in (None, True, ".") else float(value)

    def genotype(self, sample):
        """The called alleles as sequences, sorted; None for a missing call."""
        text = self.calls.get(sample, {}).get("GT", ".")
        indices = text.replace("|", "/").split("/")
        if any(index == "." for index in indices):
            return None
        alleles = [self.ref] + self.alts
        return tuple(sorted(alleles[int(index)] for index in indices))

    def gq(self, sample):
        value = self.calls.get(sample, {}).get("GQ", ".")
        return None if value == "." else int(value)

    def carries_alternative(self, sample):
        genotype = self.genotype(sample)
        return genotype is not None and any(allele != self.ref for allele in genotype)


class Records:
    """A VCF's data lines keyed by (contig, position, REF), **parsed only when asked for**.

    Holding every record parsed — a dictionary per sample per record — took 18.5 GB on two
    63-sample tomato VCFs, more than the dev container's 16 GB. The raw lines take about a
    twentieth of that, and each comparison parses a record, uses it and lets it go.
    """

    def __init__(self, samples, lines):
        self.samples, self.lines = samples, lines

    def __contains__(self, key):
        return key in self.lines

    def __len__(self):
        return len(self.lines)

    def get(self, key):
        line = self.lines.get(key)
        return None if line is None else Record(line, self.samples)

    def items(self):
        for key, line in self.lines.items():
            yield key, Record(line, self.samples)

    def values(self):
        for _, record in self.items():
            yield record


def read_vcf(path):
    samples, lines = [], {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            fields = line.split("\t", 4)
            key = (fields[0], int(fields[1]), fields[3])
            # Two records at one position with one REF would collide; keep both, told apart by
            # their order in the file, which is the same in both runs when nothing else differs.
            while key in lines:
                key = key + ("again",)
            lines[key] = line
    return samples, Records(samples, lines)


def show(label, examples):
    for example in examples[:EXAMPLES]:
        print(f"    {label}: {example}")
    if len(examples) > EXAMPLES:
        print(f"    ... and {len(examples) - EXAMPLES} more")


def compare_doubled(base, doubled):
    samples, base_records = base
    _, doubled_records = doubled
    missing = [r for k, r in base_records.items() if k not in doubled_records]
    new = [r for k, r in doubled_records.items() if k not in base_records]
    evidence_fell, penalties_grew, genotype_changed = [], [], []
    filter_changed = []
    for key, before in base_records.items():
        after = doubled_records.get(key)
        if after is None:
            continue
        if before.filter != after.filter:
            filter_changed.append(f"{before.where()} {before.filter} -> {after.filter}")
        if before.qual is not None and after.qual is not None and after.qual < before.qual:
            raw_before = before.qual + before.float_info("ABPEN") + before.float_info("SPPEN")
            raw_after = after.qual + after.float_info("ABPEN") + after.float_info("SPPEN")
            text = (
                f"{before.where()}{' STR' if before.is_str() else ''} QUAL {before.qual:g} -> "
                f"{after.qual:g}; before penalties {raw_before:g} -> {raw_after:g}; "
                f"ABPEN {before.float_info('ABPEN'):g} -> {after.float_info('ABPEN'):g}, "
                f"SPPEN {before.float_info('SPPEN'):g} -> {after.float_info('SPPEN'):g}"
            )
            (evidence_fell if raw_after < raw_before else penalties_grew).append(text)
        for sample in samples:
            gq = before.gq(sample)
            if gq is None or gq < MIN_CONFIDENT_GQ:
                continue
            if before.genotype(sample) != after.genotype(sample):
                genotype_changed.append(
                    f"{before.where()}{' STR' if before.is_str() else ''} {sample} "
                    f"{before.call_text[sample]} -> {after.call_text.get(sample)} "
                    f"(ALT now {','.join(after.alts) or '.'})"
                )
    print(f"  {len(base_records)} records in the base, {len(doubled_records)} doubled; "
          f"{len(missing)} gone, {len(new)} new, {len(filter_changed)} with FILTER changed")
    print(f"  QUAL fell at {len(evidence_fell) + len(penalties_grew)}: the quality before "
          f"penalties fell at {len(evidence_fell)}, only the penalties grew at "
          f"{len(penalties_grew)}")
    print(f"  confident genotypes (GQ >= {MIN_CONFIDENT_GQ}) changed: {len(genotype_changed)}")
    show("gone", [r.where() for r in missing])
    show("FILTER", filter_changed)
    show("evidence fell", evidence_fell)
    show("penalties grew", penalties_grew)
    show("genotype", genotype_changed)
    return bool(missing or evidence_fell or penalties_grew or genotype_changed)


def compare_empty(base, with_empty, empty_sample):
    samples, base_records = base
    _, empty_records = with_empty
    missing = [r.where() for k, r in base_records.items() if k not in empty_records]
    new = [r for k, r in empty_records.items() if k not in base_records]
    by_field = Counter()
    examples = []
    for key, before in base_records.items():
        after = empty_records.get(key)
        if after is None:
            continue
        differing = []
        for index, name in enumerate(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]):
            if before.fixed[index] != after.fixed[index]:
                differing.append(name)
        if before.info_text != after.info_text:
            for name in sorted(set(before.info) | set(after.info)):
                if before.info.get(name) != after.info.get(name):
                    differing.append(f"INFO/{name}")
        for sample in samples:
            if before.call_text[sample] != after.call_text.get(sample):
                differing.append(f"{sample}")
        if differing:
            by_field.update(differing)
            examples.append(f"{before.where()} [{', '.join(differing)}]\n        {before.line}\n"
                            f"        {after.line}")
    new_carrying = [r.where() for r in new
                    if any(r.carries_alternative(s) for s in samples)]
    print(f"  {len(base_records)} records in the base; {len(missing)} gone, {len(new)} new "
          f"({len(new_carrying)} where a non-empty sample carries an alternative), "
          f"{len(examples)} changed")
    for field, count in by_field.most_common():
        print(f"    {field}: {count}")
    show("gone", missing)
    show("new, carried", new_carrying)
    show("changed", examples)
    empty_called = sum(1 for r in empty_records.values() if r.genotype(empty_sample) is not None)
    print(f"  records where the empty sample has a called genotype: {empty_called}")
    return bool(missing or new_carrying or examples or empty_called)


def compare_strand(base, swapped):
    samples, base_records = base
    _, swapped_records = swapped
    missing = [r.where() for k, r in base_records.items() if k not in swapped_records]
    new = [r.where() for k, r in swapped_records.items() if k not in base_records]
    genotype_changed, sppen_changed, other = [], [], Counter()
    for key, before in base_records.items():
        after = swapped_records.get(key)
        if after is None:
            continue
        for sample in samples:
            if before.genotype(sample) != after.genotype(sample):
                genotype_changed.append(f"{before.where()} {sample} {before.call_text[sample]} -> "
                                        f"{after.call_text.get(sample)}")
        if before.info.get("SPPEN") != after.info.get("SPPEN"):
            sppen_changed.append(f"{before.where()} SPPEN {before.info.get('SPPEN')} -> "
                                 f"{after.info.get('SPPEN')}, QUAL {before.qual} -> {after.qual}")
        if before.line != after.line:
            other["any field"] += 1
    print(f"  {len(base_records)} records; {len(missing)} gone, {len(new)} new; "
          f"genotypes changed {len(genotype_changed)}, SPPEN changed {len(sppen_changed)}, "
          f"records differing in any field {other['any field']}")
    show("gone", missing)
    show("new", new)
    show("genotype", genotype_changed)
    show("SPPEN", sppen_changed)
    return bool(missing or new or genotype_changed or sppen_changed)


def gq_band(gq):
    if gq is None:
        return "no GQ"
    return "GQ < 20" if gq < MIN_CONFIDENT_GQ else "GQ >= 20"


def minimal_allele(pos, ref, alt):
    """An allele trimmed to its shortest spelling, so that one variant written over two record
    spans compares equal: common trailing bases dropped, then common leading ones."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return (pos, ref, alt)


def carried_alleles(records, sample, min_qual):
    """Every alternative allele `sample` is called carrying, in minimal spelling, with its copy
    count, the call's GQ, and the record — over records of at least `min_qual`."""
    carried = {}
    for record in records.values():
        if record.qual is None or record.qual < min_qual:
            continue
        genotype = record.genotype(sample)
        if genotype is None:
            continue
        for alt in set(genotype):
            if alt == record.ref:
                continue
            key = (record.chrom,) + minimal_allele(record.pos, record.ref, alt)
            # The few fields the report prints, not the record: a parsed record holds every
            # sample's call, and this keeps one per allele the sample carries.
            seen = CarriedAt(record.call_text[sample], record.qual, record.info.get("AF"),
                             record.where(), record.is_str())
            carried[key] = (genotype.count(alt), record.gq(sample), seen)
    return carried


def compare_alone(cohort, alone, sample):
    for min_qual, heading in ((0.0, "every record"), (1.0, "records of QUAL >= 1, the default "
                                                           "--min-site-quality")):
        in_cohort = carried_alleles(cohort[1], sample, min_qual)
        by_itself = carried_alleles(alone[1], sample, min_qual)
        classes = Counter()
        examples = {}
        gq_shift = []
        for key in set(in_cohort) | set(by_itself):
            a, b = in_cohort.get(key), by_itself.get(key)
            record = (a or b)[2]
            kind = "STR" if record.is_str else "generic"
            if a is None or b is None:
                copies, gq, _ = a or b
                side = "carried only alone" if a is None else "carried only in the cohort"
                label = f"{side}, {gq_band(gq)}"
            elif a[0] == b[0]:
                label = "same copies"
                if a[1] is not None and b[1] is not None:
                    gq_shift.append(b[1] - a[1])
            else:
                label = (f"different copies, cohort {gq_band(a[1])}, alone {gq_band(b[1])}")
            label = f"{kind}: {label}"
            classes[label] += 1
            describe = lambda side: (f"{side[2].call} QUAL {side[2].qual} "
                                     f"AF {side[2].af} at {side[2].where[:60]}"
                                     if side else "-")
            examples.setdefault(label, []).append(
                f"{key[0]}:{key[1]} {key[2]}>{key[3]}\n        cohort {describe(a)}"
                f"\n        alone  {describe(b)}")
        print(f"  {heading}: alternative alleles {sample} carries in either run")
        for label, count in sorted(classes.items()):
            print(f"    {label}: {count}")
        if gq_shift:
            gq_shift.sort()
            n = len(gq_shift)
            print(f"    GQ alone minus GQ in cohort, same copies: median {gq_shift[n // 2]}, "
                  f"10th percentile {gq_shift[n // 10]}, 90th {gq_shift[9 * n // 10]}, "
                  f"lower alone at {sum(1 for d in gq_shift if d < 0)} of {n}")
        for label in sorted(examples):
            if not label.endswith("same copies"):
                show(label, examples[label])
    return False


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(2)
    relation, base_path, other_path = sys.argv[1:4]
    base, other = read_vcf(base_path), read_vcf(other_path)
    if relation == "doubled":
        broken = compare_doubled(base, other)
    elif relation == "empty":
        broken = compare_empty(base, other, sys.argv[4])
    elif relation == "strand":
        broken = compare_strand(base, other)
    elif relation == "alone":
        broken = compare_alone(base, other, sys.argv[4])
    else:
        print(f"unknown relation {relation}")
        sys.exit(2)
    sys.exit(1 if broken else 0)


if __name__ == "__main__":
    main()
