"""Count the genotypes a sample's own reads do not support, in one or more pvc VCFs.

    uv run --no-project python scripts/genotype_read_support_check.py <vcf[.gz]>...

Two counts, each split into biallelic and multi-allelic records and into GQ below and at or
above 20, over sample calls with at least two reads in AD:

- **names an alternative with no reads**: the genotype carries an alternative allele whose AD is
  zero while another alternative in the record has reads;
- **no called allele has a read**: every allele the genotype carries has AD zero, though the
  sample has reads (a 0/0 call on reads that all show an alternative counts here).

Neither is always wrong: at a few reads the cohort's frequencies can outweigh one sample's reads.
What should be rare is either at GQ >= 20. On 2026-10-31 the tomato panel showed 1,620 of the
second at GQ >= 20 in one run, every one from reads whose error was minted near one — the defect
read_likelihoods.md §3.3's amendment of that day fixed. Run it on every benchmark VCF; a count
that grows between two builds is worth a look before anything else.
"""

import gzip
import sys
from collections import Counter


def records(path):
    opener = gzip.open if path.endswith((".gz", ".bgz")) else open
    with opener(path, "rt") as handle:
        samples = []
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = fields[9:]
                continue
            yield fields, samples


def check(path):
    counts = Counter()
    for fields, samples in records(path):
        alternatives = fields[4].split(",")
        kind = "multi-allelic" if len(alternatives) > 1 else "biallelic"
        keys = fields[8].split(":")
        for value in fields[9:]:
            call = dict(zip(keys, value.split(":")))
            genotype, depths = call.get("GT", "."), call.get("AD", ".")
            if "." in genotype or depths == "." or "." in depths:
                continue
            reads = [int(depth) for depth in depths.split(",")]
            if sum(reads) < 2:
                continue
            called = {int(allele) for allele in genotype.replace("|", "/").split("/")}
            try:
                confident = float(call.get("GQ", "0")) >= 20
            except ValueError:
                confident = False
            band = "GQ>=20" if confident else "GQ<20"
            counts[(kind, "calls with two or more reads", "")] += 1
            alternative_without_reads = any(
                reads[allele] == 0
                and any(depth for other, depth in enumerate(reads) if other not in (0, allele))
                for allele in called - {0}
            )
            if alternative_without_reads:
                counts[(kind, "names an alternative with no reads", band)] += 1
            if all(reads[allele] == 0 for allele in called):
                counts[(kind, "no called allele has a read", band)] += 1
    return counts


def main():
    for path in sys.argv[1:]:
        print(path)
        for (kind, what, band), count in sorted(check(path).items()):
            print(f"  {kind}\t{what}\t{band}\t{count}")


if __name__ == "__main__":
    main()
