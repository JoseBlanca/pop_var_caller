"""Give HipSTR's records a QUAL, from its own genotype posterior.

HipSTR writes `.` in QUAL and puts its confidence in FORMAT `Q`, the posterior probability of the
called genotype. Both scorers read QUAL — the tract scorer drops records without one, and vcfeval
sweeps it for the best-F1 rows — so this writes QUAL = −10·log10(1 − Q), Phred-scaled and capped at
99, from the single sample's Q. Nothing else in the record changes.

    python3 hipstr_qual.py in.vcf.gz out.vcf   (then bgzip and index out.vcf)
"""

import gzip
import math
import sys


def phred(q):
    if q >= 1.0:
        return 99.0
    return min(99.0, -10.0 * math.log10(max(1e-12, 1.0 - q)))


def main(src, dst):
    with gzip.open(src, "rt") as fin, open(dst, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 9 and f[5] == ".":
                keys = f[8].split(":")
                values = f[9].split(":")
                if "Q" in keys and keys.index("Q") < len(values):
                    try:
                        f[5] = f"{phred(float(values[keys.index('Q')])):.2f}"
                    except ValueError:
                        pass
            fout.write("\t".join(f) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
