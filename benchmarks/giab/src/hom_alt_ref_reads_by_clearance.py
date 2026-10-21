"""At truth 1/1 indels, where do the REF-looking reads end? (alignment reference bias test)

Usage: uv run --with pysam python -I benchmarks/giab/src/hom_alt_ref_reads_by_clearance.py HOM_ALT_TSV ALIGNMENT REF OUT_TSV

For each site of indel_error_by_repeat.py hom_alt table, each read crossing the slide with >= 1
aligned bases each side is classed ALT / REF / other as before, and we record its "clearance":
the smaller of (slide start - read's first aligned base) and (read's last aligned base - slide
end), in bases. Output: one row per (indel length bin, class, clearance bin) with read counts.
"""
import collections
import csv
import sys

import pysam

hom, aln, ref, out = sys.argv[1:5]
bam = pysam.AlignmentFile(aln, reference_filename=ref)
acc = collections.Counter()
for r in csv.DictReader(open(hom), delimiter="\t"):
    if int(r["spanning"]) < 20 or int(r["show_alt"]) == 0:
        continue
    chrom, s, e = r["chrom"], int(r["span_start"]), int(r["span_end"])
    want = len(r["alt"]) - len(r["ref"])
    L = abs(want)
    lbin = "1" if L == 1 else "2-3" if L <= 3 else "4-6" if L <= 6 else "7+"
    for a in bam.fetch(chrom, max(0, s - 10), e + 10):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.is_duplicate or a.is_qcfail:
            continue
        if a.mapping_quality < 20:
            continue
        if a.reference_start > s - 1 or a.reference_end < e + 1:
            continue
        net = 0
        rpos = a.reference_start
        for op, ln in a.cigartuples:
            if op in (0, 7, 8):
                rpos += ln
            elif op == 1:
                if s - 1 <= rpos <= e:
                    net += ln
            elif op in (2, 3):
                if s - 1 <= rpos <= e:
                    net -= ln
                rpos += ln
        cls = "alt" if net == want else "ref" if net == 0 else "other"
        clearance = min(s - a.reference_start, a.reference_end - e)
        cbin = ("01-02" if clearance < 3 else "03-04" if clearance < 5 else "05-09" if clearance < 10 else "10-14" if clearance < 15 else "15-24" if clearance < 25
                else "25+")
        acc[(lbin, cbin, cls)] += 1
with open(out, "w") as fh:
    fh.write("indel_len\tclearance\tclass\treads\n")
    for k in sorted(acc):
        fh.write("\t".join(k) + f"\t{acc[k]}\n")
