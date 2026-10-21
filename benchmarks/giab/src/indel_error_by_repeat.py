"""Measure how often a read shows an indel the sample does not carry, by repeat context.

Usage: uv run --with pysam --with numpy python -I benchmarks/giab/src/indel_error_by_repeat.py SAMPLE ALIGNMENT BED TRUTH_VCF REF OUT_TSV

Ground: every reference "stretch" inside the sample's confident regions with no truth variant
within 15 bp. A stretch is the longest exact tandem repeat (period 1..6) covering a position, or a
single base where there is none. For each stretch we count reads that cross it with at least 5
aligned bases on each side and no clip, and among them the reads carrying an insertion or
deletion inside the stretch (or touching its edges). Those reads are errors by construction: the
truth says the sample carries the reference there.

For each spanning read we also record what pvc would charge it if it were a REF read at an indel
locus covering this stretch: the lowest base quality over the stretch plus one base each side,
floored by the mapping-quality error (`minted_ln_read_error`), i.e. the generic path's per-read
error before calibration.

Second table (OUT_TSV.hom_alt.tsv): at truth 1/1 indels, the share of spanning reads that show the
reference instead, which is the ALT->REF error rate a 1/1 genotype has to explain.
"""

import collections
import gzip
import math
import sys

import numpy as np
import pysam

FLANK = 5
CLEAR = 15
MIN_MAPQ = 20  # pvc admits reads at mapping quality 20 and above (read filtering default)


def load_bed(path):
    out = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 3 and not line.startswith(("#", "track")):
                out[p[0]].append((int(p[1]), int(p[2])))
    return out


def load_truth(path):
    """All truth variant spans (0-based half open on the reference), and the 1/1 indels."""
    spans = collections.defaultdict(list)
    homalt = []
    with pysam.VariantFile(path) as vf:
        for rec in vf:
            start = rec.pos - 1
            end = start + len(rec.ref)
            spans[rec.chrom].append((start, end))
            gt = rec.samples[0]["GT"]
            alts = rec.alts or ()
            if len(alts) == 1 and gt in ((1, 1),) and len(alts[0]) != len(rec.ref):
                homalt.append((rec.chrom, rec.pos, rec.ref, alts[0]))
    for c in spans:
        spans[c].sort()
    return spans, homalt


def stretches(seq, offset):
    """Yield (start, end, period, unit, copies) over seq, each base in exactly one stretch.

    Greedy left to right: at each position take the period whose exact tandem run starting here
    is longest in bases (at least two full copies; period 1 needs a run of 2). Ties go to the
    smaller period. A base with no repeat starting at it is a period-0 stretch of length 1.
    """
    n = len(seq)
    i = 0
    while i < n:
        best = (1, 0, 1)  # (length, period, copies)
        for p in range(1, 7):
            if i + 2 * p > n:
                break
            j = i + p
            while j < n and seq[j] == seq[j - p]:
                j += 1
            length = j - i
            copies = length // p
            if copies >= 2 and length > best[0]:
                best = (length, p, copies)
        length, p, copies = best
        if p == 0:
            yield offset + i, offset + i + 1, 0, "", 1
            i += 1
        else:
            yield offset + i, offset + i + length, p, seq[i:i + p], copies
            i += length


def context_label(period, copies, length):
    if period == 0:
        return "none"
    if period == 1:
        return f"homopolymer_{min(length, 12)}"
    return f"period{period}_copies{min(copies, 8)}"


def ln_err(bq, mapq):
    e_b = -bq / 10 * math.log(10)
    e_m = -mapq / 10 * math.log(10)
    return max(e_b, e_m)


def read_indels(read):
    """(ref_pos_0based_before_event, kind, length) of every I/D in the read, plus clip flags."""
    out = []
    rpos = read.reference_start
    for op, ln in read.cigartuples:
        if op in (0, 7, 8):
            rpos += ln
        elif op == 1:
            out.append((rpos, "I", ln))
        elif op in (2, 3):
            out.append((rpos, "D", ln))
            rpos += ln
    return out


def spans_cleanly(read, start, end):
    """Aligned bases cover [start-FLANK, end+FLANK) with no clip in between (clips only at ends)."""
    if read.reference_start > start - FLANK or read.reference_end < end + FLANK:
        return False
    return True


def main():
    sample, aln, bed, truth, ref, out = sys.argv[1:7]
    regions = load_bed(bed)
    spans, homalt = load_truth(truth)
    fasta = pysam.FastaFile(ref)
    bam = pysam.AlignmentFile(aln, reference_filename=ref)

    # per context: stretches, spanning reads, indel reads (by kind: whole-unit / other),
    # and the sum of ln(charged) over spanning reads and count, for the geometric mean.
    acc = collections.defaultdict(lambda: [0, 0, 0, 0, 0.0, 0])
    for chrom, ivs in regions.items():
        tspans = spans.get(chrom, [])
        for (rs, re_) in ivs:
            seq = fasta.fetch(chrom, rs, re_).upper()
            sts = [s for s in stretches(seq, rs)]
            # keep stretches clear of truth variants and of the region edge
            keep = []
            ti = 0
            for (s, e, p, u, c) in sts:
                if s - CLEAR < rs or e + CLEAR > re_:
                    continue
                while ti < len(tspans) and tspans[ti][1] < s - CLEAR:
                    ti += 1
                clash = False
                k = ti
                while k < len(tspans) and tspans[k][0] <= e + CLEAR:
                    if tspans[k][1] >= s - CLEAR:
                        clash = True
                        break
                    k += 1
                if not clash:
                    keep.append((s, e, p, u, c))
            if not keep:
                continue
            # reads over the region once; for each read, test stretches it spans
            reads = []
            for r in bam.fetch(chrom, rs, re_):
                if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
                    continue
                if r.is_qcfail or r.mapping_quality < MIN_MAPQ:
                    continue
                reads.append(r)
            reads.sort(key=lambda r: r.reference_start)
            if not reads:
                continue
            starts = np.array([r.reference_start for r in reads])
            ends = np.array([r.reference_end for r in reads])
            # every indel event of every read: anchor position, read index, length
            ev_pos, ev_read, ev_len = [], [], []
            for idx, r in enumerate(reads):
                for (pos, kind, ln) in read_indels(r):
                    ev_pos.append(pos)
                    ev_read.append(idx)
                    ev_len.append(ln)
            order = np.argsort(ev_pos, kind="stable")
            ev_pos = np.array(ev_pos, dtype=np.int64)[order]
            ev_read = np.array(ev_read, dtype=np.int64)[order]
            ev_len = np.array(ev_len, dtype=np.int64)[order]
            qual_maps = {}
            for n_st, (s, e, p, u, c) in enumerate(keep):
                label = context_label(p, c, e - s)
                a = acc[label]
                a[0] += 1
                lo = np.searchsorted(starts, s - 1000, "left")
                hi = np.searchsorted(starts, s - FLANK, "right")
                cand = np.nonzero(ends[lo:hi] >= e + FLANK)[0] + lo
                a[1] += len(cand)
                if len(cand) == 0:
                    continue
                i0 = np.searchsorted(ev_pos, s - 1, "left")
                i1 = np.searchsorted(ev_pos, e, "right")
                if i1 > i0:
                    spanning = set(cand.tolist())
                    seen = set()
                    for k in range(i0, i1):
                        ri = int(ev_read[k])
                        if ri in spanning and ri not in seen:
                            seen.add(ri)
                            if p > 0 and ev_len[k] % p == 0:
                                a[2] += 1
                            else:
                                a[3] += 1
                # charged error, sampled: every 7th stretch, up to 25 reads
                if n_st % 7:
                    continue
                for ri in cand[:25]:
                    ri = int(ri)
                    m = qual_maps.get(ri)
                    if m is None:
                        r = reads[ri]
                        q = r.query_qualities
                        m = {rp: q[qp] for qp, rp in r.get_aligned_pairs(matches_only=True)}
                        qual_maps[ri] = m
                    mn = min((m.get(x, 99) for x in range(s - 1, e + 1)), default=99)
                    if mn == 99:
                        continue
                    a[4] += ln_err(mn, reads[ri].mapping_quality)
                    a[5] += 1

    with open(out, "w") as fh:
        fh.write("sample\tcontext\tstretches\tspanning_reads\tindel_reads_whole_unit\t"
                 "indel_reads_other\tindel_rate\tcharged_geomean\tcharged_reads\n")
        for label in sorted(acc):
            st, sp, w, o, lsum, n = acc[label]
            rate = (w + o) / sp if sp else float("nan")
            gm = math.exp(lsum / n) if n else float("nan")
            fh.write(f"{sample}\t{label}\t{st}\t{sp}\t{w}\t{o}\t{rate:.6g}\t{gm:.4g}\t{n}\n")

    # hom-alt indels: share of spanning reads showing the reference allele
    with open(out + ".hom_alt.tsv", "w") as fh:
        fh.write("sample\tchrom\tpos\tref\talt\tspan_start\tspan_end\tcontext\t"
                 "spanning\tshow_alt\tshow_ref\tshow_other\n")
        for (chrom, pos, refa, alta) in homalt:
            if not any(s <= pos - 1 < e for (s, e) in regions.get(chrom, [])):
                continue
            # the slide: extend over the repeat the indel's unit sits in
            if len(alta) > len(refa):
                unit = alta[len(refa):] if alta.startswith(refa) else alta[1:]
            else:
                unit = refa[len(alta):] if refa.startswith(alta) else refa[1:]
            L = len(unit)
            s0 = pos  # 0-based first base after anchor
            win = fasta.fetch(chrom, max(0, s0 - 200), s0 + 400).upper()
            off = max(0, s0 - 200)
            # slide right while the reference repeats the unit's rotation
            e = s0
            seqlen = len(win)
            # right slide for insertion: ref[e] == unit[(e-s0) % L]
            k = 0
            while (e - off) < seqlen and win[e - off] == unit[k % L]:
                e += 1
                k += 1
            if len(alta) < len(refa):
                e = max(e, s0 + L)
            # left slide
            s = s0
            k = 0
            while s - 1 - off >= 0 and win[s - 1 - off] == unit[(L - 1 - k) % L]:
                s -= 1
                k += 1
            span = e - s
            per = L
            ctx = context_label(per, span // per if per else 1, span)
            spanning = alt_n = ref_n = other = 0
            want = len(alta) - len(refa)
            for r in bam.fetch(chrom, max(0, s - 10), e + 10):
                if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
                    continue
                if r.is_qcfail or r.mapping_quality < MIN_MAPQ:
                    continue
                if not spans_cleanly(r, s, e):
                    continue
                spanning += 1
                net = 0
                for (p, kind, ln) in read_indels(r):
                    if s - 1 <= p <= e:
                        net += ln if kind == "I" else -ln
                if net == want:
                    alt_n += 1
                elif net == 0:
                    ref_n += 1
                else:
                    other += 1
            fh.write(f"{sample}\t{chrom}\t{pos}\t{refa}\t{alta}\t{s}\t{e}\t{ctx}\t"
                     f"{spanning}\t{alt_n}\t{ref_n}\t{other}\n")


if __name__ == "__main__":
    main()
