"""Scores the phase a pvc VCF writes (GT with |, PS) against the founder simulator's truth
(neighbour_evidence_phase_output.md §7.2).

usage: uv run python phase_vcf_score.py <refined.vcf.gz> <cell dir>

Records are matched to the truth by contig and position, alleles by their sequence; only SNP
records — every allele one base, in the call and in the truth — are scored, since an indel or a
repeat can be spelled differently in the two, which a match by sequence counts as wrong. A phased
column is *right* where the truth is heterozygous for the same two alleles and copy 0 of the truth
carries the allele written first. A *flip* is a link between two consecutive phased, truly
heterozygous records of one sample in one block whose rightness differs. Reported: true
heterozygotes phased, phased columns whose genotype is wrong, flips inside blocks, block N50 in
records and in kb, and the flips by deep and shallow samples (a mixed-depth cell's first samples
are the deep ones, as founder_grid.py names them).
"""
import gzip
import os
import sys
from collections import Counter, defaultdict


def records(path):
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                yield 'header', fields[9:]
                continue
            yield 'record', fields


def deep_count(name):
    depth = name.split('_d', 1)[1]
    return int(depth.split('x')[1]) if '+' in depth else 0


def n50(sizes):
    sizes = sorted(sizes, reverse=True)
    total, acc = sum(sizes), 0
    for size in sizes:
        acc += size
        if acc >= total / 2:
            return size
    return 0


def main():
    vcf, cell = sys.argv[1], sys.argv[2]
    name = os.path.basename(os.path.normpath(cell))
    deep = deep_count(name)
    truth = {}
    for kind, fields in records(os.path.join(cell, 'sim', 'truth.vcf.gz')):
        if kind == 'record':
            alleles = [fields[3]] + fields[4].split(',')
            truth[(fields[0], int(fields[1]))] = (
                alleles, [tuple(int(a) for a in g.split('|')) for g in fields[9:]])
    true_hets = Counter()
    phased = defaultdict(list)  # sample -> [(contig, pos, right, ps)]
    wrong_genotype = 0
    unphased_true_het = 0
    for kind, fields in records(vcf):
        if kind == 'header':
            samples = fields
            continue
        key = (fields[0], int(fields[1]))
        if key not in truth:
            continue
        alleles = [fields[3]] + fields[4].split(',')
        true_alleles, genotypes = truth[key]
        if any(len(allele) != 1 for allele in alleles + true_alleles):
            continue
        keys = fields[8].split(':')
        ps_at = keys.index('PS') if 'PS' in keys else None
        for sample, column in enumerate(fields[9:]):
            values = column.split(':')
            true_pair = [true_alleles[a] for a in genotypes[sample]]
            het = true_pair[0] != true_pair[1]
            true_hets[sample] += het
            gt = values[0]
            if '|' not in gt:
                unphased_true_het += het
                continue
            first, second = (alleles[int(a)] for a in gt.split('|'))
            if not het or sorted([first, second]) != sorted(true_pair):
                wrong_genotype += 1
                continue
            phased[sample].append((key[0], key[1], first == true_pair[0], int(values[ps_at])))
    links = flips = 0
    by_group = Counter()
    by_group_flips = Counter()
    sizes, spans = [], []
    for sample, rows in phased.items():
        blocks = defaultdict(list)
        for contig, pos, right, ps in rows:
            blocks[(contig, ps)].append(pos)
        for positions in blocks.values():
            sizes.append(len(positions))
            spans.append(max(positions) - min(positions))
        group = 'deep' if sample < deep else ('shallow' if deep else 'all')
        for (c0, p0, r0, b0), (c1, p1, r1, b1) in zip(rows, rows[1:]):
            if c0 == c1 and b0 == b1:
                links += 1
                flips += r0 != r1
                by_group[group] += 1
                by_group_flips[group] += r0 != r1
    het_total = sum(true_hets.values())
    phased_total = sum(len(rows) for rows in phased.values())
    groups = '; '.join(f'{g} {by_group_flips[g]} of {by_group[g]}' for g in sorted(by_group))
    print(f'| {name} | {phased_total:,} of {het_total:,} ({100 * phased_total / max(het_total, 1):.0f}%) '
          f'| {wrong_genotype:,} | {flips:,} of {links:,} ({100 * flips / max(links, 1):.2f}%) '
          f'| {n50(sizes)} ({n50(spans) / 1000:.0f}) | {groups} |')


if __name__ == '__main__':
    main()
