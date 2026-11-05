"""Scores the phase a pvc VCF writes against the founder simulator's truth, splitting the links by
where they came from (read phase spec §7.1).

usage: uv run python read_phase_score.py <cell dir> <with reads.vcf.gz> [<population only.vcf.gz>]

The VCF never says which evidence made a link, so the source is read off a second run of the same
cohort by a build whose pass uses no reads (the step-8 build of branch `neighbour-pass`): a link of
the first run is *population* when the second run phases the same two records in one block with
the same relation, and *made by reads* otherwise; a link of the second run whose two records the
first run does not phase together is *broken by reads*, and it is scored for whether it was right.
Without the second VCF every link is counted as made by reads, which is the case for a run without
refinement.

Scoring is phase_vcf_score.py's: SNP records matched to the truth by position, a phased column
*right* where the truth is heterozygous for the same two alleles and copy 0 of the truth carries the
allele written first; a link is two consecutive phased, truly heterozygous records of one sample in
one block, *wrong* where their rightness differs. Rows split deep and shallow samples in a
mixed-depth cell (the first `count` samples of a name ending `+<depth>x<count>` are deep).
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
    depth = name.split('_d', 1)[1].split('_s', 1)[0]
    return int(depth.split('x')[1]) if '+' in depth else 0


def truth_of(cell):
    truth = {}
    for kind, fields in records(os.path.join(cell, 'sim', 'truth.vcf.gz')):
        if kind == 'record':
            alleles = [fields[3]] + fields[4].split(',')
            truth[(fields[0], int(fields[1]))] = (
                alleles, [tuple(int(a) for a in g.split('|')) for g in fields[9:]])
    return truth


def phased_columns(vcf, truth):
    """Per sample, its phased truly heterozygous SNP columns in order, as
    (contig, pos, ps, right, copy-0 base), and the counts of true heterozygotes and of phased
    columns with a wrong genotype."""
    phased = defaultdict(list)
    true_hets = 0
    wrong_genotype = 0
    for kind, fields in records(vcf):
        if kind == 'header':
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
            true_hets += het
            gt = values[0]
            if '|' not in gt or ps_at is None:
                continue
            first, second = (alleles[int(a)] for a in gt.split('|'))
            if not het or sorted([first, second]) != sorted(true_pair):
                wrong_genotype += 1
                continue
            phased[sample].append((key[0], key[1], values[ps_at], first == true_pair[0], first))
    return phased, true_hets, wrong_genotype


def links_of(phased):
    """Per sample, its links: {(contig, pos, pos'): (wrong, relation)} — relation, whether the two
    records' copy 0 bases are both the truth's copy 0 or both not, read off the written orientation
    alone as 'same' or 'opposite' rightness."""
    out = {}
    for sample, rows in phased.items():
        for (c0, p0, b0, r0, _), (c1, p1, b1, r1, _) in zip(rows, rows[1:]):
            if c0 == c1 and b0 == b1:
                out[(sample, c0, p0, p1)] = r0 != r1
    return out


def blocks_of(phased):
    """Per sample and record, (block, rightness): two records are phased together with the same
    relation when they share a block and the rightness relation is the same."""
    out = {}
    for sample, rows in phased.items():
        for contig, pos, ps, right, _ in rows:
            out[(sample, contig, pos)] = ((contig, ps), right)
    return out


def main():
    cell, with_reads = sys.argv[1], sys.argv[2]
    population_only = sys.argv[3] if len(sys.argv) > 3 else None
    name = os.path.basename(os.path.normpath(cell))
    deep = deep_count(name)
    truth = truth_of(cell)
    phased, true_hets, wrong_genotype = phased_columns(with_reads, truth)
    links = links_of(phased)
    group = (lambda sample: 'deep' if sample < deep else 'shallow') if deep else (lambda _: 'all')
    count = Counter()
    if population_only:
        other, _, _ = phased_columns(population_only, truth)
        other_blocks = blocks_of(other)
        mine_blocks = blocks_of(phased)
        for (sample, contig, p0, p1), wrong in links.items():
            a, b = other_blocks.get((sample, contig, p0)), other_blocks.get((sample, contig, p1))
            # Same block in the population-only run, with the same relation between the two.
            together = a and b and a[0] == b[0] and (a[1] == b[1]) == (not wrong)
            source = 'population' if together else 'reads'
            count[source, group(sample), 'links'] += 1
            count[source, group(sample), 'wrong'] += wrong
        for (sample, contig, p0, p1), wrong in links_of(other).items():
            a, b = mine_blocks.get((sample, contig, p0)), mine_blocks.get((sample, contig, p1))
            if not (a and b and a[0] == b[0]):
                count['broken', group(sample), 'links'] += 1
                count['broken', group(sample), 'wrong'] += wrong
    else:
        for (sample, *_), wrong in links.items():
            count['reads', group(sample), 'links'] += 1
            count['reads', group(sample), 'wrong'] += wrong
    phased_total = sum(len(rows) for rows in phased.values())
    groups = sorted({key[1] for key in count})
    for source in ('population', 'reads', 'broken'):
        cells = []
        for g in groups:
            n, w = count[source, g, 'links'], count[source, g, 'wrong']
            if n:
                cells.append(f'{g} {w} of {n:,}' + (f' ({100 * w / n:.2f}%)' if n >= 300 else ' (under 300)'))
        if cells:
            print(f'| {name} | {source} | {"; ".join(cells)} |')
    print(f'| {name} | phased | {phased_total:,} of {true_hets:,} true heterozygotes, '
          f'{wrong_genotype:,} phased with a wrong genotype |')


if __name__ == '__main__':
    main()
