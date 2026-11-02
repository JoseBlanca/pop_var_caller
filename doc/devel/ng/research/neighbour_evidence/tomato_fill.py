"""The neighbour-evidence pass on a real cohort: how much missing data it fills, and how often a fill
is wrong, measured by hiding genotypes the reads settled.

usage:
  uv run python tomato_fill.py missing <calls.vcf[.gz]> [<refined.vcf[.gz]>]
  uv run python tomato_fill.py mask <calls.vcf[.gz]> <masked.vcf> <hidden.tsv> [--share 0.02]
                                    [--min-depth 6] [--min-gq 30] [--seed 1]
  uv run python tomato_fill.py score <hidden.tsv> <refined.vcf[.gz]>

missing  Sample-sites at PASS records written ./., before and after the pass; with a second VCF,
         by record kind (SNP, indel, STR) and both counts.
mask     Hides a share of the genotypes calling made from enough reads: biallelic PASS SNPs where
         the sample has at least --min-depth reads and GQ at least --min-gq. A hidden genotype is
         written as calling writes a sample with no reads — GT ./., GQ ., DP 0, AD and PL zeros —
         and its GP is kept, so other samples copying it see what its reads said, as they would
         had it been missing for want of reads only at this one record. hidden.tsv keeps the
         genotypes hidden, the truth the score is measured against.
score    For each hidden genotype, what the pass wrote: a genotype, right or wrong, or ./.; by the
         written GQ, and by whether the truth carries the alternative.

What the masked score measures is a fill at a record where the sample had reads — a site with
coverage, among neighbours of the same sample; a site with no reads in a low-coverage accession can
be harder (a repeat, a duplication). Read it as the pass's accuracy where the template evidence is
real, not as a forecast for every ./. in the cohort.
"""
import collections, gzip, random, sys


def open_text(path):
    return gzip.open(path, 'rt') if path.endswith(('.gz', '.bgz')) else open(path)


def records(path):
    """(fields, format keys, sample names) for every PASS record."""
    samples = None
    for line in open_text(path):
        if line.startswith('##'):
            continue
        f = line.rstrip('\n').split('\t')
        if line.startswith('#'):
            samples = f[9:]
            continue
        if f[6] != 'PASS':
            continue
        yield f, f[8].split(':'), samples


def kind(f):
    if 'STR' in f[7].split(';'):
        return 'STR'
    alleles = [f[3]] + f[4].split(',')
    return 'SNP' if all(len(a) == 1 for a in alleles) else 'indel'


def is_missing(gt):
    return '.' in gt


def missing(paths):
    for path in paths:
        counts = collections.Counter()
        for f, keys, samples in records(path):
            k = kind(f)
            for column in f[9:]:
                counts[k, 'all'] += 1
                counts[k, 'missing'] += is_missing(column.split(':', 1)[0])
        total = sum(counts[k, 'all'] for k in ('SNP', 'indel', 'STR'))
        gone = sum(counts[k, 'missing'] for k in ('SNP', 'indel', 'STR'))
        print(f'{path}: {gone:,} of {total:,} sample-genotypes missing ({100 * gone / max(total, 1):.1f}%)')
        for k in ('SNP', 'indel', 'STR'):
            n, m = counts[k, 'all'], counts[k, 'missing']
            if n:
                print(f'  {k}: {m:,} of {n:,} ({100 * m / n:.1f}%)')


def mask(source, out_path, hidden_path, share, min_depth, min_gq, seed):
    draw = random.Random(seed)
    hidden = 0
    with open(out_path, 'w') as out, open(hidden_path, 'w') as truth:
        truth.write('contig\tpos\tref\talt\tsample\tgt\tdepth\tgq\n')
        samples = None
        for line in open_text(source):
            if line.startswith('#'):
                out.write(line)
                if not line.startswith('##'):
                    samples = line.rstrip('\n').split('\t')[9:]
                continue
            f = line.rstrip('\n').split('\t')
            if f[6] != 'PASS' or len(f[3]) != 1 or len(f[4]) != 1:
                out.write(line)
                continue
            keys = f[8].split(':')
            at = {k: i for i, k in enumerate(keys)}
            for s, column in enumerate(f[9:]):
                values = column.split(':')
                if len(values) < len(keys) or is_missing(values[at['GT']]):
                    continue
                if int(values[at['DP']]) < min_depth or values[at['GQ']] == '.' \
                        or int(values[at['GQ']]) < min_gq or draw.random() >= share:
                    continue
                truth.write(f"{f[0]}\t{f[1]}\t{f[3]}\t{f[4]}\t{samples[s]}\t{values[at['GT']]}\t"
                            f"{values[at['DP']]}\t{values[at['GQ']]}\n")
                values[at['GT']] = './.'
                values[at['GQ']] = '.'
                values[at['DP']] = '0'
                values[at['AD']] = '0,0'
                values[at['PL']] = '0,0,0'
                f[9 + s] = ':'.join(values)
                hidden += 1
            out.write('\t'.join(f) + '\n')
    print(f'hid {hidden:,} genotypes')


def norm(gt):
    return '/'.join(sorted(gt.replace('|', '/').split('/')))


def score(hidden_path, refined):
    truth = {}
    for line in open(hidden_path):
        c = line.rstrip('\n').split('\t')
        if c[0] == 'contig':
            continue
        truth.setdefault((c[0], int(c[1]), c[2], c[3]), {})[c[4]] = norm(c[5])
    seen = collections.Counter()
    by_gq = collections.Counter()
    found = set()
    for f, keys, samples in records(refined):
        key = (f[0], int(f[1]), f[3], f[4])
        if key not in truth:
            continue
        found.add(key)
        at = {k: i for i, k in enumerate(keys)}
        for s, column in enumerate(f[9:]):
            want = truth[key].get(samples[s])
            if want is None:
                continue
            values = column.split(':')
            carried = 'non-ref' if want != '0/0' else 'ref'
            if is_missing(values[at['GT']]):
                seen['missing', carried] += 1
                continue
            right = norm(values[at['GT']]) == want
            seen['right' if right else 'wrong', carried] += 1
            gq = int(values[at['GQ']])
            band = next(b for b in ((10, 19), (20, 29), (30, 49), (50, 98), (99, 99), (0, 9))
                        if b[0] <= gq <= b[1])
            by_gq[band, 'n'] += 1
            by_gq[band, 'wrong'] += not right
    lost = sum(len(v) for k, v in truth.items() if k not in found)
    total = sum(len(v) for v in truth.values())
    filled = sum(seen['right', c] + seen['wrong', c] for c in ('ref', 'non-ref'))
    wrong = seen['wrong', 'ref'] + seen['wrong', 'non-ref']
    print(f'hidden genotypes: {total:,}')
    print(f'  filled: {filled:,} ({100 * filled / max(total, 1):.1f}%), wrong {wrong:,}'
          f' ({100 * wrong / max(filled, 1):.2f}% of fills)')
    print(f"  left ./.: {seen['missing', 'ref'] + seen['missing', 'non-ref']:,}; "
          f'record no longer written: {lost:,}')
    print('| truth | filled right | filled wrong | left ./. |')
    print('|---|---|---|---|')
    for c in ('ref', 'non-ref'):
        print(f"| {c} | {seen['right', c]:,} | {seen['wrong', c]:,} | {seen['missing', c]:,} |")
    print()
    print('| written GQ | fills | wrong |')
    print('|---|---|---|')
    for band in ((0, 9), (10, 19), (20, 29), (30, 49), (50, 98), (99, 99)):
        n, w = by_gq[band, 'n'], by_gq[band, 'wrong']
        print(f"| {band[0]}–{band[1]} | {n:,} | {w:,}" + (f' ({100 * w / n:.1f}%)' if n else '') + ' |')


def option(args, name, default, kind):
    if name in args:
        return kind(args[args.index(name) + 1])
    return default


if __name__ == '__main__':
    command, args = sys.argv[1], sys.argv[2:]
    if command == 'missing':
        missing(args)
    elif command == 'mask':
        mask(args[0], args[1], args[2], option(args, '--share', 0.02, float),
             option(args, '--min-depth', 6, int), option(args, '--min-gq', 30, int),
             option(args, '--seed', 1, int))
    elif command == 'score':
        score(args[0], args[1])
    else:
        sys.exit(__doc__)
