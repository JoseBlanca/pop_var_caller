"""How well the neighbour check catches wrong fills: ng_neighbour_pass_probe's output scored against
the founder simulator's truth, at SNP records matching a truth record.

usage: uv run python score_neighbour_check.py <truth.vcf.gz> <probe.tsv>

For each sample-site the probe filled (calling wrote ./., the pass called it with GQ >= 10), the
sample's neighbour check in that window gives a disagreement rate, disagreed / checked. For each
limit on that rate, the fills a run would keep and how many of them are wrong; a sample whose check
counted nothing is kept only with no limit.
"""
import gzip, sys

LIMITS = [None, 0.20, 0.10, 0.05, 0.03, 0.02, 0.01]

truth, samples = {}, None
for line in gzip.open(sys.argv[1], 'rt'):
    if line.startswith('##'):
        continue
    f = line.rstrip('\n').split('\t')
    if line.startswith('#'):
        samples = f[9:]
        continue
    if len(f[3]) == 1 and len(f[4]) == 1:
        truth[(f[0], int(f[1]), f[3], f[4])] = {
            s: '/'.join(sorted(g.replace('|', '/').split('/'))) for s, g in zip(samples, f[9:])}

missing = 0
fills = []  # (rate or None, wrong)
for line in open(sys.argv[2]):
    f = line.rstrip('\n').split('\t')
    if f[0] == 'contig':
        continue
    key = (f[0], int(f[1]), f[2], f[3])
    if key not in truth or f[5] != 'flat':
        continue
    missing += 1
    if f[7] != 'called':
        continue
    checked, disagreed = int(f[10]), int(f[11])
    rate = disagreed / checked if checked else None
    fills.append((rate, '/'.join(sorted(f[8].split('/'))) != truth[key][f[4]]))

print(f'missing before: {missing:,}')
print('| limit on disagreement | fills kept | of them wrong |')
print('|---|---|---|')
for limit in LIMITS:
    kept = [wrong for rate, wrong in fills
            if limit is None or (rate is not None and rate <= limit)]
    wrong = sum(kept)
    share = f' ({100 * wrong / len(kept):.1f}%)' if kept else ''
    print(f"| {'none' if limit is None else limit} | {len(kept):,} | {wrong}{share} |")
