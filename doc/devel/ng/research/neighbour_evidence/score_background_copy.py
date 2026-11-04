"""How well the background-copy probability catches wrong fills: ng_neighbour_pass_probe's output
scored against the founder simulator's truth, at SNP records matching a truth record.

usage: uv run python score_background_copy.py <truth.vcf.gz> <probe.tsv>

A fill is a sample-site calling wrote ./. and the pass called (GQ >= 10). For each limit on the
probability that the sample copies the background template there (no sampled haplotype to copy),
the fills a run would keep and how many of them are wrong.
"""
import gzip, sys

LIMITS = [None, 0.5, 0.3, 0.2, 0.1, 0.05]
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
missing, fills = 0, []
for line in open(sys.argv[2]):
    f = line.rstrip('\n').split('\t')
    if f[0] == 'contig':
        continue
    key = (f[0], int(f[1]), f[2], f[3])
    if key not in truth or f[5] != 'flat':
        continue
    missing += 1
    if f[7] == 'called':
        fills.append((float(f[12]), '/'.join(sorted(f[8].split('/'))) != truth[key][f[4]]))
print(f'missing before: {missing:,}')
print('| limit on background copy | fills kept | of them wrong |')
print('|---|---|---|')
for limit in LIMITS:
    kept = [wrong for background, wrong in fills if limit is None or background <= limit]
    share = f' ({100 * sum(kept) / len(kept):.1f}%)' if kept else ''
    print(f"| {'none' if limit is None else limit} | {len(kept):,} | {sum(kept)}{share} |")
