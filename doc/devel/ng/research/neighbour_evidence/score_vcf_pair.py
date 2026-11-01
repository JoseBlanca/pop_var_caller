"""Score two pvc VCFs of one cohort — calling alone, and with the neighbour-evidence pass — against
the founder simulator's truth.vcf, at SNP records whose position and alleles match a truth record.
One unit is one sample at one record. Rows split by whether calling (the first VCF) gave the
sample a genotype there; a record the second VCF dropped counts its samples as 0/0 calls only
where the truth says so (the pass dropped it for being called 0/0 everywhere).

usage: uv run python score_vcf_pair.py truth.vcf calling.vcf.gz refined.vcf.gz
"""
import collections, gzip, sys


def read(path, snps_only=True):
    opener = gzip.open if path.endswith('.gz') else open
    out, samples = {}, None
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('##'):
                continue
            f = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                samples = f[9:]
                continue
            if snps_only and not (len(f[3]) == 1 and all(len(a) == 1 for a in f[4].split(','))):
                continue
            alleles = [f[3]] + f[4].split(',')
            gts = {}
            for s, col in zip(samples, f[9:]):
                gt = col.split(':')[0].replace('|', '/')
                gts[s] = None if '.' in gt else tuple(sorted(alleles[int(i)] for i in gt.split('/')))
            out[(f[0], int(f[1]), f[3])] = gts
    return out, samples


truth, samples = read(sys.argv[1])
calling, _ = read(sys.argv[2])
refined, _ = read(sys.argv[3])
count = collections.Counter()
for key, before in calling.items():
    if key not in truth:
        continue
    after = refined.get(key)
    for s in samples:
        t = truth[key][s]
        kind = 'nonref' if any(a != key[2] for a in t) else 'ref'
        row = 'called by calling' if before[s] is not None else 'no-call by calling'
        count[(row, kind, 'n')] += 1
        count[(row, kind, 'calling wrong')] += before[s] is not None and before[s] != t
        a = after[s] if after is not None else (key[2], key[2])
        count[(row, kind, 'pass no-call')] += a is None
        count[(row, kind, 'pass wrong')] += a is not None and a != t
print('calling\ttruth\tsample-sites\tcalling wrong\tpass wrong\tpass no-call')
for row in ('called by calling', 'no-call by calling'):
    for kind in ('nonref', 'ref'):
        g = lambda k: count[(row, kind, k)]
        print(f"{row}\t{kind}\t{g('n')}\t{g('calling wrong')}\t{g('pass wrong')}\t{g('pass no-call')}")
