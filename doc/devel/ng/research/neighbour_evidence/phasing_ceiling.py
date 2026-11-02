"""What perfect phasing would be worth to the neighbour-evidence pass: the ceiling of any phasing
trick, measured on founder_grid.py cells with the Python reference walk (hmm.haplotypes).

usage: uv run python phasing_ceiling.py <cell dir> [<cell dir> ...]

Per target sample, its 16 templates are chosen as pvc chooses them — the other samples that
disagree least where both have reads, by the allele shares of their written GP — and the target is
re-called three ways, the templates differing only in what they carry:

  pvc        each template sample's allele shares from its GP, as the pass reads them
  genotypes  each template sample's true genotype, unphased: shares 0, ½ or 1
  haplotypes each template sample's two true haplotypes as two templates, 32 in all

plus the background template (the cohort's allele frequencies) in each. `genotypes` against
`haplotypes` is what phasing alone is worth; `pvc` against `genotypes` what uncertain genotypes
cost. Biallelic SNP PASS records matching a truth record only; segment and run lengths fixed at the
simulator's (300 kb, 2 Mb), F from the parameters file — no fitting, no background limit, no
evidence-nearby rule, so the numbers compare the three arms and are not the pass's own. A sample at
a record with no reads is filled when the posterior's GQ is at least 10.
"""
import gzip, math, multiprocessing, os, sys, tomllib
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hmm import background_shares, haplotypes  # noqa: E402

K = 16
COMPARED = 500
SEGMENT_BASES = 300_000.0
RUN_BASES = 2_000_000.0
ARMS = ('pvc', 'genotypes', 'haplotypes')


def load(cell):
    truth = {}
    with gzip.open(os.path.join(cell, 'sim', 'truth.vcf.gz'), 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f[3]) == 1 and len(f[4]) == 1:
                truth[(f[0], int(f[1]), f[3], f[4])] = [tuple(int(a) for a in g.split('|'))
                                                       for g in f[9:]]
    records = []
    with gzip.open(os.path.join(cell, 'calls_off.vcf.gz'), 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            f = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                samples = f[9:]
                continue
            key = (f[0], int(f[1]), f[3], f[4])
            if f[6] != 'PASS' or key not in truth:
                continue
            keys = f[8].split(':')
            pl_at, gp_at = keys.index('PL'), keys.index('GP')
            af = float(dict(kv.split('=', 1) for kv in f[7].split(';') if '=' in kv)['AF'])
            lik, shares, called = [], [], []
            for column in f[9:]:
                v = column.split(':')
                pl = [int(x) for x in v[pl_at].split(',')]
                lik.append(10 ** (-(np.array(pl) - min(pl)) / 10))
                gp = [float(x) for x in v[gp_at].split(',')]
                alt = gp[1] / 2 + gp[2]
                shares.append(np.array([1 - alt, alt]))
                called.append(None if '.' in v[0] else tuple(sorted(int(a) for a in v[0].split('/'))))
            records.append(dict(contig=f[0], pos=int(f[1]), af=af, lik=lik, shares=shares,
                                called=called, truth=truth[key]))
    params = tomllib.load(open(os.path.join(cell, 'params.toml'), 'rb'))
    inbreeding = {row['sample']: row['inbreeding_coefficient']['value']
                  for row in params['inbreeding']['by_sample']}
    return samples, records, [inbreeding[s] for s in samples]


def informative(lik):
    return lik.max() - lik.min() > 1e-9


def choose(records, target, n):
    rows = [r for r in records if informative(r['lik'][target])]
    if len(rows) > COMPARED:
        rows = [rows[i * len(rows) // COMPARED] for i in range(COMPARED)]
    scores = []
    for other in range(n):
        if other == target:
            continue
        shared = [r for r in rows if informative(r['lik'][other])]
        d = sum(0.5 * np.abs(r['shares'][target] - r['shares'][other]).sum() for r in shared)
        mean = (d + 1) / (len(shared) + 2)
        scores.append((mean + 2 * math.sqrt(mean * (1 - mean) / (len(shared) + 3)), other))
    return [o for _, o in sorted(scores)[:K]]


def template_rows(record, templates, arm):
    if arm == 'pvc':
        rows = [record['shares'][t] for t in templates]
    elif arm == 'genotypes':
        rows = [np.array([1 - sum(record['truth'][t]) / 2, sum(record['truth'][t]) / 2])
                for t in templates]
    else:
        rows = [np.eye(2)[a] for t in templates for a in record['truth'][t]]
    rows.append(background_shares([record['af']]))
    return np.array(rows)


def solve(args):
    cell, target = args
    samples, records, inbreeding = DATA[cell]
    templates = choose(records, target, len(samples))
    by_contig = {}
    for r in records:
        by_contig.setdefault(r['contig'], []).append(r)
    counts = {}
    for arm in ARMS:
        c = dict(reads=0, calling_wrong=0, pass_wrong=0, missing=0, filled=0, filled_wrong=0)
        for rows in by_contig.values():
            sites = [(r['lik'][target], None, None) for r in rows]
            posteriors, _ = haplotypes(sites, [r['pos'] for r in rows],
                                       [template_rows(r, templates, arm) for r in rows],
                                       inbreeding[target], RUN_BASES, SEGMENT_BASES)
            for r, post in zip(rows, posteriors):
                truth = tuple(sorted(r['truth'][target]))
                best = int(np.argmax(post))
                genotype = [(0, 0), (0, 1), (1, 1)][best]
                if r['called'][target] is not None:
                    c['reads'] += 1
                    c['calling_wrong'] += r['called'][target] != truth
                    c['pass_wrong'] += genotype != truth
                else:
                    c['missing'] += 1
                    gq = -10 * math.log10(max(1 - post[best], 1e-10))
                    if gq >= 10:
                        c['filled'] += 1
                        c['filled_wrong'] += genotype != truth
        counts[arm] = c
    return cell, counts


DATA = {}


def main():
    cells = sys.argv[1:]
    for cell in cells:
        DATA[cell] = load(cell)
    jobs = [(cell, t) for cell in cells for t in range(len(DATA[cell][0]))]
    totals = {cell: {arm: {} for arm in ARMS} for cell in cells}
    with ProcessPoolExecutor(6, mp_context=multiprocessing.get_context('fork')) as pool:
        for cell, counts in pool.map(solve, jobs):
            for arm, c in counts.items():
                for k, v in c.items():
                    totals[cell][arm][k] = totals[cell][arm].get(k, 0) + v
    print('| cell | templates | wrong where samples have reads, calling → pass | missing | filled (GQ ≥ 10) | of them wrong |')
    print('|---|---|---|---|---|---|')
    for cell in cells:
        name = os.path.basename(os.path.normpath(cell))
        for arm in ARMS:
            c = totals[cell][arm]
            share = f" ({100 * c['filled_wrong'] / c['filled']:.1f}%)" if c['filled'] else ''
            print(f"| {name} | {arm} | {c['calling_wrong']:,} → {c['pass_wrong']:,} of {c['reads']:,} | "
                  f"{c['missing']:,} | {c['filled']:,} | {c['filled_wrong']:,}{share} |", flush=True)


if __name__ == '__main__':
    main()
