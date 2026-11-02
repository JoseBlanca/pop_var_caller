"""What perfect phasing would be worth to the neighbour-evidence pass: the ceiling of any phasing
trick, measured on founder_grid.py cells with the Python reference walk (hmm.haplotypes).

usage: uv run python phasing_ceiling.py <cell dir> [<cell dir> ...]

Per target sample, its 16 templates are chosen as pvc chooses them — the other samples that
disagree least where both have reads, by the allele shares of their written GP — and the target is
re-called three ways, the templates differing only in what they carry:

  pvc        each template sample's allele shares from its GP, as the pass reads them
  genotypes  each template sample's true genotype, unphased: shares 0, ½ or 1
  haplotypes each template sample's two true haplotypes as two templates, 32 in all

  subtraction  each template sample's two haplotypes as phased from the data, by subtraction (the
             owner's idea, 2026-11-02): one copy follows the sample's closest relatives along the
             chromosome — a one-copy copying walk over their haplotypes — and the other copy is
             what the sample's own genotype leaves once that copy is taken out, from its reads
             and, without reads, its inbreeding and the cohort's frequencies. Round one copies
             the relatives' GP; round two their round-one haplotypes. The target is left out of
             every step that phases its templates, so its reads do not come back to it.

  beagle     each template sample's two haplotypes as Beagle 5.5 phased them from calling's GT,
             the target left out of the run (beagle_phase.py); only where its runs exist

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
ARMS = ('pvc', 'genotypes', 'haplotypes', 'subtraction', 'beagle')
MISMATCH = 0.01


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


def ranking(records, target, n):
    """Every other sample, best template first, as pvc ranks them."""
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
    return [o for _, o in sorted(scores)]


def choose(records, target, n, exclude=()):
    return [o for o in RANKS[id(records)][target] if o not in exclude][:K]


def one_copy_phase(rows, sample, template_alt, inbreeding):
    """The subtraction for one sample over one contig's records `rows`. `template_alt` is a
    (sites, T) array of each template's probability of carrying the alternative. Copy one copies
    one template at a time (switching as one copy of the pass's walk does); copy two is drawn
    given copy one — the same allele with probability F, else from the cohort's frequency — and
    the sample's genotype likelihood ties the two. Returns (copy one's, copy two's) probability of
    carrying the alternative at each record."""
    n, T = template_alt.shape
    q = np.stack([1 - template_alt, template_alt], axis=2) * (1 - MISMATCH) + MISMATCH / 2
    f = min(max(inbreeding, 1e-4), 1 - 1e-4)
    E = np.empty((n, 2))
    second = np.empty(n)  # P(copy two alt | copy one alt-or-not) mixed below
    given = np.empty((n, 2))
    for i, r in enumerate(rows):
        l0, l1, l2 = r['lik'][sample]
        L = np.array([[l0, l1], [l1, l2]])
        pi = np.array([1 - r['af'], r['af']])
        for a in range(2):
            cond = (1 - f) * pi + f * np.eye(2)[a]
            joint = L[a] * cond
            E[i, a] = joint.sum()
            given[i, a] = joint[1] / joint.sum()
    e = np.einsum('tka,ta->tk', q, E)
    pos = np.array([r['pos'] for r in rows], dtype=float)
    fwd = np.empty((n, T))
    prev = None
    for i in range(n):
        if prev is None:
            cur = np.full(T, 1 / T)
        else:
            keep = math.exp(-(pos[i] - pos[i - 1]) / SEGMENT_BASES)
            cur = keep * prev + (1 - keep) / T
        cur = cur * e[i]
        cur /= cur.sum()
        fwd[i] = prev = cur
    bwd = np.ones((n, T))
    for i in range(n - 2, -1, -1):
        keep = math.exp(-(pos[i + 1] - pos[i]) / SEGMENT_BASES)
        b = bwd[i + 1] * e[i + 1]
        cur = keep * b + (1 - keep) * b.sum() / T
        bwd[i] = cur / cur.sum()
    post = fwd * bwd
    post /= post.sum(axis=1, keepdims=True)
    # copy one's allele given the template and the sample's own evidence, averaged over templates
    alt_given_k = q[:, :, 1] * E[:, 1][:, None] / np.einsum('tka,ta->tk', q, E)
    first = (post * alt_given_k).sum(axis=1)
    second = (1 - first) * given[:, 0] + first * given[:, 1]
    return first, second


def subtraction_haplotypes(cell, target, wanted):
    """Phased haplotypes, by contig, of the samples `wanted`, built without `target`."""
    samples, records, inbreeding = DATA[cell]
    n = len(samples)
    by_contig = CONTIGS[cell]
    rank = RANKS[id(records)]
    # round one: every sample but the target, copying its relatives' GP
    round_one = {}
    needed = set(wanted)
    for w in wanted:
        needed.update(o for o in rank[w] if o != target)
    for t in needed:
        if t == target:
            continue
        templates = [o for o in rank[t] if o != target][:K]
        out = {}
        for contig, rows in by_contig.items():
            alt = np.array([[r['shares'][o][1] for o in templates] for r in rows])
            out[contig] = one_copy_phase(rows, t, alt, inbreeding[t])
        round_one[t] = out
    # round two: the wanted samples, copying their relatives' round-one haplotypes
    result = {}
    for t in wanted:
        templates = [o for o in rank[t] if o != target][:K]
        out = {}
        for contig, rows in by_contig.items():
            alt = np.column_stack([round_one[o][contig][c] for o in templates for c in (0, 1)])
            out[contig] = one_copy_phase(rows, t, alt, inbreeding[t])
        result[t] = out
    return result


def beagle_haplotypes(cell, target):
    """Sample name -> {(contig, position): (allele, allele)} from the run without `target`, or None."""
    path = os.path.join(cell, 'beagle', f'out_{target}.vcf.gz')
    if not os.path.exists(path):
        return None
    out = {}
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            f = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                names = f[9:]
                out = {name: {} for name in names}
                continue
            for name, gt in zip(names, f[9:]):
                out[name][(f[0], int(f[1]))] = tuple(int(a) for a in gt.split('|'))
    return out


def template_rows(record, templates, arm, phased=None, index=None, beagle=None, samples=None):
    if arm == 'beagle':
        rows = [np.eye(2)[a] for t in templates
                for a in beagle[samples[t]][(record['contig'], record['pos'])]]
        rows.append(background_shares([record['af']]))
        return np.array(rows)
    if arm == 'subtraction':
        rows = [np.array([1 - phased[t][record['contig']][c][index], phased[t][record['contig']][c][index]])
                for t in templates for c in (0, 1)]
        rows.append(background_shares([record['af']]))
        return np.array(rows)
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
    by_contig = CONTIGS[cell]
    phased = subtraction_haplotypes(cell, target, templates)
    beagle = beagle_haplotypes(cell, target)
    counts = {}
    for arm in ARMS:
        if arm == 'beagle' and beagle is None:
            continue
        c = dict(reads=0, calling_wrong=0, pass_wrong=0, missing=0, filled=0, filled_wrong=0)
        for rows in by_contig.values():
            sites = [(r['lik'][target], None, None) for r in rows]
            posteriors, _ = haplotypes(sites, [r['pos'] for r in rows],
                                       [template_rows(r, templates, arm, phased, i, beagle, samples)
                                        for i, r in enumerate(rows)],
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
CONTIGS = {}
RANKS = {}


def main():
    cells = sys.argv[1:]
    for cell in cells:
        DATA[cell] = load(cell)
        samples, records, _ = DATA[cell]
        CONTIGS[cell] = {}
        for r in records:
            CONTIGS[cell].setdefault(r['contig'], []).append(r)
        RANKS[id(records)] = [ranking(records, t, len(samples)) for t in range(len(samples))]
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
            if not c:
                continue
            share = f" ({100 * c['filled_wrong'] / c['filled']:.1f}%)" if c['filled'] else ''
            print(f"| {name} | {arm} | {c['calling_wrong']:,} → {c['pass_wrong']:,} of {c['reads']:,} | "
                  f"{c['missing']:,} | {c['filled']:,} | {c['filled_wrong']:,}{share} |", flush=True)


if __name__ == '__main__':
    main()
