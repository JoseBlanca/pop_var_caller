"""The phased pass on one small window, computed in Python for the Rust test that must match it to
1e-9 (phasing spec §8.1): `src/neighbour_evidence/testdata/phased_window.txt`.

usage: uv run --with numpy python phaser_reference.py > ../../../../../src/neighbour_evidence/testdata/phased_window.txt

The window is drawn here — samples copying a few founder haplotypes, partly autozygous, reads at
about 2× — and written as each record's PL and GP, then solved as `refine_window_in_rounds` solves
it: round 1 with the samples' GP shares as templates, later rounds with the haplotypes of the round
before (one of each sample, scored as `choose_haplotypes` scores them), a path drawn from the walk
with the same SplitMix64 stream and state order as `random.rs`, the haplotypes re-estimated along it,
twenty core records of each sample hidden, and the result the posterior averaged over rounds 2 to 5
(the prediction times the reads at hidden records). The walk is `hmm.haplotypes`, which the Rust
walk already matches to 1e-9. The parameters are fixed — the Rust side freezes each sample's fit at
the same grid point — and every input goes through the f32 the Rust records store.
"""
import math, os, sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hmm import haplotypes, stays  # noqa: E402

SAMPLES, RECORDS, FOUNDERS = 8, 300, 3
SPACING = 500
SEGMENT, RUN, SHARE, MISMATCH = 300_000.0, 300_000.0, 0.5, 0.01
CHOSEN = 1 * 6 + 3            # SEGMENT_GRID[1], RUN_GRID[3], segment-major
ROUNDS, WINDOW_KEY, K, HIDDEN = 5, 11, 16, 20
MASK = (1 << 64) - 1
SEED = 0x6E65_6967_6862_6F75


def finalise(z):
    z = ((z ^ (z >> 30)) * 0xBF58_476D_1CE4_E5B9) & MASK
    z = ((z ^ (z >> 27)) * 0x94D0_49BB_1331_11EB) & MASK
    return z ^ (z >> 31)


class Uniforms:
    def __init__(self, purpose, window, sample, round_):
        state = SEED
        for value in (purpose, window, sample, round_):
            state = finalise(state ^ value)
        self.state = state

    def next(self):
        self.state = (self.state + 0x9E37_79B9_7F4A_7C15) & MASK
        return (finalise(self.state) >> 11) / float(1 << 53)


def pick(weights, u):
    total = sum(weights)
    target = u * total
    cumulative, last_positive = 0.0, 0
    for index, weight in enumerate(weights):
        if weight > 0:
            last_positive = index
        cumulative += weight
        if cumulative > target:
            return index
    return last_positive


def f32(x):
    return float(np.float32(x))


# ---- the window ----------------------------------------------------------------------------------
draw = np.random.default_rng(20261102)
founders = draw.random((FOUNDERS, RECORDS)) < 0.4
afs = founders.mean(axis=0) * 0.8 + 0.1
pl, gp, present = [], [], []
for s in range(SAMPLES):
    copies = [draw.integers(FOUNDERS), draw.integers(FOUNDERS)]
    autozygous = draw.random() < 0.5
    pls, gps, here = [], [], []
    for i in range(RECORDS):
        for c in (0, 1):
            if draw.random() < 0.01:
                copies[c] = draw.integers(FOUNDERS)
        if draw.random() < 0.002:
            autozygous = not autozygous
        alleles = [int(founders[copies[0], i]), int(founders[copies[0] if autozygous else copies[1], i])]
        depth = draw.poisson(2.0)
        alt = sum(int(draw.random() < 0.5 if alleles[0] != alleles[1] else alleles[0] == 1)
                  for _ in range(depth))
        ref = depth - alt
        ll = [ref * math.log10(0.99) + alt * math.log10(0.01),
              depth * math.log10(0.5),
              ref * math.log10(0.01) + alt * math.log10(0.99)]
        top = max(ll)
        pls.append([int(round(-10 * (x - top))) for x in ll])
        p = afs[i]
        prior = [(1 - p) ** 2, 2 * p * (1 - p), p * p]
        post = [pr * 10 ** (x - top) for pr, x in zip(prior, ll)]
        tot = sum(post)
        gps.append([round(x / tot, 3) for x in post])
        here.append(draw.random() > 0.02)
    pl.append(pls); gp.append(gps); present.append(here)
positions = [1 + i * SPACING for i in range(RECORDS)]

# ---- as the Rust records store it ----------------------------------------------------------------
def row(s, i):
    if not present[s][i]:
        return 'absent'
    return 'flat' if len(set(pl[s][i])) == 1 else 'informative'


def likelihood(s, i):
    if not present[s][i]:
        return np.ones(3)
    low = min(pl[s][i])
    return np.array([f32(10.0 ** (-(x - low) / 10.0)) for x in pl[s][i]])


def shares(s, i):
    if not present[s][i]:
        return np.array([f32(0.5), f32(0.5)])
    g = gp[s][i]
    q = [g[0] + g[1] / 2, g[1] / 2 + g[2]]
    t = sum(q)
    return np.array([f32(x / t) for x in q]) if t > 0 else np.array([f32(0.5), f32(0.5)])


def background(i):
    p = round(afs[i], 6)
    q = [max(0.0, 1 - p), max(0.0, p)]
    t = sum(q)
    return np.array([f32(x / t) for x in q])


LIK = [[likelihood(s, i) for i in range(RECORDS)] for s in range(SAMPLES)]
SHARES = [[shares(s, i) for i in range(RECORDS)] for s in range(SAMPLES)]
BACK = [background(i) for i in range(RECORDS)]
ROW = [[row(s, i) for i in range(RECORDS)] for s in range(SAMPLES)]
SITES = [i for i in range(RECORDS) if sum(ROW[s][i] != 'absent' for s in range(SAMPLES)) >= 2]


def bound(d, n):
    mean = (d + 1) / (n + 2)
    return mean + 2 * math.sqrt(mean * (1 - mean) / (n + 3))


def compared(s):
    informative = [i for i in SITES if ROW[s][i] == 'informative']
    if len(informative) <= 500:
        return informative
    return [informative[j * len(informative) // 500] for j in range(500)]


def choose_samples(s):
    rows = compared(s)
    if not rows:
        return None
    scored = []
    for o in range(SAMPLES):
        if o == s:
            continue
        shared = [i for i in rows if ROW[o][i] == 'informative']
        d = sum(0.5 * float(np.abs(SHARES[s][i] - SHARES[o][i]).sum()) for i in shared)
        scored.append((bound(d, len(shared)), o))
    scored.sort()
    return [o for _, o in scored[:K]]


def choose_haplotypes(s, previous, store):
    rows = compared(s)
    offers = []
    for o in range(SAMPLES):
        if o == s:
            continue
        shared = [i for i in rows if ROW[o][i] == 'informative']
        totals = [0.0, 0.0]
        for i in shared:
            post = previous[s][i]
            lacking = [post[2], post[0]]   # P(genotype without allele 0), without allele 1
            for c in (0, 1):
                h = store[o][i][c]
                totals[c] += h[0] * lacking[0] + h[1] * lacking[1]
        scores = [bound(t, len(shared)) for t in totals]
        offers.append((scores[1], o, 1) if scores[1] < scores[0] else (scores[0], o, 0))
    offers.sort()
    return [(o, c) for _, o, c in offers[:K]]


def hidden_records(s):
    candidates = [i for i in SITES if ROW[s][i] == 'informative']
    u = Uniforms(1, WINDOW_KEY, s, 0)
    count = min(HIDDEN, len(candidates))
    for d in range(count):
        remaining = len(candidates) - d
        picked = d + min(int(u.next() * remaining), remaining - 1)
        candidates[d], candidates[picked] = candidates[picked], candidates[d]
    return sorted(candidates[:count])


def draw_path(liks, q_sites, uniforms):
    """The forward walk, then one path drawn backwards and each copy's haplotype along it."""
    n, T = len(SITES), q_sites[0].shape[0]
    fwd_d, fwd_p, Ms = [], [], []
    for j in range(n):
        q = q_sites[j]
        l0, l1, l2 = liks[j]
        M = np.array([[l0, l1], [l1, l2]])
        Ms.append(M)
        e_pair = q @ M @ q.T
        e_desc = q @ np.diag(M)
        if j == 0:
            d = np.full(T, SHARE / T)
            p = np.full((T, T), (1 - SHARE) / (T * T))
        else:
            gap = positions[SITES[j]] - positions[SITES[j - 1]]
            s_d, s_n = stays(gap, SHARE, RUN)
            keep = math.exp(-gap / SEGMENT)
            prev_d, prev_p = fwd_d[-1], fwd_p[-1]
            mass_d, mass_p = prev_d.sum(), prev_p.sum()
            d_in = keep * prev_d + (1 - keep) * mass_d / T
            r_, c_ = prev_p.sum(axis=1), prev_p.sum(axis=0)
            p_in = (keep * keep * prev_p + keep * (1 - keep) * (r_[:, None] + c_[None, :]) / T
                    + (1 - keep) ** 2 * mass_p / (T * T))
            d = s_d * d_in + (1 - s_n) * mass_p / T
            p = s_n * p_in + (1 - s_d) * mass_d / (T * T)
        d = d * e_desc
        p = p * e_pair
        total = d.sum() + p.sum()
        fwd_d.append(d / total); fwd_p.append(p / total)
    states = [None] * n
    for j in range(n - 1, -1, -1):
        if j == n - 1:
            w_d, w_p = fwd_d[j], fwd_p[j]
        else:
            gap = positions[SITES[j + 1]] - positions[SITES[j]]
            s_d, s_n = stays(gap, SHARE, RUN)
            keep = math.exp(-gap / SEGMENT)
            kind, k = states[j + 1]
            if kind == 'p':
                k1, k2 = k
                w_d = fwd_d[j] * (1 - s_d) / (T * T)
                w_p = np.empty((T, T))
                for a in range(T):
                    for b in range(T):
                        t = (1 - keep) * (1 - keep) / (T * T)
                        if a == k1:
                            t += keep * (1 - keep) / T
                        if b == k2:
                            t += keep * (1 - keep) / T
                        if a == k1 and b == k2:
                            t += keep * keep
                        w_p[a, b] = fwd_p[j][a, b] * s_n * t
            else:
                w_d = np.array([fwd_d[j][a] * s_d * ((1 - keep) / T + (keep if a == k else 0.0))
                                for a in range(T)])
                w_p = fwd_p[j] * (1 - s_n) / T
        index = pick(list(w_d) + list(w_p.ravel()), uniforms.next())
        states[j] = ('d', index) if index < T else ('p', divmod(index - T, T))
    haps = []
    for j, (kind, k) in enumerate(states):
        q, M = q_sites[j], Ms[j]
        joint = np.outer(q[k[0]], q[k[1]]) * M if kind == 'p' else np.diag(q[k] * np.diag(M))
        total = joint.sum()
        if total <= 0:
            haps.append((np.full(2, 0.5), np.full(2, 0.5)))
            continue
        joint = joint / total
        haps.append((joint.sum(axis=1), joint.sum(axis=0)))
    return haps


# ---- the rounds ----------------------------------------------------------------------------------
round_one = [choose_samples(s) for s in range(SAMPLES)]
hidden = [hidden_records(s) if round_one[s] is not None else [] for s in range(SAMPLES)]
store = [{i: (SHARES[s][i].astype(float), SHARES[s][i].astype(float)) for i in SITES}
         for s in range(SAMPLES)]
previous = [dict() for _ in range(SAMPLES)]
sums = [None] * SAMPLES
for round_ in range(1, ROUNDS + 1):
    new_store, new_previous, posts = [None] * SAMPLES, [dict() for _ in range(SAMPLES)], [None] * SAMPLES
    for s in range(SAMPLES):
        if round_one[s] is None:
            continue
        liks = [np.ones(3) if i in hidden[s] else LIK[s][i] for i in SITES]
        if round_ == 1:
            rows = [[SHARES[t][i] for t in round_one[s]] + [BACK[i]] for i in SITES]
        else:
            pairs = choose_haplotypes(s, previous, store)
            rows = [[store[o][i][c] for o, c in pairs] + [BACK[i]] for i in SITES]
        raw = [np.array(r, dtype=float) for r in rows]
        posteriors, _ = haplotypes([(l, None, None) for l in liks], [positions[i] for i in SITES],
                                   raw, SHARE, RUN, SEGMENT, MISMATCH)
        posts[s] = {i: post for i, post in zip(SITES, posteriors)}
        new_previous[s] = posts[s]
        blurred = [r * (1 - MISMATCH) + MISMATCH / 2 for r in raw]
        haps = draw_path(liks, blurred, Uniforms(2, WINDOW_KEY, s, round_))
        new_store[s] = {i: ((np.full(2, 0.5), np.full(2, 0.5)) if ROW[s][i] == 'absent' else h)
                        for i, h in zip(SITES, haps)}
    for s in range(SAMPLES):
        if round_one[s] is None:
            continue
        store[s] = new_store[s]
        previous[s] = new_previous[s]
        if round_ >= 2:
            if sums[s] is None:
                sums[s] = {i: posts[s][i].copy() for i in SITES}
            else:
                for i in SITES:
                    sums[s][i] += posts[s][i]

# ---- the fixture ---------------------------------------------------------------------------------
print(f'window {SAMPLES} {RECORDS} {SHARE} {CHOSEN} {ROUNDS} {WINDOW_KEY}')
for i in range(RECORDS):
    columns = []
    for s in range(SAMPLES):
        if not present[s][i]:
            columns.append('.')
        else:
            columns.append(','.join(map(str, pl[s][i])) + '|' + ','.join(f'{x:.3f}' for x in gp[s][i]))
    print('record', positions[i], f'{round(afs[i], 6)}', ' '.join(columns))
for s in range(SAMPLES):
    print('hidden', s, ' '.join(map(str, hidden[s])))
for s in range(SAMPLES):
    if sums[s] is None:
        continue
    for i in SITES:
        if ROW[s][i] == 'absent':
            continue
        mean = sums[s][i] / (ROUNDS - 1)
        if i in hidden[s]:
            mean = mean * LIK[s][i]
            mean = mean / mean.sum()
        print('posterior', i, s, ' '.join(f'{x:.17g}' for x in mean))
