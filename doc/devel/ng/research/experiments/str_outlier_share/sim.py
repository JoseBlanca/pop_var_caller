"""Phase 2/3 prototype: simulate repeat-tract reads and fit the stray-read share truth-free.

The generator, per stratum (one period, one reference repeat count):
  - allele-length shares m over offsets -9..+9 repeats (default: a homopolymer stratum's fitted shares);
  - each sample's genotype drawn from a Dirichlet-multinomial with concentration 1 over m;
  - reads per (tract, sample) Poisson(depth); each read from one of the two copies at random;
  - a read is a stray with probability `stray` (uniform over offsets -12..+12), otherwise it slips with
    probability `slip`, shorter with probability `shorter`, by k repeats with chance (1-f) f^(k-1);
  - a share `noisy` of tracts is noisy in every sample: a fraction `noisy_frac` of their reads comes from
    one other offset, fixed per tract, 1-4 repeats from the reference.

The fit (what estimate-parameters would do, never shown the truth): maximise, over slip, shorter,
fall-off, stray and the 19 shares, the sum over (tract, sample) of
  ln sum_g prior(g | m) prod_reads [ (1 - stray) S(read - allele mixture) + stray / 25 ].

usage: uv run --with numpy --with scipy tmp/p2/sim.py
"""
import itertools, sys, time
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

SPAN = 9                      # allele offsets -9..9
OFFS = np.arange(-12, 13)     # read offsets recorded, 25 buckets
NB = len(OFFS)
ALLELES = np.arange(-SPAN, SPAN + 1)
GENOS = [(a, b) for a in range(len(ALLELES)) for b in range(a, len(ALLELES))]
# a homopolymer stratum's fitted shares (period 1, ~15 repeats), rounded; centre = reference
DEFAULT_SHARES = np.array([1e-4, 1e-4, 1e-4, 2e-4, 2e-4, 5e-4, 1e-3, 3e-3, 2.5e-2, 0.90, 4e-2, 1.5e-2, 6e-3, 3e-3, 2e-3, 1e-3, 5e-4, 3e-4, 2e-4])
DEFAULT_SHARES /= DEFAULT_SHARES.sum()


def stutter_matrix(slip, shorter, falloff):
    """P(read offset | allele offset), without strays: rows alleles (19), columns read buckets (25)."""
    M = np.zeros((len(ALLELES), NB))
    for i, a in enumerate(ALLELES):
        for j, x in enumerate(OFFS):
            d = x - a
            if d == 0:
                M[i, j] = 1 - slip
            elif 1 <= abs(d) <= 10:
                side = shorter if d < 0 else 1 - shorter
                M[i, j] = slip * side * (1 - falloff) * falloff ** (abs(d) - 1)
    return M


def emission(slip, shorter, falloff, stray):
    return (1 - stray) * stutter_matrix(slip, shorter, falloff) + stray / NB


def genotype_log_prior(m, conc=1.0):
    out = np.empty(len(GENOS))
    for g, (a, b) in enumerate(GENOS):
        out[g] = np.log(m[a] * (conc * m[a] + 1) / (conc + 1)) if a == b else np.log(2 * conc * m[a] * m[b] / (conc + 1))
    return out


def simulate(rng, tracts, samples, depth, slip=0.04, shorter=0.66, falloff=0.23, stray=0.003, noisy=0.0, noisy_frac=0.3, shares=DEFAULT_SHARES):
    """Counts per (tract, sample) over the 25 read buckets, shape (tracts*samples, 25)."""
    E = emission(slip, shorter, falloff, stray)
    counts = np.zeros((tracts, samples, NB))
    noisy_tract = rng.random(tracts) < noisy
    other = rng.choice([-4, -3, -2, -1, 1, 2, 3, 4], size=tracts) + SPAN + 3  # bucket index of the extra length
    for t in range(tracts):
        for s in range(samples):
            # Dirichlet-multinomial with concentration 1: draw a sample-specific frequency, then two copies
            p = rng.dirichlet(shares * 1.0 + 1e-12)
            a, b = rng.choice(len(ALLELES), size=2, p=p)
            n = rng.poisson(depth)
            if n == 0:
                continue
            k = rng.binomial(n, noisy_frac) if noisy_tract[t] else 0
            copies = rng.integers(0, 2, size=n - k)
            na = int((copies == 0).sum())
            counts[t, s] += rng.multinomial(na, E[a]) + rng.multinomial(n - k - na, E[b])
            counts[t, s, other[t]] += k
    return counts.reshape(-1, NB)


def unpack(x):
    sig = lambda v: 1 / (1 + np.exp(-v))
    slip, shorter, falloff, stray = sig(x[0]), sig(x[1]), sig(x[2]), sig(x[3])
    m = np.exp(np.append(x[4:], 0.0)); m /= m.sum()
    return slip, shorter, falloff, stray, m


def neg_ll(x, C):
    slip, shorter, falloff, stray, m = unpack(x)
    logE = np.log(emission(slip, shorter, falloff, stray))
    # per genotype, the read log-probability under the half-half mixture of its two alleles
    G = np.array([np.logaddexp(logE[a], logE[b]) - np.log(2) for a, b in GENOS])  # (190, 25)
    L = C @ G.T + genotype_log_prior(m)                                            # (rows, 190)
    return -logsumexp(L, axis=1).sum()


def fit(C, fixed_stray=None):
    logit = lambda p: np.log(p / (1 - p))
    m0 = np.full(len(ALLELES), 0.005); m0[SPAN] = 1 - 0.005 * (len(ALLELES) - 1)
    x0 = np.concatenate([[logit(0.05), 0.0, logit(0.3), logit(0.01)], np.log(m0[:-1] / m0[-1])])
    if fixed_stray is not None:
        f = lambda y: neg_ll(np.concatenate([y[:3], [logit(fixed_stray)], y[3:]]), C)
        y0 = np.delete(x0, 3)
        o = minimize(f, y0, method="L-BFGS-B", options=dict(maxiter=500))
        x = np.concatenate([o.x[:3], [logit(fixed_stray)], o.x[3:]])
        return unpack(x), -o.fun
    o = minimize(neg_ll, x0, args=(C,), method="L-BFGS-B", options=dict(maxiter=500))
    return unpack(o.x), -o.fun


if __name__ == "__main__":
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    reads_budget = 600_000  # keep each fit's read count comparable to a stratum's census cap
    print("samples depth noisy | true stray -> fitted stray | slip (true 0.04) shorter (0.66) falloff (0.23) | seconds")
    for noisy in (0.0, 0.03):
        for samples, depth in itertools.product((1, 10, 100), (3, 10, 30, 100, 300)):
            tracts = int(min(2000, max(50, reads_budget / (samples * depth))))
            if samples * tracts * depth > 3_000_000:
                continue
            for stray in (0.001, 0.003, 0.01):
                C = simulate(rng, tracts, samples, depth, stray=stray, noisy=noisy)
                t0 = time.time()
                (slip, shorter, falloff, est, m), ll = fit(C)
                print(f"{samples:5d} {depth:5d} {noisy:5.2f} | {stray:.3f} -> {est:.4f} | {slip:.4f} {shorter:.3f} {falloff:.3f} | {tracts} tracts {time.time()-t0:5.1f}s", flush=True)
