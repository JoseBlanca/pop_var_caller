"""A heavier stutter tail: one geometric against a mixture of two, with the smoothed per-stratum stray
share, on HG002's real reads.

usage: uv run --with numpy --with scipy tmp/p3/tail_fit.py COV PERIOD KAPPA > tmp/p3/tail_COV_pPERIOD.tsv

Slip sizes, given a slip in one direction:
  geometric  P(k) = (1 - f) f^(k-1)
  mixture    P(k) = (1 - w) (1 - f) f^(k-1) + w (1 - g) g^(k-1),  g > f  (w and g fitted per stratum)
Evidence and the stray share's smoothing (second-difference penalty KAPPA on its logit) as in
smooth_fit.py. Both shapes are fitted on even tracts and scored on odd ones, then the mixture on all.
"""
import os, sys
from collections import defaultdict
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

sys.path.insert(0, "/Users/jose/devel/pop_var_caller-str-group-b/tmp/p2")
import sim

W = "/Users/jose/devel/pop_var_caller-str-group-b"
cov, period, KAPPA = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
MAX_TRACTS, MIN_TRACTS = 2000, 100
sig = lambda v: 1 / (1 + np.exp(-v))
logit = lambda p: np.log(p / (1 - p))
D = sim.OFFS[None, :] - sim.ALLELES[:, None]            # read offset minus allele offset, (19, 25)
K = np.abs(D)


def stutter(slip, shorter, f, w, g):
    side = np.where(D < 0, shorter, 1 - shorter)
    size = (1 - w) * (1 - f) * f ** np.maximum(K - 1, 0) + w * (1 - g) * g ** np.maximum(K - 1, 0)
    M = np.where(K == 0, 1 - slip, np.where(K <= 10, slip * side * size, 0.0))
    return M


GENO_A = np.array([a for a, b in sim.GENOS]); GENO_B = np.array([b for a, b in sim.GENOS])


def nll(x, C, O, mixture):
    slip, shorter, f, stray = sig(x[0]), sig(x[1]), sig(x[2]), sig(x[3])
    if mixture:
        w, g = sig(x[4]), f + (1 - f) * sig(x[5])   # g between f and 1
        m_x = x[6:]
    else:
        w, g = 0.0, 0.5
        m_x = x[4:]
    m = np.exp(np.append(m_x, 0.0)); m /= m.sum()
    logE = np.log((1 - stray) * stutter(slip, shorter, f, w, g) + stray / sim.NB)
    G = np.logaddexp(logE[GENO_A], logE[GENO_B]) - np.log(2)
    L = C @ G.T + sim.genotype_log_prior(m)
    return -logsumexp(L, axis=1).sum() - O.sum() * np.log(stray / sim.NB)


# ---- evidence (as smooth_fit.py)
spell = defaultdict(lambda: defaultdict(int)); meta = {}
for line in open(f"{W}/tmp/p1/tracts_{cov}.tsv"):
    f = line.rstrip("\n").split("\t")
    if f[0] == "contig" or f[14] != "complete" or len(f[4]) != period: continue
    key = (f[0], f[1], f[2]); meta[key] = len(f[7]); spell[key][(len(f[13]) - len(f[7]), f[22])] += int(f[17])
strata = defaultdict(list)
for key, sp in spell.items():
    n = sum(sp.values()); counts = np.zeros(sim.NB); other = 0
    for (d, flank), k in sp.items():
        fd = sum(int(x) for x in flank.split(",")) if flank else 0
        if fd or d % period:
            if k / n < 0.2 and not os.environ.get('NO_OTHER'): other += k
            continue
        counts[int(np.clip(d // period, -12, 12)) + 12] += k
    if counts.sum() > 0: strata[meta[key] // period].append((counts, other))
rng = np.random.default_rng(1)
reps = sorted(r for r, ts in strata.items() if len(ts) >= MIN_TRACTS)
data = {}
for r in reps:
    ts = strata[r]
    if len(ts) > MAX_TRACTS: ts = [ts[i] for i in rng.choice(len(ts), MAX_TRACTS, replace=False)]
    C = np.array([c for c, _ in ts]); O = np.array([o for _, o in ts], float)
    data[r] = dict(even=(C[0::2], O[0::2]), odd=(C[1::2], O[1::2]), all=(C, O))


def fit_period(which, mixture):
    m0 = np.full(len(sim.ALLELES), 0.005); m0[sim.SPAN] = 1 - 0.005 * (len(sim.ALLELES) - 1)
    head = [logit(0.05), 0.0, logit(0.1), logit(0.003)] + ([logit(0.05), 0.0] if mixture else [])
    x0 = np.concatenate([head, np.log(m0[:-1] / m0[-1])])
    X = {r: x0.copy() for r in reps}
    for sweep in range(6):
        moved = 0.0
        for i, r in enumerate(reps):
            def obj(x):
                v = nll(x, *data[r][which], mixture)
                if KAPPA:
                    tt = {**{q: X[q][3] for q in reps}, r: x[3]}
                    for j in range(max(1, i - 1), min(len(reps) - 1, i + 2)):
                        v += KAPPA * (tt[reps[j - 1]] - 2 * tt[reps[j]] + tt[reps[j + 1]]) ** 2
                return v
            new = minimize(obj, X[r], method="L-BFGS-B", options=dict(maxiter=300)).x
            moved = max(moved, abs(new[3] - X[r][3])); X[r] = new
        if moved < 0.01: break
    return X


print(f"# {cov} period {period} kappa {KAPPA}: strata {reps}")
held = {}
for name, mixture in (("geometric", False), ("mixture", True)):
    X = fit_period("even", mixture)
    held[name] = sum(nll(X[r], *data[r]["odd"], mixture) for r in reps)
    print(f"# {name:10s} held-out -lnL on odd tracts {held[name]:12.1f}", flush=True)
print(f"# mixture minus geometric, held-out: {held['mixture'] - held['geometric']:.1f} (negative = mixture better)")
for name, mixture in (("geometric", False), ("mixture", True)):
    X = fit_period("all", mixture)
    print(f"## {name}")
    print("period\trepeats\ttracts\tstray\tslip\tshorter\tfalloff\ttail_share\ttail_falloff")
    for r in reps:
        x = X[r]; f = sig(x[2])
        w, g = (sig(x[4]), f + (1 - f) * sig(x[5])) if mixture else (0.0, float("nan"))
        print(f"{period}\t{r}\t{len(data[r]['all'][0])}\t{sig(x[3]):.5f}\t{sig(x[0]):.4f}\t{sig(x[1]):.3f}\t{f:.3f}\t{w:.4f}\t{g:.3f}", flush=True)
