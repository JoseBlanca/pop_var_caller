"""The stray share fitted per stratum and smoothed along reference repeat count, on HG002's real reads.

usage: uv run --with numpy --with scipy tmp/p3/smooth_fit.py COV PERIOD > tmp/p3/smooth_COV_pPERIOD.tsv

Per stratum (this period, one reference repeat count, >= 100 tracts): sim.py's read model, plus reads
calling sends to the outlier term — a non-whole number of repeats off, or a flank-deletion spelling
carrying under 20% of its tract's reads — each scored at stray/25 under every genotype.

The share's logit is penalised on its second difference along repeat count, strength kappa, fitted by
block coordinate descent (each stratum's own numbers given its neighbours' shares). kappa is chosen
by fitting on even tracts and scoring the odd ones; 'shared' is one share for the whole period.
"""
import sys
from collections import defaultdict
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, "/Users/jose/devel/pop_var_caller-str-group-b/tmp/p2")
import sim

W = "/Users/jose/devel/pop_var_caller-str-group-b"
cov, period = sys.argv[1], int(sys.argv[2])
MAX_TRACTS, MIN_TRACTS = 2000, 100
logit = lambda p: np.log(p / (1 - p))

# ---- evidence: per tract, whole-repeat bucket counts and the count of outlier-only reads
spell = defaultdict(lambda: defaultdict(int)); meta = {}
for line in open(f"{W}/tmp/p1/tracts_{cov}.tsv"):
    f = line.rstrip("\n").split("\t")
    if f[0] == "contig" or f[14] != "complete" or len(f[4]) != period: continue
    key = (f[0], f[1], f[2])
    meta[key] = len(f[7])
    spell[key][(len(f[13]) - len(f[7]), f[22])] += int(f[17])
strata = defaultdict(list)
for key, sp in spell.items():
    n = sum(sp.values()); counts = np.zeros(sim.NB); other = 0
    for (d, flank), k in sp.items():
        fd = sum(int(x) for x in flank.split(",")) if flank else 0
        if fd or d % period:
            if k / n < 0.2: other += k     # too rare to be an allele: calling scores it as an outlier
            continue                        # an allele spelled with a flank deletion: left out, as today
        counts[int(np.clip(d // period, -12, 12)) + 12] += k
    if counts.sum() > 0:
        strata[meta[key] // period].append((counts, other))
rng = np.random.default_rng(1)
reps = sorted(r for r, ts in strata.items() if len(ts) >= MIN_TRACTS)
data = {}
for r in reps:
    ts = strata[r]
    if len(ts) > MAX_TRACTS: ts = [ts[i] for i in rng.choice(len(ts), MAX_TRACTS, replace=False)]
    C = np.array([c for c, _ in ts]); O = np.array([o for _, o in ts], float)
    data[r] = dict(even=(C[0::2], O[0::2]), odd=(C[1::2], O[1::2]), all=(C, O))


def nll(x, C, O):
    """Stratum negative log-likelihood; x as sim.unpack."""
    stray = sim.unpack(x)[3]
    return sim.neg_ll(x, C) - O.sum() * np.log(stray / sim.NB)


def fit_period(which, kappa, shared=False):
    """Block coordinate descent; returns per-stratum parameter vectors."""
    m0 = np.full(len(sim.ALLELES), 0.005); m0[sim.SPAN] = 1 - 0.005 * (len(sim.ALLELES) - 1)
    x0 = np.concatenate([[logit(0.05), 0.0, logit(0.3), logit(0.003)], np.log(m0[:-1] / m0[-1])])
    X = {r: x0.copy() for r in reps}
    for sweep in range(6):
        moved = 0.0
        if shared:  # one share: optimise it given every stratum's other numbers, then each stratum's others
            f = lambda t: sum(nll(np.concatenate([X[r][:3], t, X[r][4:]]), *data[r][which]) for r in reps)
            t = minimize(f, [X[reps[0]][3]], method="L-BFGS-B").x
            for r in reps: X[r][3] = t[0]
        for i, r in enumerate(reps):
            th = lambda r2: X[r2][3]
            def obj(x):
                v = nll(x, *data[r][which])
                if kappa and not shared:
                    tt = {**{q: th(q) for q in reps}, r: x[3]}
                    for j in range(max(1, i - 1), min(len(reps) - 1, i + 2)):
                        v += kappa * (tt[reps[j - 1]] - 2 * tt[reps[j]] + tt[reps[j + 1]]) ** 2
                return v
            if shared:
                y = minimize(lambda z: obj(np.concatenate([z[:3], [X[r][3]], z[3:]])), np.delete(X[r], 3), method="L-BFGS-B", options=dict(maxiter=200)).x
                new = np.concatenate([y[:3], [X[r][3]], y[3:]])
            else:
                new = minimize(obj, X[r], method="L-BFGS-B", options=dict(maxiter=200)).x
            moved = max(moved, abs(new[3] - X[r][3])); X[r] = new
        if moved < 0.01: break
    return X


print(f"# {cov} period {period}: strata {reps}")
results = {}
for label, kappa, shared in (("independent", 0.0, False), ("kappa 1", 1.0, False), ("kappa 10", 10.0, False), ("kappa 100", 100.0, False), ("shared", 0.0, True)):
    X = fit_period("even", kappa, shared)
    held = sum(nll(X[r], *data[r]["odd"]) for r in reps)
    results[label] = (held, X)
    print(f"# {label:12s} held-out -lnL on odd tracts {held:12.1f}", flush=True)
best = min(results, key=lambda k: results[k][0])
print(f"# best by held-out: {best}")
kappa = {"independent": 0, "kappa 1": 1, "kappa 10": 10, "kappa 100": 100}.get(best, 0)
X = fit_period("all", float(kappa), best == "shared")
print("period\trepeats\ttracts\tstray\tslip\tshorter\tfalloff\toutlier_only_reads")
for r in reps:
    slip, shorter, falloff, stray, _ = sim.unpack(X[r])
    C, O = data[r]["all"]
    print(f"{period}\t{r}\t{len(C)}\t{stray:.5f}\t{slip:.4f}\t{shorter:.3f}\t{falloff:.3f}\t{int(O.sum())}", flush=True)
