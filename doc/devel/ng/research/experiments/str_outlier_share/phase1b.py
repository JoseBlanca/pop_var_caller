"""Three tract classes for reads beyond one repeat, and what the clean class's reads are."""
import sys
sys.argv, covs = sys.argv[:1], sys.argv[1:]
sys.path.insert(0, "/Users/jose/devel/pop_var_caller-str-group-b/tmp/p1")
import numpy as np
from scipy.optimize import minimize
from scipy.stats import binom
import phase1 as P
def sig(x): return 1 / (1 + np.exp(-x))
def ll3(k, n, x):
    w = np.exp(np.append(x[:2], 0)); w /= w.sum(); r = sig(x[2:])
    comps = [np.log(w[i]) + k * np.log(r[i]) + (n - k) * np.log1p(-r[i]) for i in range(3)]
    L = np.logaddexp.reduce(comps, axis=0)
    return L.sum(), [np.exp(c - L) for c in comps], w, r
for cov in covs:
    rows = P.classify(cov)
    for cls in ("p1", "p2", "p3+"):
        R = rows[cls]; n = np.array([r[0] for r in R], float); k = np.array([r[2] for r in R], float)
        best = None
        for a in (-1, 1):
            for b in (-3, -1):
                x0 = [a, b, np.log(0.003/0.997), np.log(0.1/0.9), np.log(0.6/0.4)]
                o = minimize(lambda x: -ll3(k, n, x)[0], x0, method="Nelder-Mead", options=dict(maxiter=20000, xatol=1e-7, fatol=1e-7))
                if best is None or o.fun < best.fun: best = o
        L, post, w, r = ll3(k, n, best.x)
        order = np.argsort(r); w, r = w[order], r[order]; post = [post[i] for i in order]
        print(f"{cov} {cls}: three classes ln L {L:.1f}: " + "  ".join(f"share {w[i]:.4f} at rate {r[i]:.4f}" for i in range(3)))
        obs = [np.sum(k == j) for j in range(6)] + [np.sum(k >= 6)]
        e = [sum(w[i] * binom.pmf(j, n, r[i]).sum() for i in range(3)) for j in range(6)]; e.append(len(n) - sum(e))
        print("    tail k=0..5,>=6  observed " + " ".join(f"{o:6d}" for o in obs) + "\n                     three classes " + " ".join(f"{x:6.1f}" for x in e))
        clean = post[0] > 0.9
        tot = {}
        for r_, c in zip(R, clean):
            if c:
                for kk, v in r_[3].items(): tot[kk] = tot.get(kk, 0) + v
        N = sum(tot.values())
        print("    clean tracts' reads: " + "  ".join(f"{c} {tot.get(c,0)/N:.5f}" for c in ("one", "two", "far", "nonwhole", "flank")))
        for i, lab in ((1, "middle"), (2, "het-like")):
            sel = post[i] > 0.5
            print(f"    {lab} tracts: {int(sel.sum())}, of them non-reference truth {sum(1 for r_, s in zip(R, sel) if s and r_[4]=='alt')}")
