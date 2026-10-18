"""Phase 1 of the outlier-share program: where reads land at HG002 truth-homozygous tracts, from the
psp (every tract, quiet ones included), and three descriptions of the per-tract counts.

usage: uv run --with numpy --with scipy tmp/p1/phase1.py COV [COV ...]

Per truth-homozygous tract, complete reads only, offsets from the true length:
  at       at the true length
  one      one repeat off (stutter's first step)
  two      exactly two repeats off
  far      three or more repeats off
  nonwhole a non-whole number of repeats off (period >= 2 only)
  flank    the read deletes flank bases beside the tract (counted here whatever its length)
'beyond' = two + far + nonwhole + flank: what a one-step stutter model cannot place cheaply.

For 'one' and for 'beyond', per period class, fits to the per-tract counts k of n:
  binomial      one rate
  beta-binomial a mean rate and a within-tract correlation rho
  two classes   a share w of noisy tracts with rate r_noisy, the rest at r_clean
and prints a tail table (tracts with k such reads, observed against each fit's expectation).
"""
import bisect, gzip, sys
from collections import defaultdict
import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln, betaln

W = "/Users/jose/devel/pop_var_caller-str-group-b"
TRUTH = "/Users/jose/devel/pop_var_caller/benchmarks/ssr_hg002/truth/HG002_GRCh38_TandemRepeats_v1.0.1_50000.vcf.gz"

truth = defaultdict(list)
with gzip.open(TRUTH, "rt") as fh:
    for l in fh:
        if l[0] == "#": continue
        f = l.split("\t"); gt = f[9].split(":")[0].replace("/", "|").split("|"); alts = f[4].split(",")
        truth[f[0]].append((int(f[1]), [0 if a in "0." else len(alts[int(a) - 1]) - len(f[3]) for a in gt[:2]]))
for c in truth: truth[c].sort()


def load(cov):
    tracts = {}
    for line in open(f"{W}/tmp/p1/tracts_{cov}.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] == "contig" or f[14] != "complete": continue
        key = (f[0], int(f[1]), int(f[2]))
        t = tracts.setdefault(key, dict(motif=f[4], ref=len(f[7]), reads=[]))
        fd = sum(int(x) for x in f[22].split(",")) if f[22] else 0
        t["reads"].append((len(f[13]) - len(f[7]) - fd, fd > 0, int(f[17])))
    return tracts


def classify(cov):
    rows = defaultdict(list)  # period class -> (n, one, beyond, cats dict)
    for (c, s, e), t in load(cov).items():
        recs = truth.get(c, []); i = bisect.bisect_left(recs, (s - 1, [])); hap = [0, 0]
        while i < len(recs) and recs[i][0] <= e:
            hap[0] += recs[i][1][0]; hap[1] += recs[i][1][1]; i += 1
        if hap[0] != hap[1]: continue
        p = len(t["motif"]); tv = hap[0]
        cats = defaultdict(int)
        for d, flank, k in t["reads"]:
            off = d - tv
            if flank: cats["flank"] += k
            elif off == 0: cats["at"] += k
            elif off % p: cats["nonwhole"] += k
            elif abs(off) == p: cats["one"] += k
            elif abs(off) == 2 * p: cats["two"] += k
            else: cats["far"] += k
        n = sum(cats.values())
        if n == 0: continue
        cls = "p1" if p == 1 else ("p2" if p == 2 else "p3+")
        kind = "ref" if tv == 0 else "alt"
        rows[cls].append((n, cats["one"], cats["two"] + cats["far"] + cats["nonwhole"] + cats["flank"], dict(cats), kind))
    return rows


def ll_binom(k, n, r): return np.sum(k * np.log(r) + (n - k) * np.log1p(-r))
def ll_bb(k, n, m, rho):
    a = m * (1 - rho) / rho; b = (1 - m) * (1 - rho) / rho
    return np.sum(betaln(k + a, n - k + b) - betaln(a, b))
def ll_two(k, n, w, r0, r1):
    l0 = np.log1p(-w) + k * np.log(r0) + (n - k) * np.log1p(-r0)
    l1 = np.log(w) + k * np.log(r1) + (n - k) * np.log1p(-r1)
    return np.sum(np.logaddexp(l0, l1)), l1 - np.logaddexp(l0, l1)


def sig(x): return 1 / (1 + np.exp(-x))


def fit(k, n):
    r = k.sum() / n.sum()
    out = {"binomial": (ll_binom(k, n, r), dict(rate=r))}
    best = minimize(lambda x: -ll_bb(k, n, sig(x[0]), sig(x[1])), [np.log(r / (1 - r)), -3], method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-6, maxiter=4000))
    out["beta-binomial"] = (-best.fun, dict(rate=sig(best.x[0]), rho=sig(best.x[1])))
    bestt = None
    for w0 in (0.01, 0.05, 0.2):
        for r1 in (0.05, 0.2, 0.4):
            x0 = [np.log(w0 / (1 - w0)), np.log(r / 2 / (1 - r / 2)), np.log(r1 / (1 - r1))]
            b = minimize(lambda x: -ll_two(k, n, sig(x[0]), sig(x[1]), sig(x[2]))[0], x0, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-7, maxiter=8000))
            if bestt is None or b.fun < bestt.fun: bestt = b
    w, r0, r1 = (sig(v) for v in bestt.x)
    if r0 > r1: w, r0, r1 = 1 - w, r1, r0
    out["two classes"] = (-bestt.fun, dict(share_noisy=w, rate_clean=r0, rate_noisy=r1))
    return out


def tail(k, n, fits, maxk=6):
    from scipy.stats import binom, betabinom
    obs = [np.sum(k == j) for j in range(maxk)] + [np.sum(k >= maxk)]
    exp = {}
    _, p = fits["binomial"]; e = [binom.pmf(j, n, p["rate"]).sum() for j in range(maxk)]; exp["binomial"] = e + [len(n) - sum(e)]
    _, p = fits["beta-binomial"]; a = p["rate"] * (1 - p["rho"]) / p["rho"]; b = (1 - p["rate"]) * (1 - p["rho"]) / p["rho"]
    e = [betabinom.pmf(j, n, a, b).sum() for j in range(maxk)]; exp["beta-binomial"] = e + [len(n) - sum(e)]
    _, p = fits["two classes"]
    e = [((1 - p["share_noisy"]) * binom.pmf(j, n, p["rate_clean"]) + p["share_noisy"] * binom.pmf(j, n, p["rate_noisy"])).sum() for j in range(maxk)]
    exp["two classes"] = e + [len(n) - sum(e)]
    print("        reads k   observed   binomial   beta-binomial   two classes")
    for j in range(maxk + 1):
        lab = f"{j}" if j < maxk else f">={maxk}"
        print(f"        {lab:>7s}  {obs[j]:9d}  {exp['binomial'][j]:9.1f}  {exp['beta-binomial'][j]:14.1f}  {exp['two classes'][j]:12.1f}")


for cov in sys.argv[1:]:
    rows = classify(cov)
    for cls in ("p1", "p2", "p3+"):
        R = rows[cls]
        n = np.array([r[0] for r in R], float)
        tot = defaultdict(int)
        for r in R:
            for kk, v in r[3].items(): tot[kk] += v
        N = n.sum()
        print(f"\n=== {cov} {cls}: {len(R)} truth-homozygous tracts ({sum(r[4]=='ref' for r in R)} reference, {sum(r[4]=='alt' for r in R)} non-reference), {int(N)} complete reads, mean {N/len(R):.1f} a tract")
        print("    reads: " + "  ".join(f"{c} {tot[c]/N:.4f}" for c in ("at", "one", "two", "far", "nonwhole", "flank")))
        for what, idx in (("one repeat off", 1), ("beyond one repeat", 2)):
            k = np.array([r[idx] for r in R], float)
            fits = fit(k, n)
            print(f"  {what}:")
            for name, (ll, par) in fits.items():
                print(f"    {name:14s} ln L {ll:11.1f}   " + "  ".join(f"{a}={b:.4g}" for a, b in par.items()))
            tail(k, n, fits)
