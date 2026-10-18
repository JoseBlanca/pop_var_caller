"""Truth-free fit of the stray-read share on HG002's real reads, per stratum (period, reference repeats).

usage: uv run --with numpy --with scipy tmp/p3/real_fit.py COV > tmp/p3/real_fit_COV.tsv

Reads: complete reads from the psp dump (tmp/p1/tracts_COV.tsv), offsets from the reference length in
whole repeats, clipped to -12..+12. Reads a non-whole number of repeats off, and reads with a flank
deletion, are not fitted (today's stutter fit does not see them either) and are counted per stratum.
The fit is sim.py's: slip, direction, fall-off, stray share and 19 allele shares, summed over genotypes.
"""
import sys
from collections import defaultdict
import numpy as np
sys.path.insert(0, "/Users/jose/devel/pop_var_caller-str-group-b/tmp/p2")
import sim

W = "/Users/jose/devel/pop_var_caller-str-group-b"
cov = sys.argv[1]
MAX_TRACTS = 3000
MIN_TRACTS = 100

tracts = {}
for line in open(f"{W}/tmp/p1/tracts_{cov}.tsv"):
    f = line.rstrip("\n").split("\t")
    if f[0] == "contig" or f[14] != "complete":
        continue
    p = len(f[4])
    if not 1 <= p <= 6:
        continue
    key = (f[0], f[1], f[2])
    t = tracts.setdefault(key, dict(p=p, reps=len(f[7]) // p, counts=np.zeros(sim.NB), nonwhole=0, flank=0))
    k = int(f[17])
    fd = sum(int(x) for x in f[22].split(",")) if f[22] else 0
    d = len(f[13]) - len(f[7])
    if fd:
        t["flank"] += k
    elif d % p:
        t["nonwhole"] += k
    else:
        t["counts"][int(np.clip(d // p, -12, 12)) + 12] += k

strata = defaultdict(list)
for t in tracts.values():
    if t["counts"].sum() > 0:
        strata[(t["p"], t["reps"])].append(t)

print("period\trepeats\ttracts\treads\tstray\tslip\tshorter\tfalloff\tnonwhole_share\tflank_share\tfar_share")
rng = np.random.default_rng(1)
for (p, reps), ts in sorted(strata.items()):
    if len(ts) < MIN_TRACTS:
        continue
    if len(ts) > MAX_TRACTS:
        ts = [ts[i] for i in rng.choice(len(ts), MAX_TRACTS, replace=False)]
    C = np.array([t["counts"] for t in ts])
    n = C.sum()
    other_nonwhole = sum(t["nonwhole"] for t in ts)
    other_flank = sum(t["flank"] for t in ts)
    # reads 3+ repeats from the tract's own most-read length: a model-free stand-in for strays
    far = sum(t["counts"][np.abs(np.arange(sim.NB) - np.argmax(t["counts"])) >= 3].sum() for t in ts)
    (slip, shorter, falloff, stray, m), ll = sim.fit(C)
    tot = n + other_nonwhole + other_flank
    print(f"{p}\t{reps}\t{len(ts)}\t{int(n)}\t{stray:.5f}\t{slip:.4f}\t{shorter:.3f}\t{falloff:.3f}\t{other_nonwhole/tot:.5f}\t{other_flank/tot:.5f}\t{far/n:.5f}", flush=True)
