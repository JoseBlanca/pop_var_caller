"""Where far reads land at truth-homozygous homopolymer tracts (clean ones: <= 35% of reads off the
allele), by reference length. Offsets in repeats from the true allele; reads >= 3 off."""
import sys
sys.argv, covs = sys.argv[:1], sys.argv[1:]
sys.path.insert(0, "/Users/jose/devel/pop_var_caller-str-group-b/tmp/p1")
from collections import Counter, defaultdict
import bisect
import phase1 as P
for cov in covs:
    by = defaultdict(Counter); n_reads = Counter()
    for (c, s, e), t in P.load(cov).items():
        if len(t["motif"]) != 1: continue
        recs = P.truth.get(c, []); i = bisect.bisect_left(recs, (s - 1, [])); hap = [0, 0]
        while i < len(recs) and recs[i][0] <= e:
            hap[0] += recs[i][1][0]; hap[1] += recs[i][1][1]; i += 1
        if hap[0] != hap[1]: continue
        reads = [(d - hap[0], k) for d, fl, k in t["reads"] if not fl]
        n = sum(k for _, k in reads)
        if n == 0 or sum(k for o, k in reads if o != 0) / n > 0.35: continue
        L = t["ref"]; b = "8-12" if L <= 12 else ("13-17" if L <= 17 else ("18-22" if L <= 22 else "23+"))
        n_reads[b] += n
        for o, k in reads:
            if abs(o) >= 2: by[b][o] += k
    print(f"== {cov}: reads per 10,000 at each offset from the true allele, clean truth-homozygous homopolymers")
    offs = [-10, -8, -6, -5, -4, -3, -2, 2, 3, 4, 5, 6, 8, 10]
    print("length  reads    " + " ".join(f"{o:>6d}" for o in offs) + "   |>=6| total")
    for b in ("8-12", "13-17", "18-22", "23+"):
        n = n_reads[b]
        far6 = sum(k for o, k in by[b].items() if abs(o) >= 6)
        print(f"{b:6s} {n:7d}  " + " ".join(f"{1e4*by[b][o]/n:6.1f}" for o in offs) + f"   {1e4*far6/n:6.1f}")
