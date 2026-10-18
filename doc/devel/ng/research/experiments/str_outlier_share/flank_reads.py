"""Flank-deletion reads at homopolymer tracts, against truth (psp dump, complete reads, edge fix on).
Per read with a flank deletion: total length change = tract change - flank bases deleted. Does it
equal one of the tract's true allele lengths? Split by whether the spelling is a minority (<20% of the
tract's reads, what calling sends to the outlier term) and by the truth genotype's kind."""
import sys, bisect
sys.argv, covs = sys.argv[:1], sys.argv[1:]
sys.path.insert(0, "/Users/jose/devel/pop_var_caller-str-group-b/tmp/p1")
from collections import Counter, defaultdict
import phase1 as P
W = "/Users/jose/devel/pop_var_caller-str-group-b"
for cov in covs:
    tr = defaultdict(lambda: defaultdict(int)); meta = {}
    for line in open(f"{W}/tmp/p1/tracts_{cov}.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] == "contig" or f[14] != "complete" or len(f[4]) != 1: continue
        key = (f[0], int(f[1]), int(f[2])); meta[key] = (len(f[7]), f[5][-1:], f[6][:1], f[4])
        tr[key][(len(f[13]) - len(f[7]), f[22])] += int(f[17])
    tally = Counter(); reads_all = 0; kinds = Counter(); side = Counter(); flankbase = Counter()
    for key, sp in tr.items():
        c, s, e = key
        recs = P.truth.get(c, []); i = bisect.bisect_left(recs, (s - 1, [])); hap = [0, 0]
        while i < len(recs) and recs[i][0] <= e:
            hap[0] += recs[i][1][0]; hap[1] += recs[i][1][1]; i += 1
        kind = "ref" if hap == [0, 0] else ("het" if hap[0] != hap[1] else ("hom_del" if hap[0] < 0 else "hom_ins"))
        n = sum(sp.values()); reads_all += n; kinds[kind] += n
        for (d, flank), k in sp.items():
            if not flank: continue
            l, r = (int(x) for x in flank.split(","))
            total = d - l - r
            minority = k / n < 0.2
            match = total in hap
            tally[(kind, "minority" if minority else "majority", "matches a true allele" if match else "matches none")] += k
            if minority:
                side["left" if l else "right"] += k
                _, lb, rb, motif = meta[key]
                flankbase[("left" if l else "right", "flank base = repeat base" if (lb if l else rb) == motif else "flank base differs")] += k
    print(f"== {cov} homopolymers: {reads_all} complete reads")
    for kind in ("ref", "het", "hom_del", "hom_ins"):
        tot = kinds[kind]
        rows = {k[1:]: v for k, v in tally.items() if k[0] == kind}
        print(f"  truth {kind:7s} ({tot:7d} reads): " + "; ".join(f"{a} {b}: {v} ({1e4*v/tot:.1f}/10k)" for (a, b), v in sorted(rows.items())))
    print("  minority flank reads by side:", dict(side), " flank base:", dict(flankbase))
