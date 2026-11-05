"""Natural cohort: biallelic SNPs of the truth, by depth group (10 samples each: 30x, 20x, 10x, 5x)
and true genotype: right / wrong / missing / no record. usage: score_nat.py calls.vcf.gz truth.vcf.gz"""
import gzip, sys
from collections import Counter
def load(path):
    out, samples = {}, None
    for line in gzip.open(path, "rt"):
        if line.startswith("##"): continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#"): samples = f[9:]; continue
        out[(f[0], int(f[1]))] = ([f[3]] + f[4].split(","), f[8].split(":"), dict(zip(samples, f[9:])))
    return out, samples
T, samples = load(sys.argv[2]); C, _ = load(sys.argv[1])
groups = (sys.argv[3] if len(sys.argv) > 3 else "30x,20x,10x,5x").split(",")
n = Counter()
for key, (al, _, cols) in T.items():
    if len(al) != 2 or any(len(a) != 1 for a in al): continue
    for i, s in enumerate(samples):
        g = groups[i // 10]
        t = sorted(al[int(x)] for x in cols[s].split("|"))
        kind = "het" if t[0] != t[1] else ("ref" if t[0] == al[0] else "alt")
        if key not in C or C[key][0][0] != al[0]:
            n[(g, kind, "no record")] += 1; continue
        cal, keys, ccols = C[key]
        gt = dict(zip(keys, ccols[s].split(":")))["GT"]
        if "." in gt: out = "missing"
        else:
            c = sorted(cal[int(x)] for x in gt.replace("|", "/").split("/"))
            out = "right" if c == t else "wrong"
        n[(g, kind, out)] += 1
print("group kind     right   wrong missing  norec")
for g in groups:
    for kind in ("het", "alt", "ref"):
        tot = sum(n[(g, kind, o)] for o in ("right", "wrong", "missing", "no record"))
        print(f"{g:>4} {kind}  " + " ".join(f"{100*n[(g, kind, o)]/max(1,tot):6.1f}%" for o in ("right", "wrong", "missing", "no record")) + f"  of {tot}")
