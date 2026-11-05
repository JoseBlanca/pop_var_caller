"""F2 offspring at SNPs where the parents are different homozygotes: right / wrong / missing by
true genotype, plus mean GP of the true genotype. usage: score.py calls.vcf.gz truth.vcf.gz"""
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
n = Counter(); gp = Counter()
for key, (al, _, cols) in T.items():
    if any(len(a) != 1 for a in al): continue
    g0, g1 = cols["sim000"].replace("|", "/"), cols["sim001"].replace("|", "/")
    if g0 not in ("0/0", "1/1") or g1 not in ("0/0", "1/1") or g0 == g1: continue
    for s in samples[2:]:
        t = sorted(al[int(x)] for x in cols[s].split("|"))
        kind = "het" if t[0] != t[1] else "hom"
        if key not in C or C[key][0][:1] != al[:1]:
            n[(kind, "no record")] += 1; continue
        cal, keys, ccols = C[key]
        v = dict(zip(keys, ccols[s].split(":")))
        if "." in v["GT"]:
            out = "missing"
        else:
            c = sorted(cal[int(x)] for x in v["GT"].replace("|", "/").split("/"))
            out = "right" if c == t else "wrong"
        n[(kind, out)] += 1
        if v.get("GP", ".") != "." and len(cal) == 2:
            p = [float(x) for x in v["GP"].split(",")]
            idx = 1 if kind == "het" else (0 if t[0] == cal[0] else 2)
            gp[kind] += p[idx]; gp[(kind, "n")] += 1
for kind in ("het", "hom"):
    tot = sum(v for (k, _), v in n.items() if k == kind)
    print(f"{kind}: " + "  ".join(f"{o} {n[(kind, o)]} ({100*n[(kind, o)]/tot:.1f}%)" for o in ("right", "wrong", "missing", "no record"))
          + f"  mean GP(true) {gp[kind]/max(1,gp[(kind,'n')]):.3f}")
