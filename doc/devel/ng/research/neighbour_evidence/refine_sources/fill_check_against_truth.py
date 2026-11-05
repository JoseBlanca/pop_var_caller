"""The pass's fill check against the truth. usage: fillcheck.py on.vcf.gz truth.vcf.gz groups
groups: comma list, 10 samples a group in sample order (or 'f2' for 2 parents + 20 offspring).
True fill error: biallelic SNPs where the sample's DP is 0 and the pass wrote a genotype, against the
truth. The check: the run's <calls>.fill_check.tsv, summed by group."""
import csv, gzip, sys
from collections import Counter
calls, truth, groups = sys.argv[1:4]
def load(path):
    out, samples = {}, None
    for line in gzip.open(path, "rt"):
        if line.startswith("##"): continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#"): samples = f[9:]; continue
        out[(f[0], int(f[1]))] = ([f[3]] + f[4].split(","), f[8].split(":"), dict(zip(samples, f[9:])))
    return out, samples
T, samples = load(truth); C, _ = load(calls)
if groups == "f2":
    group_of = {s: ("parents" if i < 2 else "F2 1x") for i, s in enumerate(samples)}
else:
    g = groups.split(","); group_of = {s: g[i // 10] for i, s in enumerate(samples)}
n = Counter()
for key, (al, _, cols) in T.items():
    if len(al) != 2 or any(len(a) != 1 for a in al) or key not in C: continue
    cal, keys, ccols = C[key]
    if cal[0] != al[0]: continue
    for s in samples:
        v = dict(zip(keys, ccols[s].split(":")))
        if v.get("DP") not in ("0", "."): continue
        if "." in v["GT"]: n[(group_of[s], "withheld")] += 1; continue
        t = sorted(al[int(x)] for x in cols[s].split("|"))
        c = sorted(cal[int(x)] for x in v["GT"].replace("|", "/").split("/"))
        n[(group_of[s], "filled")] += 1
        n[(group_of[s], "wrong")] += c != t
check = Counter()
for r in csv.DictReader(open(calls.replace(".vcf.gz", ".fill_check.tsv")), delimiter="\t"):
    gname = group_of[r["sample"]]
    check[(gname, "checked")] += int(r["sites_checked"]); check[(gname, "wrong")] += int(r["fills_wrong_there"])
print("group     true: filled  wrong  rate   | check: checked wrong  rate")
for gname in dict.fromkeys(group_of.values()):
    f, w = n[(gname, "filled")], n[(gname, "wrong")]
    c, cw = check[(gname, "checked")], check[(gname, "wrong")]
    print(f"{gname:>8}  {f:>12} {w:>6} {100*w/max(f,1):5.1f}%  | {c:>13} {cw:>5} {100*cw/max(c,1):5.1f}%")
