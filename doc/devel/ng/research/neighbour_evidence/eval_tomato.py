"""Score the neighbour-evidence models on tomato1: a downsampled condition against the full-depth
calls of the same accessions.

usage: uv run --with numpy python eval_tomato.py <calls dir> <test name> <deep names file> [K]

Reads <calls>/<name>.rows.tsv (pvc's per-sample likelihood rows, branch neighbour-hmm),
<calls>/<name>.vcf (only for the sample order) and <calls>/<name>.parameters.toml (each
sample's genome-wide F). Truth is the `full` condition's baseline call where its posterior is
at least 0.99.
"""

import re
import sys
import time
from collections import defaultdict

import numpy as np

import hmm


def sample_names(vcf):
    with open(vcf) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    raise SystemExit(f"no header in {vcf}")


def inbreeding(toml):
    text = open(toml).read()
    return {m.group(1): float(m.group(2)) for m in re.finditer(
        r'sample = "([^"]+)", inbreeding_coefficient = \{ value = ([0-9.e-]+)', text)}


def load(calls, name):
    names = sample_names(f"{calls}/{name}.vcf")
    f_by_name = inbreeding(f"{calls}/{name}.parameters.toml")
    records = []
    with open(f"{calls}/{name}.rows.tsv") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            contig, start, end, kind, alleles, depths = fields[:6]
            per_sample = []
            for cell in fields[8:]:
                if cell == ".":
                    per_sample.append(None)
                    continue
                ll, conc = cell.split("|")
                ll = np.array([float(x) for x in ll.split(";")])
                lik = np.exp(ll - ll.max())
                ind, ibd = hmm.prior_rows([float(x) for x in conc.split(";")])
                per_sample.append((lik, ind, ibd))
            records.append(dict(
                key=(int(contig), int(start), int(end)), kind=kind,
                alleles=alleles.split(","), depth=[int(x) for x in depths.split(",")],
                rows=per_sample))
    records.sort(key=lambda r: r["key"])
    return names, f_by_name, records


def genotype_strings(alleles, g):
    j, k = hmm.genotype_pairs(len(alleles))[g]
    return tuple(sorted((alleles[j], alleles[k])))


def allele_share(post, n_alleles):
    share = np.zeros(n_alleles)
    for g, (j, k) in enumerate(hmm.genotype_pairs(n_alleles)):
        share[j] += post[g] / 2
        share[k] += post[g] / 2
    return share


def baseline_shares(names, f_by_name, records):
    """Every sample's allele share at every record under pvc's own prior (None if set aside)."""
    shares = []
    for j, name in enumerate(names):
        fj = f_by_name.get(name, 0.0)
        col = []
        for r in records:
            row = r["rows"][j]
            if row is None:
                col.append(None)
            else:
                col.append(allele_share(hmm.posterior(row[0], (1 - fj) * row[1] + fj * row[2]),
                                        len(r["alleles"])))
        shares.append(col)
    return shares


def models_for_sample(names, f_by_name, records, shares, s, k_templates):
    """Posterior per record under each model; at a record where pvc set the sample aside the
    baseline is None (pvc writes ./.) and the neighbour models see a flat likelihood under a
    prior built from the other samples' shares."""
    f = f_by_name.get(names[s], 0.0)
    sites = []
    for i, r in enumerate(records):
        row = r["rows"][s]
        if row is None:
            A = len(r["alleles"])
            others = [shares[j][i] for j in range(len(names)) if j != s and shares[j][i] is not None]
            alpha = 1e-3 + 2 * np.sum(others, axis=0) if others else np.ones(A)
            ind, ibd = hmm.prior_rows(alpha)
            row = (np.ones(len(ind)), ind, ibd)
        sites.append(row)
    positions = [r["key"][0] * 1e10 + r["key"][1] for r in records]
    out = {"baseline": [None if r["rows"][s] is None else p
                        for r, p in zip(records, hmm.baseline(sites, f))]}
    t = time.time()
    (az, _, _), (f_fit, run) = hmm.fit_autozygosity(sites, positions, f)
    out["autozygosity"] = az
    out["_time_autozygosity"] = time.time() - t
    out["_autozygosity_fit"] = (f, f_fit, run)
    if k_templates and len(names) > 1:
        t = time.time()
        own = [allele_share(p, len(r["alleles"])) for p, r in zip(hmm.baseline(sites, f), records)]

        def agreement(j):
            return -sum(np.abs(a - b).sum() for a, b in zip(own, shares[j]) if b is not None)

        chosen = sorted((j for j in range(len(names)) if j != s), key=agreement, reverse=True)[:k_templates]
        tmpl = []
        for i, r in enumerate(records):
            A = len(r["alleles"])
            tmpl.append(np.stack([shares[j][i] if shares[j][i] is not None else np.full(A, 1 / A)
                                  for j in chosen]))
        best = None
        for segment in (3e4, 3e5):
            post, ll = hmm.haplotypes(sites, positions, tmpl, f, run, segment)
            if best is None or ll > best[1]:
                best = (post, ll, segment)
        out["haplotypes"] = best[0]
        out["_time_haplotypes"] = time.time() - t
    return out


def main():
    calls, test, deep_file = sys.argv[1:4]
    k_templates = int(sys.argv[4]) if len(sys.argv) > 4 else 16
    # A deep file of "name" lines scores one cohort call; "file_base name" lines with a test name
    # like single_deep3 score one single-sample call per accession, single_deep3_<file_base>.
    deep_lines = [line.split() for line in open(deep_file) if line.strip()]
    deep = [fields[-1] for fields in deep_lines]
    t_names, t_f, t_records = load(calls, "full")
    single = test.startswith("single_")
    if not single:
        names, f_by_name, records = load(calls, test)
        by_key = {r["key"]: n for n, r in enumerate(records)}
        shares = baseline_shares(names, f_by_name, records)
    tally = defaultdict(lambda: defaultdict(int))
    for fields, name in zip(deep_lines, deep):
        if single:
            names, f_by_name, records = load(calls, f"{test}_{fields[0]}")
            by_key = {r["key"]: n for n, r in enumerate(records)}
            shares = baseline_shares(names, f_by_name, records)
        ts, s = t_names.index(name), names.index(name)
        truth = {}
        tf = t_f.get(name, 0.0)
        for r in t_records:
            row = r["rows"][ts]
            if row is None:
                continue
            p = hmm.posterior(row[0], (1 - tf) * row[1] + tf * row[2])
            g = int(np.argmax(p))
            if p[g] >= 0.99:
                truth[r["key"]] = (r["kind"], r["alleles"], genotype_strings(r["alleles"], g))
        out = models_for_sample(names, f_by_name, records, shares, s, k_templates)
        fit = out["_autozygosity_fit"]
        print(f"# {name}: F genome-wide {fit[0]:.2f}, HMM {fit[1]:.2f}, run {fit[2]:.0f}; "
              f"seconds autozygosity {out['_time_autozygosity']:.1f} "
              f"haplotypes {out.get('_time_haplotypes', 0):.1f}; records {len(records)}", file=sys.stderr)
        for key, (kind, t_alleles, t_gt) in truth.items():
            ref = (t_alleles[0], t_alleles[0])
            nonref = t_gt != ref
            het = t_gt[0] != t_gt[1]
            i = by_key.get(key)
            depth = records[i]["depth"][s] if i is not None else -1
            if i is None:
                dclass = "no record"
            elif records[i]["rows"][s] is None:
                dclass = "set aside (./.)"
            else:
                dclass = "0 reads" if depth == 0 else ("1-3 reads" if depth <= 3 else "4+ reads")
            kclass = "STR" if kind != "G" else "SNP/indel"
            for model in ("baseline", "autozygosity", "haplotypes"):
                if model not in out:
                    continue
                if i is None:
                    called = ref
                elif out[model][i] is None:
                    called = None  # pvc's ./.
                else:
                    called = genotype_strings(records[i]["alleles"], int(np.argmax(out[model][i])))
                wrong = called != t_gt
                for group in ("all", dclass, kclass):
                    c = tally[(model, group)]
                    c["sites"] += 1
                    c["wrong"] += wrong
                    c["nonref"] += nonref
                    c["nonref_wrong"] += wrong and nonref
                    c["het"] += het
                    c["het_wrong"] += wrong and het
    print("model\tgroup\tsites\twrong\tnonref_wrong\tnonref\thet_wrong\thet")
    for (model, group), c in sorted(tally.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        print(f"{model}\t{group}\t{c['sites']}\t{c['wrong']}\t{c['nonref_wrong']}\t{c['nonref']}"
              f"\t{c['het_wrong']}\t{c['het']}")


if __name__ == "__main__":
    main()
