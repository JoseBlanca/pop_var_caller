"""Does the fill check pvc reports estimate each sample's real fill error rate? On founder_grid.py
cells: each is called again with `--refine-genotypes on`, writing
<cell>/refined_<label>.fill_check.tsv beside its VCF, and each sample's reported share is set
against its fills' error rate measured on the simulator's truth.

usage: uv run python validate_fill_check.py <pvc> <label> <cell dir> [<cell dir> ...]

Per cell: the fills over all samples, true wrong and the check's estimate of wrong (each sample's
reported share times its fills, over samples with a share); per sample, in a second table, the
samples with 30 fills or more. Scoring as founder_grid.py: SNP records matching a truth record. The
check covers every record kind, the truth only SNPs, so the two are compared as rates.
"""
import csv, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from founder_grid import read_vcf, refine_call  # noqa: E402


def refine(pvc, cell, label):
    return refine_call(pvc, cell, label, 6)


def true_fills(cell, refined_path):
    truth, samples = read_vcf(os.path.join(cell, 'sim', 'truth.vcf.gz'))
    calling, _ = read_vcf(os.path.join(cell, 'calls_off.vcf.gz'))
    refined, _ = read_vcf(refined_path)
    fills = {s: [0, 0] for s in samples}
    for key, before in calling.items():
        after = refined.get(key)
        if key not in truth or after is None:
            continue
        for s in samples:
            if before[s][0] is None and after[s][0] is not None:
                fills[s][0] += 1
                fills[s][1] += after[s][0] != truth[key][s][0]
    return fills


def main():
    pvc, label, cells = sys.argv[1], sys.argv[2], sys.argv[3:]
    per_sample = []
    print('| cell | fills (SNPs) | truly wrong | samples with a share | their fills | truly wrong there | estimated wrong there |')
    print('|---|---|---|---|---|---|---|')
    for cell in cells:
        name = os.path.basename(os.path.normpath(cell))
        refined, table = refine(pvc, cell, label)
        fills = true_fills(cell, refined)
        reported = {r['sample']: r for r in csv.DictReader(open(table), delimiter='\t')}
        n = sum(f[0] for f in fills.values())
        wrong = sum(f[1] for f in fills.values())
        with_share = [s for s in fills if reported[s]['wrong_share'] != '.']
        n_share = sum(fills[s][0] for s in with_share)
        wrong_share = sum(fills[s][1] for s in with_share)
        estimate = sum(fills[s][0] * float(reported[s]['wrong_share']) for s in with_share)
        print(f'| {name} | {n:,} | {wrong:,} ({100 * wrong / max(n, 1):.2f}%) | {len(with_share)} of {len(fills)} | '
              f'{n_share:,} | {wrong_share:,} ({100 * wrong_share / max(n_share, 1):.2f}%) | '
              f'{estimate:.0f} ({100 * estimate / max(n_share, 1):.2f}%) |', flush=True)
        for s in with_share:
            if fills[s][0] >= 30:
                per_sample.append((name, s, fills[s][0], fills[s][1], float(reported[s]['wrong_share']),
                                   int(reported[s]['sites_checked'])))
    print()
    print('Samples with 30 fills or more: the true fill error rate against the reported share.')
    print()
    print('| reported share | samples | their fills | truly wrong |')
    print('|---|---|---|---|')
    for low, high in [(0, 0.005), (0.005, 0.01), (0.01, 0.02), (0.02, 0.05), (0.05, 0.1), (0.1, 1.01)]:
        rows = [r for r in per_sample if low <= r[4] < high]
        n = sum(r[2] for r in rows)
        w = sum(r[3] for r in rows)
        print(f'| {100 * low:g}–{100 * high:g}% | {len(rows)} | {n:,} | {w:,}'
              + (f' ({100 * w / n:.2f}%)' if n else '') + ' |')


if __name__ == '__main__':
    main()
