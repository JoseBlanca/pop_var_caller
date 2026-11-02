"""Two builds of pvc's neighbour-evidence pass on the same founder_grid.py cells: each re-calls the
cell's pass-off VCF (refine-genotypes, the same pass as a run with it on), and both are scored
against the simulator's truth, deep and shallow samples apart.

usage: uv run python compare_templates.py <label>=<pvc> <label>=<pvc> -- <cell dir> [<cell dir> ...]

The deep samples of a mixed-depth cell (name ending d<shallow>+<deep>x<count>) are the first
<count>; in a cell of one depth every sample is in one group. Scoring as founder_grid.py: SNP
records whose position and alleles match a truth record, one sample at one record a unit.
"""
import collections, os, re, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from founder_grid import read_vcf  # noqa: E402

THREADS = 6


def refine(pvc, cell, label):
    out = os.path.join(cell, f'refined_{label}.vcf.gz')
    if not os.path.exists(out):
        subprocess.run([pvc, 'refine-genotypes', '--vcf', os.path.join(cell, 'calls_off.vcf.gz'),
                        '--parameters', os.path.join(cell, 'params.toml'), '--output', out,
                        '--log', os.path.join(cell, f'refined_{label}.log'),
                        '--threads', str(THREADS), '--force'],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return out


def score(cell, refined_path, deep_count):
    truth, samples = read_vcf(os.path.join(cell, 'sim', 'truth.vcf.gz'))
    calling, _ = read_vcf(os.path.join(cell, 'calls_off.vcf.gz'))
    refined, _ = read_vcf(refined_path)
    group = {s: ('deep' if i < deep_count else 'shallow') for i, s in enumerate(samples)}
    count = collections.Counter()
    for key, before in calling.items():
        if key not in truth:
            continue
        after = refined.get(key)
        for s in samples:
            t = truth[key][s][0]
            g = group[s]
            b = before[s][0]
            a = after[s][0] if after is not None else (key[2], key[2])
            if b is not None:
                count[g, 'reads'] += 1
                count[g, 'calling_wrong'] += b != t
                count[g, 'pass_wrong'] += a is not None and a != t
            else:
                count[g, 'missing'] += 1
                count[g, 'filled'] += a is not None
                count[g, 'filled_wrong'] += a is not None and a != t
    return count


def main():
    split = sys.argv.index('--')
    builds = [arg.split('=', 1) for arg in sys.argv[1:split]]
    print('| cell | samples | build | with reads: wrong, calling → pass | missing | filled | of them wrong | missing after |')
    print('|---|---|---|---|---|---|---|---|')
    for cell in sys.argv[split + 1:]:
        name = os.path.basename(os.path.normpath(cell))
        mixed = re.search(r'\+\d+(?:\.\d+)?x(\d+)$', name)
        deep_count = int(mixed.group(1)) if mixed else 0
        for label, pvc in builds:
            c = score(cell, refine(pvc, cell, label), deep_count)
            for g in ('deep', 'shallow'):
                if not c[g, 'reads'] + c[g, 'missing']:
                    continue
                share = f" ({100 * c[g, 'filled_wrong'] / c[g, 'filled']:.1f}%)" if c[g, 'filled'] else ''
                print(f"| {name} | {g if deep_count else 'all'} | {label} | "
                      f"{c[g, 'calling_wrong']:,} → {c[g, 'pass_wrong']:,} of {c[g, 'reads']:,} | "
                      f"{c[g, 'missing']:,} | {c[g, 'filled']:,} | {c[g, 'filled_wrong']:,}{share} | "
                      f"{c[g, 'missing'] - c[g, 'filled']:,} |", flush=True)


if __name__ == '__main__':
    main()
