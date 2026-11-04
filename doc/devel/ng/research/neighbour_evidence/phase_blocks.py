"""The prototype check of neighbour_evidence_phase_output.md §7.1: block rules on phaser.py's dumped
rounds (with --extra-draws), scored against the simulator's truth.

usage: uv run --with numpy python phase_blocks.py <dump dir> <cell dir> [lean ...]

The dumps come from phaser.py --rounds 5 --haplotypes 16 --one-per-sample --extra-draws 8 --dump DIR;
EXTRA=n uses only the first n extra draws of each round.

Per sample and contig: the running consensus of rounds 2-5 (§2.1); a record is a candidate when it
is truly heterozygous (the proxy for "written het") and the consensus leans by >= lean; consecutive
candidates are linked unless farther apart than 300 kb (the prototype's S). A link rule then keeps or
breaks each link; blocks of one record are written unphased.

Rules (each round's paths read unaligned, by their own relative orientation of the two records):
  agree      every round's path agrees with the consensus's relative orientation
  agree+lean the same, and every round's path leans >= 0.5 at both records (the adopted rule)
  draws      all paths agree: the round's path and its extra draws, every round
  draws+lean the same, and every path leans at both records
  share>=x   at least x of all paths agree (x = 0.97, 0.9)
Calibration: links binned by the share of all paths that disagree; observed flips per bin.
"""
import os, pickle, sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phasing_ceiling as pc  # noqa: E402

SEGMENT = 300_000


FLIP_COST = 2.0


def load(dump, cell, r):
    with open(os.path.join(dump, f'{cell}.round{r}.pickle'), 'rb') as fh:
        return pickle.load(fh)


def lean(a, b):
    return np.abs(a - b) >= 0.5


def align_states(ref, other):
    """Spec §2.1 steps 1-3. Returns per-record swap flags (bool array)."""
    rf, rs = ref
    f, s = other
    both = np.flatnonzero(lean(rf, rs) & lean(f, s))
    flip = np.zeros(len(f), dtype=bool)
    if len(both) == 0:
        return flip
    agree = (rf[both] > rs[both]) == (f[both] > s[both])
    cost = np.array([0.0 if agree[0] else 1.0, 1.0 if agree[0] else 0.0])
    back = []
    for a in agree[1:]:
        stay, change = cost, cost[::-1] + FLIP_COST
        choice = (change < stay).astype(int)
        back.append(choice)
        cost = np.minimum(stay, change) + np.array([0.0 if a else 1.0, 1.0 if a else 0.0])
    state = 0 if cost[0] <= cost[1] else 1
    states = [state]
    for choice in reversed(back):
        if choice[state]:
            state = 1 - state
        states.append(state)
    states.reverse()
    current = states[0]
    pos = 0
    for rec in range(len(f)):
        while pos < len(both) and both[pos] <= rec:
            current = states[pos]
            pos += 1
        flip[rec] = current == 1
    return flip


def apply(other, flip):
    f, s = other
    return np.where(flip, s, f), np.where(flip, f, s)

RULES = ['agree', 'agree+lean', 'draws', 'draws+lean', 'share>=0.97', 'share>=0.9']
BINS = [(0.0, 0.0), (0.0, 0.03), (0.03, 0.1), (0.1, 0.25), (0.25, 1.01)]


def deep_count(name):
    depth = name.split('_d', 1)[1]
    return int(depth.split('x')[1]) if '+' in depth else 0


def n50(sizes):
    sizes = sorted(sizes, reverse=True)
    total, acc = sum(sizes), 0
    for z in sizes:
        acc += z
        if acc >= total / 2:
            return z
    return 0


def main():
    dump, cell = sys.argv[1], sys.argv[2]
    leans = [float(x) for x in sys.argv[3:]] or [0.5]
    name = os.path.basename(os.path.normpath(cell))
    deep = deep_count(name)
    _, records, _ = pc.load(cell)
    called = {(r['contig'], r['pos']): r['called'] for r in records}
    data = {k: load(dump, name, k) for k in range(2, 6)}
    contigs = data[5]['contigs']
    samples = list(data[5]['haplotypes'].keys())
    print(f'## {name}  ({len(samples)} samples, {deep} deep)')
    for lean_min in leans:
        stats = {rule: dict(het=0, phased=0, links=0, flips=0, sizes=[], spans=[],
                            groups=Counter(), group_flips=Counter()) for rule in RULES}
        calib = Counter()
        calib_flips = Counter()
        for s in samples:
            group = 'deep' if s < deep else 'shallow'
            for ci, contig in enumerate(contigs):
                rounds = [tuple(np.asarray(x, dtype=float) for x in data[k]['haplotypes'][s][contig])
                          for k in range(2, 6)]
                paths = []
                for k in range(2, 6):
                    paths.append(rounds[k - 2])
                    for f, t in data[k]['extra_draws'][s][contig][:int(os.environ.get('EXTRA', '99'))]:
                        paths.append((np.asarray(f, dtype=float), np.asarray(t, dtype=float)))
                sums = [rounds[0][0].copy(), rounds[0][1].copy()]
                n = 1
                for other in rounds[1:]:
                    a = apply(other, align_states((sums[0] / n, sums[1] / n), other))
                    sums[0] += a[0]
                    sums[1] += a[1]
                    n += 1
                c0, c1 = sums[0] / n, sums[1] / n
                rows = data[5]['records'][ci]
                positions = np.array([p for p, _ in rows])
                truth = [t[s] for _, t in rows]
                het = np.array([a != b for a, b in truth])
                first_alt = np.array([a == 1 for a, b in truth])
                fill = np.array([called[(contig, p)][s] is None for p, _ in rows])
                cand = np.flatnonzero(het & (np.abs(c0 - c1) >= lean_min))
                cons = c0 > c1
                correct = cons == first_alt
                round_o = np.array([f > t for f, t in rounds])
                round_lean = np.array([np.abs(f - t) >= 0.5 for f, t in rounds])
                path_o = np.array([f > t for f, t in paths])
                path_lean = np.array([np.abs(f - t) >= 0.5 for f, t in paths])
                for rule in RULES:
                    stats[rule]['het'] += int(het.sum())
                # links between consecutive candidates
                keep = {rule: [] for rule in RULES}
                for j in range(1, len(cand)):
                    i, p = cand[j], cand[j - 1]
                    near = positions[i] - positions[p] <= SEGMENT
                    rel = cons[i] == cons[p]
                    r_agree = ((round_o[:, i] == round_o[:, p]) == rel)
                    r_lean = round_lean[:, i] & round_lean[:, p]
                    p_agree = ((path_o[:, i] == path_o[:, p]) == rel)
                    p_lean = path_lean[:, i] & path_lean[:, p]
                    share = p_agree.mean()
                    flip = correct[i] != correct[p]
                    if near:
                        disagree = 1 - share
                        for lo, hi in BINS:
                            if (lo == hi == 0 and disagree == 0) or (lo < disagree <= hi):
                                calib[(lo, hi)] += 1
                                calib_flips[(lo, hi)] += flip
                                break
                    keep['agree'].append(near and r_agree.all())
                    keep['agree+lean'].append(near and r_agree.all() and r_lean.all())
                    keep['draws'].append(near and p_agree.all())
                    keep['draws+lean'].append(near and p_agree.all() and p_lean.all())
                    keep['share>=0.97'].append(near and share >= 0.97)
                    keep['share>=0.9'].append(near and share >= 0.9)
                for rule in RULES:
                    st = stats[rule]
                    blocks, current = [], [cand[0]] if len(cand) else []
                    for j, kept in enumerate(keep[rule], start=1):
                        if kept:
                            current.append(cand[j])
                        else:
                            blocks.append(current)
                            current = [cand[j]]
                    if current:
                        blocks.append(current)
                    for b in blocks:
                        if len(b) < 2:
                            continue  # a block of one record is written unphased
                        st['phased'] += len(b)
                        st['sizes'].append(len(b))
                        st['spans'].append(int(positions[b[-1]] - positions[b[0]]))
                        for p, i in zip(b, b[1:]):
                            flip = correct[i] != correct[p]
                            st['links'] += 1
                            st['flips'] += flip
                            g = group if deep else 'all'
                            if fill[i] or fill[p]:
                                g += ', touches a fill'
                            st['groups'][g] += 1
                            st['group_flips'][g] += flip
        print(f'\n### lean >= {lean_min}\n')
        print('| link rule | true hets phased | flips inside blocks | N50 records | N50 kb | phased in blocks >= 10 | by group: flips of links |')
        print('|---|---|---|---|---|---|---|')
        for rule in RULES:
            st = stats[rule]
            big = sum(z for z in st['sizes'] if z >= 10)
            groups = '; '.join(f'{g} {st["group_flips"][g]} of {st["groups"][g]}'
                               for g in sorted(st['groups']))
            print(f'| {rule} | {st["phased"]:,} of {st["het"]:,} ({100 * st["phased"] / st["het"]:.0f}%) | '
                  f'{st["flips"]:,} of {st["links"]:,} ({100 * st["flips"] / max(st["links"], 1):.2f}%) | '
                  f'{n50(st["sizes"])} | {n50(st["spans"]) / 1000:.0f} | '
                  f'{100 * big / max(st["phased"], 1):.0f}% | {groups} |')
        print('\nCalibration — links (within 300 kb, before any rule) by the share of all paths that '
              'disagree with the consensus:\n')
        print('| paths disagreeing | links | flips |')
        print('|---|---|---|')
        for lo, hi in BINS:
            label = 'none' if lo == hi == 0 else f'({lo}, {min(hi, 1.0)}]'
            k = calib[(lo, hi)]
            print(f'| {label} | {k:,} | {calib_flips[(lo, hi)]:,} ({100 * calib_flips[(lo, hi)] / max(k, 1):.2f}%) |')
    sys.stdout.flush()


if __name__ == '__main__':
    main()
