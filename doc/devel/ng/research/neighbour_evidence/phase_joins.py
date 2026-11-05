"""The join test of neighbour_evidence_phase_output.md §2.4 on two overlapping windows of each contig,
phased separately by phaser.py --region (window A = [0, 600 kb), window B = [400 kb, 1 Mb)).

usage: uv run --with numpy python phase_joins.py <window dumps dir> <cell dir> [k]

The window dumps are <dir>/win_0-600000 and <dir>/win_400000-1000000, from phaser.py --region.

A writes its records below 500 kb; B's left margin is [400, 500) kb, its core from 500 kb. Per sample
and contig: A's last written block; of its records in B's margin that B phases, the k nearest the
core boundary; the join holds if there are k, they lie in one block of B that reaches B's first
phased core record, and the relative orientation of A and B is the same at all k. Scored: the share
of joins that hold, and, where one holds, whether the link from A's last written phased record to
B's first phased core record flips against the truth.
"""
import os, sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from phase_blocks import load, align_states, apply  # noqa: E402

BOUNDARY, MARGIN_START = 500_000, 400_000
SEGMENT = 300_000


def window(dump, name, s, ci, rule):
    data = {k: load(dump, name, k) for k in range(2, 6)}
    contig = data[5]['contigs'][ci]
    rounds = [tuple(np.asarray(x, dtype=float) for x in data[k]['haplotypes'][s][contig]) for k in range(2, 6)]
    paths = []
    for k in range(2, 6):
        paths.append(rounds[k - 2])
        paths += [(np.asarray(f, float), np.asarray(t, float)) for f, t in data[k]['extra_draws'][s][contig]]
    if rule == 'agree+lean':
        paths = rounds
    sums = [rounds[0][0].copy(), rounds[0][1].copy()]
    n = 1
    for other in rounds[1:]:
        a = apply(other, align_states((sums[0] / n, sums[1] / n), other))
        sums[0] += a[0]; sums[1] += a[1]; n += 1
    c0, c1 = sums[0] / n, sums[1] / n
    rows = data[5]['records'][ci]
    pos = np.array([p for p, _ in rows])
    truth = [t[s] for _, t in rows]
    het = np.array([a != b for a, b in truth])
    cons = c0 > c1
    correct = cons == np.array([a == 1 for a, b in truth])
    cand = np.flatnonzero(het & (np.abs(c0 - c1) >= 0.5))
    po = np.array([f > t for f, t in paths])
    pl = np.array([np.abs(f - t) >= 0.5 for f, t in paths])
    blocks, cur = [], [cand[0]] if len(cand) else []
    for j in range(1, len(cand)):
        i, p = cand[j], cand[j - 1]
        ok = (pos[i] - pos[p] <= SEGMENT
              and (((po[:, i] == po[:, p]) == (cons[i] == cons[p])).all())
              and (pl[:, i] & pl[:, p]).all())
        if ok:
            cur.append(i)
        else:
            blocks.append(cur); cur = [i]
    if cur:
        blocks.append(cur)
    blocks = [b for b in blocks if len(b) >= 2]
    info = {}
    for bi, b in enumerate(blocks):
        for i in b:
            info[int(pos[i])] = (bi, bool(cons[i]), bool(correct[i]))
    links = [(p, i) for b in blocks for p, i in zip(b, b[1:])]
    inside = len(links), sum(correct[i] != correct[p] for p, i in links)
    return info, inside


def main():
    root, cell = sys.argv[1], sys.argv[2]
    k = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    name = os.path.basename(os.path.normpath(cell))
    da, db = f'{root}/win_0-600000', f'{root}/win_400000-1000000'
    meta = load(da, name, 5)
    samples, ncontigs = list(meta['haplotypes'].keys()), len(meta['contigs'])
    out = {}
    for rule in ('agree+lean', 'draws+lean'):
        c = Counter()
        for s in samples:
            for ci in range(ncontigs):
                a, ia = window(da, name, s, ci, rule)
                b, ib = window(db, name, s, ci, rule)
                c['links_inside'] += ia[0] + ib[0]
                c['flips_inside'] += ia[1] + ib[1]
                written = sorted(p for p in a if p < BOUNDARY)
                core_b = sorted(p for p in b if p >= BOUNDARY)
                if not written or not core_b:
                    c['nothing to join'] += 1
                    continue
                c['joins tried'] += 1
                last_block = a[written[-1]][0]
                shared = [p for p in written if a[p][0] == last_block and p >= MARGIN_START and p in b]
                shared = shared[-k:]
                if len(shared) < k:
                    c['too few shared'] += 1
                    continue
                first_core = core_b[0]
                b_block = b[first_core][0]
                if any(b[p][0] != b_block for p in shared):
                    c['not one block of B'] += 1
                    continue
                rel = {a[p][1] == b[p][1] for p in shared}
                if len(rel) != 1:
                    c['orientations disagree'] += 1
                    continue
                c['joined'] += 1
                same = rel.pop()
                # with B swapped to agree, B's correctness at a record flips when not `same`
                ok_a = a[written[-1]][2]
                cons_b, correct_b = b[first_core][1], b[first_core][2]
                truth_alt_b = cons_b if correct_b else not cons_b
                orient_b = cons_b if same else not cons_b
                truth_ok_b = orient_b == truth_alt_b
                c['join flips'] += ok_a != truth_ok_b
        out[rule] = c
    print(f'## {name}, k = {k}')
    print('| link rule | joins tried | joined | too few shared | not one block of B | orientations disagree | flips at joins | flips inside windows |')
    print('|---|---|---|---|---|---|---|---|')
    for rule, c in out.items():
        t = c['joins tried']
        print(f"| {rule} | {t} | {c['joined']} ({100 * c['joined'] / max(t, 1):.0f}%) | {c['too few shared']} | "
              f"{c['not one block of B']} | {c['orientations disagree']} | {c['join flips']} of {c['joined']} | "
              f"{c['flips_inside']} of {c['links_inside']} ({100 * c['flips_inside'] / max(c['links_inside'], 1):.2f}%) |")


if __name__ == '__main__':
    main()
