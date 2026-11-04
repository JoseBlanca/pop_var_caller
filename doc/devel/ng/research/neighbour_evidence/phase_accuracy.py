"""How well the phasing inside the pass phases: the haplotypes phaser.py --dump wrote, against the
founder simulator's true phase.

usage: uv run --with numpy python phase_accuracy.py <dump dir> <cell name> <rounds>

For each sample at each record its two haplotypes are each copy's probability of carrying the
alternative. A record is *phased* for a sample where the two differ by at least 0.5 — one copy
leans to the alternative, the other to the reference — and the truth is heterozygous there. Its
orientation is right where the copy leaning to the alternative is the truth's first copy carrying
it. Reported, as for Beagle (beagle_switch_errors.py): the share of truly heterozygous records
phased, and how often the orientation flips between consecutive phased records of a sample on a
contig (a switch error; a random phase flips half the time). Also the phased records where the truth
is homozygous.

Two estimates: the last round's haplotypes alone, and a consensus over rounds 2 to the last — each
round's orientation aligned to round 2's by the fewest flips (a two-state path over the records both
phase, a mismatch costing 1 and a flip 2), then the aligned haplotypes averaged.
"""
import os, pickle, sys

import numpy as np

FLIP_COST = 2.0


def load(dump, cell, round_):
    with open(os.path.join(dump, f'{cell}.round{round_}.pickle'), 'rb') as fh:
        return pickle.load(fh)


def leaning(first, second):
    return np.abs(first - second) >= 0.5


def align(reference, other):
    """Flip `other`'s copies stretch by stretch to agree with `reference`, by the fewest flips."""
    ref_first, ref_second = reference
    first, second = other
    both = np.flatnonzero(leaning(ref_first, ref_second) & leaning(first, second))
    flip = np.zeros(len(first), dtype=bool)
    if len(both) == 0:
        return first, second
    agree = (ref_first[both] > ref_second[both]) == (first[both] > second[both])
    # Two states: kept (0) or flipped (1); a site costs 1 where the state disagrees with it.
    cost = np.array([0.0 if agree[0] else 1.0, 1.0 if agree[0] else 0.0])
    back = []
    for a in agree[1:]:
        stay, change = cost, cost[::-1] + FLIP_COST
        choice = (change < stay).astype(int)
        best = np.minimum(stay, change)
        back.append(choice)
        cost = best + np.array([0.0 if a else 1.0, 1.0 if a else 0.0])
    state = int(np.argmin(cost))
    states = [state]
    for choice in reversed(back):
        if choice[state]:
            state = 1 - state
        states.append(state)
    states.reverse()
    # Each record takes the state of the last site both phase at or before it.
    current = states[0]
    position = 0
    for record in range(len(first)):
        while position < len(both) and both[position] <= record:
            current = states[position]
            position += 1
        flip[record] = current == 1
    return np.where(flip, second, first), np.where(flip, first, second)


def score(haplotypes, records):
    het = phased = wrong_hom = pairs = switches = 0
    for sample, by_contig in haplotypes.items():
        for contig_index, (first, second) in enumerate(by_contig):
            truth = records[contig_index]
            previous = None
            for i, (_, alleles) in enumerate(truth):
                a, b = alleles[sample]
                lean = abs(first[i] - second[i]) >= 0.5
                if a == b:
                    wrong_hom += lean
                    continue
                het += 1
                if not lean:
                    continue
                phased += 1
                orientation = (first[i] > second[i]) == (a == 1)
                if previous is not None:
                    pairs += 1
                    switches += orientation != previous
                previous = orientation
    return het, phased, wrong_hom, pairs, switches


def main():
    dump, cell, rounds = sys.argv[1], sys.argv[2], int(sys.argv[3])
    data = {r: load(dump, cell, r) for r in range(2, rounds + 1)}
    records = data[rounds]['records']
    as_lists = lambda hap: {s: [tuple(np.asarray(c) for c in hap[s][contig])
                                for contig in data[rounds]['contigs']] for s in hap}
    last = as_lists(data[rounds]['haplotypes'])
    rounds_haps = [as_lists(data[r]['haplotypes']) for r in range(2, rounds + 1)]
    consensus = {}
    for sample in last:
        per_contig = []
        for contig_index in range(len(records)):
            reference = rounds_haps[0][sample][contig_index]
            firsts, seconds = [reference[0]], [reference[1]]
            for other in rounds_haps[1:]:
                f, s = align(reference, other[sample][contig_index])
                firsts.append(f)
                seconds.append(s)
            per_contig.append((np.mean(firsts, axis=0), np.mean(seconds, axis=0)))
        consensus[sample] = per_contig
    print('| cell | estimate | truly het | phased | switch errors between consecutive phased | phased where truth homozygous |')
    print('|---|---|---|---|---|---|')
    for label, haps in (('last round', last), (f'consensus of rounds 2–{rounds}', consensus)):
        het, phased, wrong_hom, pairs, switches = score(haps, records)
        print(f'| {cell} | {label} | {het:,} | {phased:,} ({100 * phased / max(het, 1):.0f}%) | '
              f'{switches:,} of {pairs:,} ({100 * switches / max(pairs, 1):.1f}%) | {wrong_hom:,} |')


if __name__ == '__main__':
    main()
