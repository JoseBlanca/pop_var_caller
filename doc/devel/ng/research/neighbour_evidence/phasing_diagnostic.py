"""How well the subtraction in phasing_ceiling.py phases: over 20 samples, the truly heterozygous sites
it phases confidently (its two copies' alternative probabilities at least 0.6 apart), and how often
the phase flips between consecutive confident ones — a switch error; random phase flips half the
time. usage: uv run --with numpy python phasing_diagnostic.py <cell dir> [...]"""
import sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import phasing_ceiling as pc
for cell in sys.argv[1:]:
    pc.DATA[cell] = pc.load(cell)
    samples, records, _ = pc.DATA[cell]
    pc.CONTIGS[cell] = {}
    for r in records:
        pc.CONTIGS[cell].setdefault(r['contig'], []).append(r)
    pc.RANKS[id(records)] = [pc.ranking(records, t, len(samples)) for t in range(len(samples))]
    wanted = list(range(min(20, len(samples))))
    phased = pc.subtraction_haplotypes(cell, -1, wanted)
    het = confident = right_alleles = pairs = switches = 0
    for t in wanted:
        for contig, rows in pc.CONTIGS[cell].items():
            first, second = phased[t][contig]
            prev = None
            for i, r in enumerate(rows):
                a, b = r['truth'][t]
                if a == b:
                    prev_keep = prev
                    continue
                het += 1
                if abs(first[i] - second[i]) < 0.6:
                    continue
                confident += 1
                o = (first[i] > second[i]) == (a == 1)
                if prev is not None:
                    pairs += 1
                    switches += o != prev
                prev = o
    print(f'{cell}: truly het sites {het:,}; phased confidently {confident:,} ({100*confident/max(het,1):.0f}%); '
          f'switch between consecutive confident het sites {switches:,} of {pairs:,} ({100*switches/max(pairs,1):.0f}%)')
