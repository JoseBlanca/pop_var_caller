"""Tabulate founder_grid.py's cells for the phase-four report.

usage: uv run python grid_table.py <work dir>

Per cell, one sample at one matched SNP record is one unit:
  with reads     wrong genotypes, calling -> pass
  without reads  sample-sites calling wrote ./.; how many the pass calls wrong, and how many of
                 those are true-reference genotypes called non-reference
  Q5 net         errors the pass removes where samples have reads, against the false
                 non-reference calls it adds anywhere (spec §11 Q5: a default must not add more
                 than it removes in any depth class)
And, over all cells, the pass's observed error rate by GQ bin, with and without reads.
"""
import collections, glob, json, os, sys

rows = []
for path in glob.glob(os.path.join(sys.argv[1], '*', 'score.json')):
    rows.append(json.load(open(path)))
rows.sort(key=lambda r: (r['founders'], r['inbreeding'], r['samples'], r['depth']))
print('| founders | F | samples | depth | with reads: wrong, calling → pass | without reads: sites, pass wrong (of them ref→non-ref) | Q5: errors removed / false non-ref added | pass s |')
print('|---|---|---|---|---|---|---|---|')
gq = collections.Counter()
for r in rows:
    c = collections.Counter(r['counts'])
    calling = c['reads_nonref_calling_wrong'] + c['reads_ref_calling_wrong']
    passing = c['reads_nonref_pass_wrong'] + c['reads_ref_pass_wrong']
    with_reads = c['reads_nonref_n'] + c['reads_ref_n']
    without = c['no_reads_nonref_n'] + c['no_reads_ref_n']
    without_wrong = c['no_reads_nonref_pass_wrong'] + c['no_reads_ref_pass_wrong']
    false_added = (c['reads_ref_pass_false_nonref'] + c['no_reads_ref_pass_false_nonref']
                   - c['reads_ref_calling_false_nonref'])
    removed = calling - passing
    seconds = r['timings'].get('call_on', 0) - r['timings'].get('call_off', 0)
    print(f"| {r['founders']} | {r['inbreeding']} | {r['samples']} | {r['depth']:g}× | "
          f"{calling} → {passing} of {with_reads:,} | {without:,}: {without_wrong} ({c['no_reads_ref_pass_false_nonref']}) | "
          f"{removed} / {false_added} | {seconds:.1f} |")
    gq.update(r['gq'])
print()
print('A GQ floor on calls made without reads (below it, ./.): per cell, Q5 net and what is kept.')
print()
print('| founders | F | samples | depth | floor | Q5: errors removed / false non-ref added | the same, a ./. at a true non-ref genotype counted wrong | no-read calls kept, wrong | no-read true non-ref called right |')
print('|---|---|---|---|---|---|---|---|---|')
for r in rows:
    c = collections.Counter(r['counts'])
    calling = c['reads_nonref_calling_wrong'] + c['reads_ref_calling_wrong']
    passing = c['reads_nonref_pass_wrong'] + c['reads_ref_pass_wrong']
    with_reads_false = c['reads_ref_pass_false_nonref'] - c['reads_ref_calling_false_nonref']
    for floor in (0, 10, 20):
        if floor == 0:
            false_nonref = c['no_reads_ref_pass_false_nonref']
            kept = c['no_reads_nonref_n'] + c['no_reads_ref_n'] - c['no_reads_nonref_pass_nocall'] - c['no_reads_ref_pass_nocall']
            wrong = c['no_reads_nonref_pass_wrong'] + c['no_reads_ref_pass_wrong']
            right = c['no_reads_right_nonref']
        else:
            false_nonref = c[f'no_reads_floor{floor}_false_nonref']
            kept = c[f'no_reads_floor{floor}_called']
            wrong = c[f'no_reads_floor{floor}_wrong']
            right = c[f'no_reads_floor{floor}_right_nonref']
        print(f"| {r['founders']} | {r['inbreeding']} | {r['samples']} | {r['depth']:g}× | {floor} | "
              f"{calling - passing} / {with_reads_false + false_nonref} | "
              f"{calling - passing + right} / {with_reads_false + false_nonref} | {kept:,}, {wrong} | {right:,} |")
print()
print('| GQ | with reads: wrong of called | without reads: wrong of called |')
print('|---|---|---|')
for low, high in [(0, 9), (10, 19), (20, 29), (30, 49), (50, 98), (99, 99)]:
    cells = []
    for row in ('reads', 'no_reads'):
        n, wrong = gq[f'{row}_gq{low}-{high}_n'], gq[f'{row}_gq{low}-{high}_wrong']
        cells.append(f'{wrong} of {n:,}' + (f' ({100 * wrong / n:.1f}%)' if n else ''))
    print(f'| {low}–{high} | {cells[0]} | {cells[1]} |')
