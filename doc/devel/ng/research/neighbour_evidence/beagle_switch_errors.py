"""How well Beagle phased a cell: over the samples of its first leave-one-out run (out_0), the truly
heterozygous sites it writes heterozygous, and how often the phase flips between consecutive such
sites (a switch error; random phase flips half the time).
usage: uv run python beagle_switch_errors.py <cell dir> [...]"""
import gzip, sys
for cell in sys.argv[1:]:
    truth = {}
    for line in gzip.open(f'{cell}/sim/truth.vcf.gz', 'rt'):
        if line.startswith('##'): continue
        f = line.rstrip('\n').split('\t')
        if line.startswith('#'): names = f[9:]; continue
        truth[(f[0], int(f[1]))] = dict(zip(names, f[9:]))
    het = called_het = pairs = switches = 0
    rows = []
    for line in gzip.open(f'{cell}/beagle/out_0.vcf.gz', 'rt'):
        if line.startswith('##'): continue
        f = line.rstrip('\n').split('\t')
        if line.startswith('#'): bn = f[9:]; continue
        rows.append(((f[0], int(f[1])), dict(zip(bn, f[9:]))))
    for s in bn:
        prev = None; contig = None
        for key, gts in rows:
            if key[0] != contig: prev = None; contig = key[0]
            t = truth[key][s]; a, b = t.split('|')
            if a == b: continue
            het += 1
            x, y = gts[s].split('|')
            if x == y: continue
            called_het += 1
            o = (x == a)
            if prev is not None:
                pairs += 1; switches += o != prev
            prev = o
    print(f'{cell.split("/")[-1]}: truly het {het:,}; Beagle het there {100*called_het/het:.0f}%; switch between consecutive {100*switches/max(pairs,1):.0f}% of {pairs:,}')
