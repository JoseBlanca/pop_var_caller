"""Score ng_neighbour_pass_probe's output against the founder simulator's truth.vcf, at SNPs whose
record matches the truth's exactly (chrom, pos, ref, alt), by what the sample's reads said."""
import sys, collections
truth_vcf, probe = sys.argv[1], sys.argv[2]
truth, samples = {}, None
for line in open(truth_vcf):
    if line.startswith('##'): continue
    f = line.rstrip('\n').split('\t')
    if line.startswith('#'): samples = f[9:]; continue
    if len(f[3]) == 1 and len(f[4]) == 1:
        truth[(f[0], int(f[1]), f[3], f[4])] = {s: '/'.join(sorted(g.replace('|', '/').split('/'))) for s, g in zip(samples, f[9:])}
truth_positions = {(c, p) for c, p, _, _ in truth}
norm = lambda gt: '/'.join(sorted(gt.split('/'))) if gt not in ('.', './.') else './.'
count = collections.Counter()
for line in open(probe):
    f = line.rstrip('\n').split('\t')
    if f[0] == 'contig': continue
    contig, pos, ref, alt, sample, row, calling, outcome, passgt = f[0], int(f[1]), f[2], f[3], f[4], f[5], norm(f[6]), f[7], f[8]
    key = (contig, pos, ref, alt)
    if key in truth:
        t = truth[key][sample]
    elif len(ref) == 1 and len(alt) == 1 and (contig, pos) not in truth_positions:
        t = '0/0'   # a SNP record where the truth has nothing
    else:
        continue
    kind = ('nonref' if t != '0/0' else 'ref') + ('' if key in truth else '_no_truth_record')
    count[(row, kind, 'n')] += 1
    count[(row, kind, 'calling_wrong')] += calling != './.' and calling != t
    count[(row, kind, 'calling_nocall')] += calling == './.'
    if outcome == 'called':
        count[(row, kind, 'pass_wrong')] += norm(passgt) != t
    elif outcome == 'withheld':
        count[(row, kind, 'pass_withheld')] += 1
    else:
        count[(row, kind, 'pass_unchanged')] += 1
        count[(row, kind, 'pass_wrong')] += calling != './.' and calling != t
print('row\ttruth\tsample-sites\tcalling wrong\tcalling ./.\tpass wrong\tpass ./. (withheld)\tpass unchanged')
for row in ('informative', 'flat'):
    for kind in ('nonref', 'ref', 'ref_no_truth_record'):
        g = lambda k: count[(row, kind, k)]
        print(f"{row}\t{kind}\t{g('n')}\t{g('calling_wrong')}\t{g('calling_nocall')}\t{g('pass_wrong')}\t{g('pass_withheld')}\t{g('pass_unchanged')}")
