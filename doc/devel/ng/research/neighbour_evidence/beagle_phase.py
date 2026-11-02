"""Beagle 5.5 as a measuring instrument for phasing_ceiling.py: phase a founder_grid.py cell's calls
once per target sample, the target left out, so its own reads never reach the haplotypes its
templates are given.

usage: uv run python beagle_phase.py prepare <cell dir> [<cell dir> ...]
       then, in the dev container (Java), run the commands written to <cell>/beagle/commands.sh

prepare writes <cell>/beagle/input.vcf — the biallelic SNP PASS records of calls_off.vcf.gz that
match a truth record (phasing_ceiling.py's records), GT only, ./. where calling wrote it — and one
Beagle command per sample, `excludesamples=` naming it, output <cell>/beagle/out_<index>.vcf.gz.
Beagle runs with its defaults (it estimates its population size and error rate) and imputes the
./. it is given, so every template haplotype is complete.
"""
import gzip, os, sys

JAR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..', 'tmp',
                                   'beagle', 'beagle.jar'))


def prepare(cell):
    out_dir = os.path.join(cell, 'beagle')
    os.makedirs(out_dir, exist_ok=True)
    truth = set()
    with gzip.open(os.path.join(cell, 'sim', 'truth.vcf.gz'), 'rt') as fh:
        for line in fh:
            if not line.startswith('#'):
                f = line.split('\t', 5)
                if len(f[3]) == 1 and len(f[4]) == 1:
                    truth.add((f[0], int(f[1]), f[3], f[4]))
    contigs, samples, rows = [], None, []
    with gzip.open(os.path.join(cell, 'calls_off.vcf.gz'), 'rt') as fh:
        for line in fh:
            if line.startswith('##contig'):
                contigs.append(line)
                continue
            if line.startswith('##'):
                continue
            f = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                samples = f[9:]
                continue
            if f[6] != 'PASS' or (f[0], int(f[1]), f[3], f[4]) not in truth:
                continue
            gts = [c.split(':', 1)[0] for c in f[9:]]
            rows.append('\t'.join([f[0], f[1], f'{f[0]}_{f[1]}', f[3], f[4], '.', 'PASS', '.', 'GT']
                                  + gts) + '\n')
    with open(os.path.join(out_dir, 'input.vcf'), 'w') as out:
        out.write('##fileformat=VCFv4.2\n')
        out.writelines(contigs)
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
                            + samples) + '\n')
        out.writelines(rows)
    with open(os.path.join(out_dir, 'commands.sh'), 'w') as commands:
        for index, sample in enumerate(samples):
            exclude = os.path.join(out_dir, f'exclude_{index}.txt')
            with open(exclude, 'w') as fh:
                fh.write(sample + '\n')
            commands.write(f'[ -e {out_dir}/out_{index}.vcf.gz ] || java -Xmx1g -jar {JAR} '
                           f'gt={out_dir}/input.vcf excludesamples={exclude} '
                           f'out={out_dir}/out_{index} nthreads=1 > /dev/null\n')
    print(f'{cell}: {len(rows):,} records, {len(samples)} samples')


if __name__ == '__main__':
    if sys.argv[1] != 'prepare':
        sys.exit(__doc__)
    for cell in sys.argv[2:]:
        prepare(cell)
