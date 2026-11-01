"""The neighbour-evidence pass measured on the founder simulator (spec §8.2, §8.3): one cohort a
cell, called by pvc with the pass off and on, scored against the simulator's truth.

usage: uv run python founder_grid.py <work dir> <cell> [<cell> ...]
       a cell is founders,inbreeding,samples,depth  e.g. 16,0.95,50,3

Each cell's files go in <work dir>/f<founders>_F<inbreeding>_n<samples>_d<depth>/; a cell whose
score.json exists is skipped. Needs, on PATH or by absolute path below: minimap2, samtools, and the
release builds of pvc and the simulator in the worktree.

Scoring: SNP records whose position and alleles match a truth record; one unit is one sample at
one record. Rows split by whether calling gave the sample a genotype there (it had reads) or wrote
`./.` (it had none). For the pass's calls, the observed error rate by GQ bin.
"""
import collections, gzip, json, os, subprocess, sys, time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))
PVC = os.path.join(ROOT, 'target', 'release', 'pvc')
SIMULATOR = os.path.join(ROOT, 'target', 'release', 'examples', 'ng_str_stutter_benchmark')
THREADS = 6
GENOME = 'genome_bp=2000000 contigs=2 segments_per_period=250 bundles=250'
GQ_BINS = [(0, 9), (10, 19), (20, 29), (30, 49), (50, 98), (99, 99)]


def run(command, log):
    with open(log, 'a') as handle:
        handle.write('$ ' + command + '\n')
        handle.flush()
        subprocess.run(command, shell=True, check=True, stdout=handle, stderr=subprocess.STDOUT,
                       executable='/bin/bash')


def read_vcf(path, with_quality=False):
    out, samples = {}, None
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('##'):
                continue
            f = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                samples = f[9:]
                continue
            if not (len(f[3]) == 1 and all(len(a) == 1 for a in f[4].split(','))):
                continue
            alleles = [f[3]] + f[4].split(',')
            keys = f[8].split(':')
            gq_at = keys.index('GQ') if 'GQ' in keys else None
            calls = {}
            for s, col in zip(samples, f[9:]):
                values = col.split(':')
                gt = values[0].replace('|', '/')
                genotype = None if '.' in gt else tuple(sorted(alleles[int(i)] for i in gt.split('/')))
                quality = None
                if with_quality and gq_at is not None and values[gq_at] != '.':
                    quality = int(values[gq_at])
                calls[s] = (genotype, quality)
            out[(f[0], int(f[1]), f[3])] = calls
    return out, samples


def score(cell_dir):
    truth, samples = read_vcf(os.path.join(cell_dir, 'sim', 'truth.vcf.gz'))
    calling, _ = read_vcf(os.path.join(cell_dir, 'calls_off.vcf.gz'))
    refined, _ = read_vcf(os.path.join(cell_dir, 'calls_on.vcf.gz'), with_quality=True)
    count = collections.Counter()
    gq = collections.Counter()
    for key, before in calling.items():
        if key not in truth:
            continue
        after = refined.get(key)
        for s in samples:
            t = truth[key][s][0]
            kind = 'nonref' if any(a != key[2] for a in t) else 'ref'
            row = 'reads' if before[s][0] is not None else 'no_reads'
            prefix = f'{row}_{kind}'
            count[prefix + '_n'] += 1
            b = before[s][0]
            count[prefix + '_calling_wrong'] += b is not None and b != t
            count[prefix + '_calling_false_nonref'] += kind == 'ref' and b is not None and b != t
            a, quality = after[s] if after is not None else ((key[2], key[2]), 99)
            count[prefix + '_pass_nocall'] += a is None
            count[prefix + '_pass_wrong'] += a is not None and a != t
            count[prefix + '_pass_false_nonref'] += kind == 'ref' and a is not None and a != t
            if a is not None and quality is not None:
                for low, high in GQ_BINS:
                    if low <= quality <= high:
                        gq[f'{row}_gq{low}-{high}_n'] += 1
                        gq[f'{row}_gq{low}-{high}_wrong'] += a != t
    return dict(count), dict(gq)


def cell(work, spec):
    founders, inbreeding, samples, depth = spec.split(',')
    name = f'f{founders}_F{inbreeding}_n{samples}_d{depth}'
    d = os.path.join(work, name)
    if os.path.exists(os.path.join(d, 'score.json')):
        return
    os.makedirs(d, exist_ok=True)
    log = os.path.join(d, 'driver.log')
    timings = {}
    began = time.time()
    if not os.path.exists(os.path.join(d, 'sim', 'truth.vcf.gz')):
        run(f'{SIMULATOR} {d}/sim samples={samples} depth={depth} founders={founders} '
            f'inbreeding={inbreeding} {GENOME} design=pcr_free threads={THREADS} seed=20261102', log)
        run(f'samtools faidx {d}/sim/reference.fa && bgzip -f {d}/sim/truth.vcf', log)
    timings['simulate'] = time.time() - began
    names = [f'sim{i:03d}' for i in range(int(samples))]
    os.makedirs(f'{d}/psps', exist_ok=True)
    catalog = f'{d}/catalog.parquet'
    began = time.time()
    if not os.path.exists(catalog):
        run(f'{PVC} repeat-catalog --reference {d}/sim/reference.fa --output {catalog} --force', log)
    timings['catalog'] = time.time() - began
    began = time.time()
    for s in names:
        if os.path.exists(f'{d}/psps/{s}.psp'):
            continue
        rg = f"'@RG\\tID:{s}\\tSM:{s}\\tLB:{s}\\tPL:ILLUMINA'"
        run(f'minimap2 -ax sr -t {THREADS} -R {rg} {d}/sim/reference.fa {d}/sim/{s}_R1.fastq.gz '
            f'{d}/sim/{s}_R2.fastq.gz 2>/dev/null | samtools sort -@ 2 -o {d}/{s}.bam - '
            f'&& samtools index {d}/{s}.bam '
            f'&& {PVC} generate-psps --reference {d}/sim/reference.fa --catalog {catalog} '
            f'--alignment {d}/{s}.bam --output-dir {d}/psps --force '
            f'&& rm -f {d}/{s}.bam {d}/{s}.bam.bai {d}/sim/{s}_R1.fastq.gz {d}/sim/{s}_R2.fastq.gz',
            log)
    timings['map_and_psps'] = time.time() - began
    psps = ' '.join(f'--psp {d}/psps/{s}.psp' for s in names)
    common = f'--reference {d}/sim/reference.fa --catalog {catalog} {psps} --threads {THREADS}'
    began = time.time()
    if not os.path.exists(f'{d}/params.toml'):
        run(f'{PVC} estimate-parameters {common} --output {d}/params.toml --log {d}/params.log --force', log)
    timings['estimate_parameters'] = time.time() - began
    for mode in ('off', 'on'):
        began = time.time()
        run(f'{PVC} call-from-psps {common} --parameters {d}/params.toml --neighbour-evidence {mode} '
            f'--output {d}/calls_{mode}.vcf.gz --log {d}/call_{mode}.log --force', log)
        timings[f'call_{mode}'] = time.time() - began
    count, gq = score(d)
    with open(log) as handle:
        tally = [line.strip() for line in handle if 'neighbour-evidence pass:' in line]
    with open(os.path.join(d, 'score.json'), 'w') as handle:
        json.dump({'cell': name, 'founders': int(founders), 'inbreeding': float(inbreeding),
                   'samples': int(samples), 'depth': float(depth), 'counts': count, 'gq': gq,
                   'timings': timings, 'pass_log': tally}, handle, indent=1)
    print(name, 'done', flush=True)


if __name__ == '__main__':
    work = os.path.abspath(sys.argv[1])
    for spec in sys.argv[2:]:
        cell(work, spec)
