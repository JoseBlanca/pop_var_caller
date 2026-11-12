"""GIAB HG002, the whole genome (chromosomes 1 to 22): accuracy and speed on a whole human sample.

The same design as `giab_chr20`, over all 22 autosomes instead of one chromosome. The alignments are
GIAB's HG002 300× novoalign reads, already downsampled to each depth with seed 42 and kept on the
machine's data disk (`--giab-downsampled`, a path with `{sample}` and `{depth}` in it; rick's are
`/media/giab/coverage_alignments/`, made from the same 300× CRAMs `--giab-alignments` names, by
`samtools view -s 42.<fraction>`, the fraction set from each file's measured mean depth). The suite
does not make them: drawing four depths from a 312 GB CRAM takes most of a day, and the files exist.
The truth is GIAB's HG002 v4.2.1 small-variant benchmark, already whole-genome.

The callers run over the whole of each autosome, as a user would; the scoring is restricted to
GIAB's confident regions. It takes days, so a full run leaves it out: run it with
`--only giab_genome`.
"""

import os
import shlex

from lib import arms, giab_sources, readfloor
from lib.measure import run_commands
from lib.steps import Step, file_identity

TITLE = "GIAB HG002, whole genome — SNPs and indels over chromosomes 1 to 22"
DESCRIPTION = """
Human (GRCh38). GIAB HG002, chromosomes 1 to 22 (2.88 Gb), from GIAB's 300× Illumina alignments
(novoalign, 2×148 bp), downsampled to 5×, 10×, 30× and 50× with seed 42. Each caller runs over
the whole of every autosome. Truth: GIAB v4.2.1 small variants, scored by RTG vcfeval and pvc_eval
inside GIAB's confident regions. The time and memory rows are one sample, a whole human genome.
pvc runs with --call-satellites (off by default), as in the other GIAB benchmarks.
"""

IN_FULL_RUN = False
SAMPLE = "HG002"
CONTIGS = [f"chr{i}" for i in range(1, 23)]
SLICE = ("chr20", 10_000_000, 10_500_000)
# 50× stands where giab_chr20 has its source's own ~60×; the full 300× is left out, as GATK alone
# would take about a week on it.
DEPTHS = ["5x", "10x", "30x", "50x"]
ARMS = ["pvc_defaults", "pvc_fitted", "gatk", "freebayes", "deepvariant"]
q = shlex.quote


def steps(ctx):
    if not ctx.giab_downsampled:
        raise SystemExit("giab_genome needs the downsampled whole-genome alignments on this "
                         "machine: --giab-downsampled (a machine profile's BENCH_GIAB_DOWNSAMPLED)")
    truth_steps, whole_truth, whole_bed = giab_sources.whole_truth(ctx, SAMPLE)
    depths = ["10x"] if ctx.slice else DEPTHS
    inputs_dir = os.path.join(ctx.results_root, "inputs", "giab_genome")
    ref_steps, GRCH38 = giab_sources.grch38(ctx)
    out = list(ref_steps) + list(truth_steps)
    more, catalog, scoring = arms.reference_steps(ctx, "GRCh38", GRCH38)
    out.extend(more)

    def aligned(depth):
        return ctx.giab_downsampled.format(sample=SAMPLE, depth=depth)

    # ---- the regions the callers run over, and the truth and confident regions cut to them ----
    calling = os.path.join(inputs_dir, "truth", "calling.bed")
    confident = os.path.join(inputs_dir, "truth", "confident.bed")
    truth = os.path.join(inputs_dir, "truth", "truth.vcf.gz")

    def cut_truth(d):
        with open(GRCH38 + ".fai") as fh:
            lengths = {f[0]: int(f[1]) for f in (l.split() for l in fh)}
        if ctx.slice:
            ground = {SLICE[0]: (SLICE[1], SLICE[2])}
        else:
            ground = {c: (0, lengths[c]) for c in CONTIGS}
        with open(calling, "w") as fh:
            for c, (s, e) in ground.items():
                fh.write(f"{c}\t{s}\t{e}\n")
        with open(whole_bed) as src, open(confident, "w") as dst:
            for line in src:
                f = line.split()
                if f[0] not in ground:
                    continue
                s, e = max(int(f[1]), ground[f[0]][0]), min(int(f[2]), ground[f[0]][1])
                if s < e:
                    dst.write(f"{f[0]}\t{s}\t{e}\n")
        m = run_commands([("truth", f"bcftools view --no-version -R {q(calling)} "
                                    f"{q(whole_truth)} -Oz -o {q(truth)} "
                                    f"&& bcftools index -t {q(truth)}")],
                         os.path.join(d, "logs"))
        return {"measurement": m}

    out.append(Step("giab_genome/truth", "input", cut_truth,
                    {"truth": file_identity(whole_truth), "bed": file_identity(whole_bed),
                     "contigs": [SLICE[0]] if ctx.slice else CONTIGS,
                     "slice": SLICE if ctx.slice else None}))

    for depth in depths:
        base = dict(benchmark="giab_genome",
                    dataset="GIAB HG002, chromosomes 1 to 22" + (" (slice, 500 kb)" if ctx.slice else ""),
                    samples=1, depth=depth, region_bp=arms.BedLength(calling))
        alignments = {SAMPLE: aligned(depth)}
        for arm in ARMS:
            step, callset = arms.arm_steps(ctx, arm, f"giab_genome/{depth}", GRCH38, catalog,
                                           alignments, calling, base,
                                           call_satellites=arms.GIAB_CALL_SATELLITES)
            out.append(step)
            out.append(arms.score_step(ctx, callset, base, truth, confident, scoring,
                                       call_sample=SAMPLE, truth_sample=SAMPLE))
        out.extend(readfloor.floor_steps(ctx, f"giab_genome/{depth}", GRCH38, alignments,
                                         calling, base))
        out.extend(readfloor.startup_steps(ctx, f"giab_genome/{depth}", GRCH38, catalog,
                                           alignments, calling, base))
    return out
