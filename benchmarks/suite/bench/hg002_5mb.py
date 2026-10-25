"""HG002 over 5 Mb: one human sample, accuracy and speed (plan §4.3).

One HG002 alignment over 1,000 GIAB-confident intervals (5.08 Mb), at its full depth (about 300×)
and at 30×, downsampled by the suite. Its time and memory rows are the "one human sample" speed
figures.
"""

import os

from lib import arms, giab_sources, inputs, readfloor

TITLE = "HG002 over 5 Mb — one human sample, accuracy, speed and memory"
DESCRIPTION = """
Human (GRCh38). GIAB HG002 over 1,000 GIAB-confident intervals (5.08 Mb), called alone at 30×
(downsampled by the suite, seed 42) and at its full depth ("full", about 300×). Truth: GIAB v4.2.1
small variants, scored by RTG vcfeval inside the intervals, as for the per-sample panel. The time
and memory rows are each caller's whole pipeline on this one sample at 6 threads: pvc both ways
(one step, and generate-psps + estimate-parameters + call-from-psps), GATK HaplotypeCaller with 6
pair-HMM threads, freebayes as 6 processes over blocks of the intervals.
"""

DEPTHS = ["30x", "full"]
ARMS = ["pvc_defaults", "pvc_fitted", "gatk", "freebayes", "deepvariant"]


def steps(ctx):
    data = ctx.data("benchmarks", "human_genome_bottle")
    full = {"HG002": os.path.join(data, "crams", "HG002_reads_selected_1000_rg.cram")}
    truth = os.path.join(data, "results", "HG002_bottle_reference_1000rg.vcf.gz")
    intervals = os.path.join(data, "HG002_bench_azar_merged_1000.bed")
    depths = ["30x"] if ctx.slice else DEPTHS
    ref_steps, GRCH38 = giab_sources.grch38(ctx)
    out = list(ref_steps)
    more, catalog, scoring = arms.reference_steps(ctx, "GRCh38", GRCH38)
    out.extend(more)
    # The data checkout's files where it has them, else fetched from GIAB (lib/giab_sources.py).
    if not os.path.exists(full["HG002"]):
        step, full["HG002"] = giab_sources.fetched(ctx, "hg002_5mb", "HG002", intervals, GRCH38)
        out.append(step)
    if os.path.exists(truth):
        step, truth = inputs.truth_step(ctx, "hg002_5mb", truth)
    else:
        more, whole, _ = giab_sources.whole_truth(ctx, "HG002")
        out.extend(more)
        step, truth = giab_sources.cut_truth(ctx, "hg002_5mb", whole, intervals)
    out.append(step)
    step, regions = inputs.regions_step(ctx, "hg002_5mb", intervals, 10)
    out.append(step)
    ladder, alignment = inputs.depth_ladder(ctx, "hg002_5mb", GRCH38, full, regions, depths)
    out.extend(ladder)
    for depth in depths:
        base = dict(benchmark="hg002_5mb", dataset="GIAB HG002, 1,000 confident intervals",
                    samples=1, depth=depth, region_bp=arms.BedLength(regions))
        for arm in ARMS:
            step, callset = arms.arm_steps(ctx, arm, f"hg002_5mb/{depth}", GRCH38, catalog,
                                           {"HG002": alignment("HG002", depth)}, regions, base)
            out.append(step)
            out.append(arms.score_step(ctx, callset, base, truth, regions, scoring,
                                       call_sample="HG002", truth_sample="HG002"))
        out.extend(readfloor.floor_steps(ctx, f"hg002_5mb/{depth}", GRCH38,
                                         {"HG002": alignment("HG002", depth)}, regions, base))
        out.extend(readfloor.startup_steps(ctx, f"hg002_5mb/{depth}", GRCH38, catalog,
                                           {"HG002": alignment("HG002", depth)}, regions, base))
    return out
