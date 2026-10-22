"""GIAB trio: joint calling of HG002, HG003 and HG004, scored per member and for Mendelian
consistency (plan §4.2).

The three are called together over 100 random regions where all three are GIAB-confident. Only
full-depth alignments exist; the suite downsamples them to 5, 10 and 30× the way the per-sample
panel was made (`samtools view -s 42.<fraction>`, the fraction from each file's measured mean depth
over the regions), and keeps them under `inputs/`.
"""

import os

from lib import arms, hostjobs, inputs, mendelian
from lib.measure import run_commands
from lib.steps import Step, file_identity

TITLE = "GIAB trio — joint calling and Mendelian consistency"
DESCRIPTION = """
Human (GRCh38). GIAB HG002 (child), HG003 (father) and HG004 (mother) called together over 100
random regions where all three are GIAB-confident (452 kb). 5×, 10× and 30× are downsampled by the
suite from the full reads (seed 42); "300x" is the full reads. Each member is scored by RTG vcfeval
against its own GIAB v4.2.1 truth. The Mendelian rows count sites where all three have a called
genotype and at least one is not reference, and how many of those the child's genotype cannot
inherit; the "GIAB truth" rows are the same count on the three truth sets merged, the floor.
"""

MEMBERS = {"mother": "HG004", "father": "HG003", "child": "HG002"}
SAMPLES = ["HG002", "HG003", "HG004"]
DEPTHS = ["5x", "10x", "30x", "300x"]
ARMS = ["pvc_defaults", "pvc_fitted", "gatk", "freebayes", "deepvariant"]


def steps(ctx):
    trio = ctx.data("benchmarks", "giab", "mendelian")
    full = {s: os.path.join(trio, "bam", f"{s}.trio_common_random_100.sorted.bam") for s in SAMPLES}
    truth = {s: os.path.join(trio, "vcf", f"{s}_GRCh38_1_22_v4.2.1_benchmark.trio_common_random_100.vcf.gz")
             for s in SAMPLES}
    src_bed = os.path.join(trio, "bed", "trio_common_evaluable.random_100.seed42.bed")
    depths = ["10x"] if ctx.slice else DEPTHS
    out, catalog, sdf = arms.reference_steps(ctx, "GRCh38", arms.GRCH38)

    step, regions = inputs.regions_step(ctx, "giab_trio", src_bed, 5)
    out.append(step)
    ladder, alignment = inputs.depth_ladder(ctx, "giab_trio", arms.GRCH38, full, regions, depths)
    out.extend(ladder)

    for depth in depths:
        base = dict(benchmark="giab_trio", dataset="GIAB trio HG002/HG003/HG004, 100 common regions",
                    samples=3, depth=depth, region_bp=arms.BedLength(regions))
        bams = {s: alignment(s, depth) for s in SAMPLES}
        for arm in ARMS:
            step, callset = arms.arm_steps(ctx, arm, f"giab_trio/{depth}", arms.GRCH38, catalog,
                                           bams, regions, base, joint=True)
            out.append(step)
            for s in SAMPLES:
                member_base = dict(base, dataset=f"GIAB trio, member {s}, 100 common regions")
                out.append(arms.score_step(ctx, callset, member_base, truth[s], regions, sdf,
                                           call_sample=s, truth_sample=s))
            out.append(_mendelian_step(ctx, callset, base, regions))
        if depth == depths[0]:
            out.append(_truth_mendelian_step(ctx, base, truth, regions))
    return out


def _mendelian_step(ctx, callset, base, regions):
    row_base = dict(base, caller=arms.ARM_LABELS[callset.arm],
                    caller_version=ctx.version(arms.ARM_TOOL[callset.arm]))

    def action(step_dir):
        hostjobs.skip_with(ctx, callset)
        counts = mendelian.count(callset.path, regions, MEMBERS["mother"], MEMBERS["father"],
                                 MEMBERS["child"], min_qual=arms.MIN_QUAL.get(callset.arm))
        return {"rows": mendelian.rows(row_base, counts)}

    kind = arms.ARM_KIND[callset.arm]
    return Step(f"{callset.step_name}/mendelian", "score", action,
                lambda: {"calls": file_identity(callset.path), "regions": file_identity(regions),
                         "version": 1},
                attrs={"caller": arms.ARM_CALLER[callset.arm], "depth": base["depth"],
                       "home": "pvc" if kind == "pvc" else "competitor"})


def _truth_mendelian_step(ctx, base, truth, regions):
    """The same count on GIAB's three truth sets, merged with absent sites as reference."""
    row_base = dict(base, depth="-", caller="GIAB truth", caller_version="GIAB v4.2.1")

    def action(step_dir):
        merged = os.path.join(step_dir, "truth_merged.vcf.gz")
        files = " ".join(truth[s] for s in SAMPLES)
        run_commands([("merge", f"bcftools merge --no-version -0 {files} -Oz -o {merged} "
                                f"&& bcftools index -t {merged}")], os.path.join(step_dir, "logs"))
        counts = mendelian.count(merged, regions, MEMBERS["mother"], MEMBERS["father"],
                                 MEMBERS["child"])
        return {"rows": mendelian.rows(row_base, counts)}

    return Step("giab_trio/truth/mendelian", "score", action,
                lambda: {"truth": [file_identity(truth[s]) for s in SAMPLES],
                         "regions": file_identity(regions), "version": 1},
                attrs={"home": "competitor"})
