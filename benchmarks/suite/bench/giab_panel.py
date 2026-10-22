"""GIAB per-sample panel: SNP and indel accuracy, one human sample at a time (plan §4.1).

HG002, HG003 and HG004, each over its own 100 random GIAB confident regions, at seven depths. The
downsampled alignments were made outside the repository (`samtools view -s 42.<fraction>`, recorded
in their `.coverage.json`); the suite uses them as they are.
"""

import os

from lib import arms, giab_sources
from lib.steps import Step, file_identity

TITLE = "GIAB per-sample panel — SNPs and indels, one sample at a time"
DESCRIPTION = """
Human (GRCh38). GIAB HG002, HG003 and HG004, each called alone over its own 100 random regions of
GIAB's confident regions (HG002 572 kb, HG003 401 kb, HG004 453 kb). Depths 5× to 50× are
downsampled from the full reads (seed 42); "300x" is the full reads (about 300×). Truth: GIAB
v4.2.1 small variants, scored by RTG vcfeval inside the sample's confident regions. "genotype
match" needs the alleles and the genotype right; "allele match" only the alleles. Records are
scored as each caller writes them (its own filters), except freebayes, kept at QUAL ≥ 20. The
best_f1 rows are each caller's best F1 over a sweep of its QUAL threshold. pvc (fitted) fits its
parameters from that one sample's data with `estimate-parameters`.
"""

SAMPLES = ["HG002", "HG003", "HG004"]
DEPTHS = ["5x", "10x", "15x", "20x", "30x", "50x", "300x"]
FULL_DEPTH = {
    "HG002": "HG002_reads_selected_100_rg.cram",
    "HG003": "HG003_bench_azar_merged_100.sorted.m5.bam",
    "HG004": "HG004_bench_azar_merged_100.sorted.m5.bam",
}
ARMS = ["pvc_defaults", "pvc_fitted", "gatk", "freebayes", "deepvariant"]


def local_alignment(ctx, sample, depth):
    root = ctx.data("benchmarks", "giab", "per_sample", "bam", depth)
    if depth == "300x":
        return os.path.join(root, FULL_DEPTH[sample])
    return os.path.join(root, f"{sample}.{depth}.seed42.bam")


def steps(ctx):
    panel = ctx.data("benchmarks", "giab", "per_sample")
    samples, depths = (["HG002"], ["10x"]) if ctx.slice else (SAMPLES, DEPTHS)
    ref_steps, reference = giab_sources.grch38(ctx)
    out = list(ref_steps)
    more, catalog, sdf = arms.reference_steps(ctx, "GRCh38", reference)
    out.extend(more)

    def bed_of(s):
        # The data checkout's copy where it has one (the laptop), else the suite's identical copy.
        local = os.path.join(panel, "bed", f"{s}_bench_azar_merged_100.bed")
        return local if os.path.exists(local) else giab_sources.data_file(
            f"giab_panel/{s}_bench_azar_merged_100.bed")

    # The regions each sample is called and scored over: its whole bed, or its first five regions
    # for the slice. The scored regions are GIAB's confident regions, so they double as the region list.
    regions_dir = os.path.join(ctx.results_root, "inputs", "giab_panel", "regions")

    def make_regions(step_dir):
        for s in samples:
            with open(bed_of(s)) as fh:
                lines = [l for l in fh if l.strip()]
            if ctx.slice:
                lines = lines[:5]
            with open(os.path.join(step_dir, f"{s}.bed"), "w") as fh:
                fh.writelines(lines)
        return {}

    out.append(Step("giab_panel/regions", "input", make_regions,
                    {"beds": [file_identity(bed_of(s)) for s in samples],
                     "slice_regions": 5 if ctx.slice else "all"}))

    # Alignments and truth: the data checkout's files, or fetched from GIAB (lib/giab_sources.py).
    fractions = giab_sources.panel_fractions()
    aligned, truths = {}, {}
    for s in samples:
        if os.path.exists(local_alignment(ctx, s, depths[0])):
            for depth in depths:
                aligned[(s, depth)] = local_alignment(ctx, s, depth)
        else:
            step, full = giab_sources.fetched(ctx, "giab_panel", s, bed_of(s), reference)
            out.append(step)
            for depth in depths:
                if depth == "300x":
                    aligned[(s, depth)] = full
                    continue
                step, path = giab_sources.downsampled(ctx, "giab_panel", s, full, reference, depth,
                                                      fractions[(s, depth)])
                out.append(step)
                aligned[(s, depth)] = path
        local_truth = os.path.join(panel, "vcf", f"{s}_GRCh38_1_22_v4.2.1_benchmark.selected_100.vcf.gz")
        if os.path.exists(local_truth):
            truths[s] = local_truth
        else:
            more, whole, _ = giab_sources.whole_truth(ctx, s)
            out.extend(more)
            step, truths[s] = giab_sources.cut_truth(ctx, f"giab_panel_{s}", whole, bed_of(s))
            out.append(step)

    for s in samples:
        regions = os.path.join(regions_dir, f"{s}.bed")
        truth = truths[s]
        for depth in depths:
            # The region length is read when the row is written, after the regions step ran.
            base = dict(benchmark="giab_panel", dataset=f"GIAB {s}, own 100 confident regions",
                        samples=1, depth=depth, region_bp=arms.BedLength(regions))
            for arm in ARMS:
                step, callset = arms.arm_steps(ctx, arm, f"giab_panel/{s}/{depth}", reference,
                                               catalog, {s: aligned[(s, depth)]}, regions,
                                               base)
                out.append(step)
                out.append(arms.score_step(ctx, callset, base, truth, regions, sdf,
                                           call_sample=s, truth_sample=s))
    return out
