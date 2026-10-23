"""Inputs the suite makes itself, once: region lists and downsampled alignments."""

import json
import os
import subprocess

from . import callers
from .measure import run_commands
from .steps import Step, file_identity


def regions_step(ctx, name, source_bed, slice_count):
    """The regions a benchmark runs over: the whole source BED, or its first `slice_count` lines
    for the slice. Returns (step, path of the regions BED)."""
    path = os.path.join(ctx.results_root, "inputs", name, "regions", "regions.bed")

    def action(step_dir):
        with open(source_bed) as fh:
            lines = [l for l in fh if l.strip() and not l.startswith(("#", "track", "browser"))]
        with open(os.path.join(step_dir, "regions.bed"), "w") as fh:
            fh.writelines(lines[:slice_count] if ctx.slice else lines)
        return {}

    step = Step(f"{name}/regions", "input", action,
                {"bed": file_identity(source_bed),
                 "slice_regions": slice_count if ctx.slice else "all"})
    return step, path


def depth_ladder(ctx, name, reference, full, regions, depths):
    """Downsampled copies of full-depth alignments, made the way the GIAB panel's were.

    `full` maps sample to its full-depth alignment. Each sample's mean depth over `regions` is
    measured once; a depth `Nx` keeps the fraction N / mean of the reads with `samtools view -s
    42.<fraction>`. The depth "300x" or "full" means the full alignment, untouched.
    Returns (steps, alignment(sample, depth) -> path).
    """
    inputs = os.path.join(ctx.results_root, "inputs", name)
    depth_file = os.path.join(inputs, "mean_depth", "mean_depth.json")
    samples = list(full)
    steps = []

    def measure_depth(step_dir):
        means = {}
        for s in samples:
            cmd = (f"samtools depth -a -b {regions} --reference {reference} {full[s]} "
                   "| awk '{t+=$3; n++} END {printf \"%.3f\", t/n}'")
            means[s] = float(subprocess.run(cmd, shell=True, capture_output=True, text=True,
                                            check=True).stdout)
        with open(os.path.join(step_dir, "mean_depth.json"), "w") as fh:
            json.dump(means, fh)
        return {"meta": means}

    steps.append(Step(f"{name}/mean_depth", "input", measure_depth,
                      lambda: {"alignments": [file_identity(full[s]) for s in samples],
                               "regions": file_identity(regions)}))

    def alignment(s, depth):
        if depth in ("300x", "full"):
            return full[s]
        return os.path.join(inputs, "bam", depth, s, f"{s}.{depth}.seed42.bam")

    for depth in depths:
        if depth in ("300x", "full"):
            continue
        for s in samples:
            def downsample(step_dir, s=s, depth=depth):
                with open(depth_file) as fh:
                    mean = json.load(fh)[s]
                fraction = min(1.0, float(depth.rstrip("x")) / mean)
                m = run_commands([callers.downsample(reference, full[s], fraction,
                                                     alignment(s, depth))],
                                 os.path.join(step_dir, "logs"))
                return {"measurement": m, "meta": {"fraction": fraction, "source_mean": mean}}

            steps.append(Step(f"{name}/bam/{depth}/{s}", "input", downsample,
                              lambda s=s: {"alignment": file_identity(full[s]), "seed": 42,
                                           "mean_depth": file_identity(depth_file)}))
    return steps, alignment


def truth_step(ctx, name, truth):
    """An indexed copy of a truth VCF, for a truth file stored without its index.

    The data checkout is mounted read-only, so the index cannot be written beside the original.
    Returns (step, path of the copy).
    """
    path = os.path.join(ctx.results_root, "inputs", name, "truth", "truth.vcf.gz")

    def action(step_dir):
        m = run_commands([("index_truth", f"bcftools view --no-version {truth} -Oz -o {path} "
                                          f"&& bcftools index -t {path}")],
                         os.path.join(step_dir, "logs"))
        return {"measurement": m}

    return Step(f"{name}/truth", "input", action, {"truth": file_identity(truth)}), path
