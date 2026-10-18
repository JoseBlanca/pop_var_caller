"""Speed and peak memory of each caller's whole pipeline on the tomato cohort (plan §4.6).

63 tomato accessions at about 3 reads a position over 80 regions of 100 kb (8 Mb of SL4.0), called
at 1, 8 and 63 samples so the table shows how each caller grows with the cohort. The one-human-sample
figures come from the HG002 5 Mb benchmark (`hg002_5mb`), whose arms are timed the same way.

Each caller runs as its documentation tells a user to (plan §3.4), at 6 threads: one process per
sample for the one-sample stages (pvc generate-psps, GATK HaplotypeCaller GVCFs), one threaded
process for the cohort stages, freebayes as 6 processes over blocks of the regions.
"""

import os
import subprocess

from lib import arms
from lib.steps import Step, file_identity

TITLE = "Speed and peak memory — tomato cohort"
DESCRIPTION = """
Tomato (S. lycopersicum, SL4.0). 63 accessions at about 3 reads a position, over 80 regions of
100 kb (8 Mb); the first 1, 8 and all 63 accessions in name order. Each caller's whole pipeline at
6 threads, inside a Linux virtual machine given 6 cores and 32 GB on the machine named below:
wall_seconds is elapsed time, cpu_seconds user plus system time summed over every process,
peak_rss_mb_all_processes the largest total resident memory of the processes running at once (the
real cost of a stage that runs one process per sample), peak_rss_mb_largest_process the largest
single process. pvc is timed both ways it can be run: call-from-alignments in one step, and
generate-psps, estimate-parameters and call-from-psps ("whole pipeline" is their sum). The repeat
catalogue (repeat-catalog), built once per reference, is timed separately under "repeat_catalog".
"""

SIZES = [1, 8, 63]
ARMS = ["pvc_defaults", "pvc_fitted", "gatk", "freebayes"]
REFERENCE = os.path.expanduser("~/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa")


def cohort(ctx):
    """Every tomato CRAM, by its read-group sample name, in name order of the files."""
    crams = ctx.data("benchmarks", "tomato1", "crams")
    out = {}
    for f in sorted(os.listdir(crams)):
        if not f.endswith(".cram"):
            continue
        path = os.path.join(crams, f)
        header = subprocess.run(["samtools", "view", "-H", path], capture_output=True, text=True,
                                check=True).stdout
        sample = next(field[3:] for line in header.splitlines() if line.startswith("@RG")
                      for field in line.split("\t") if field.startswith("SM:"))
        out[sample] = path
    return out


def regions_step(ctx):
    src = ctx.data("benchmarks", "tomato1", "regions.bed")
    path = os.path.join(ctx.results_root, "inputs", "tomato", "regions", "regions.bed")

    def action(d):
        with open(src) as fh:
            lines = [l for l in fh if l.strip() and not l.startswith(("#", "track"))]
        with open(path, "w") as fh:
            fh.writelines(lines[:1] if ctx.slice else lines)
        return {}

    return Step("tomato/regions", "input", action,
                {"bed": file_identity(src), "slice_regions": 1 if ctx.slice else "all"}), path


def steps(ctx, sizes=None):
    sizes = sizes or ([2] if ctx.slice else SIZES)
    out, catalog, _ = arms.reference_steps(ctx, "SL4.0", REFERENCE)
    out = [s for s in out if not s.name.startswith("refs/")]  # no truth here, so no vcfeval copy
    step, regions = regions_step(ctx)
    out.append(step)
    samples = cohort(ctx)
    names = list(samples)
    for size in sizes:
        chosen = {s: samples[s] for s in names[:size]}
        base = dict(benchmark="speed", dataset="tomato cohort, SL4.0, 80 regions of 100 kb",
                    samples=size, depth="about 3x", region_bp=arms.BedLength(regions))
        for arm in ARMS:
            step, _ = arms.arm_steps(ctx, arm, f"speed/tomato_{size}", REFERENCE, catalog, chosen,
                                     regions, base, joint=True)
            out.append(step)
    return out
