"""GIAB chr22, one sample at a time: the colleague's benchmark, run when pointed at its data
(plan §4.8).

Not part of a full run and not runnable on the laptop, whose disk does not hold the data. Run it
with `--only external_chr22 --external-config <file.json>`, where the file names everything:

    {
      "reference": "/path/GRCh38.fna",
      "samples": {
        "HG002": {
          "truth": "/path/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
          "confident": "/path/HG002_chr22_confident.bed",
          "alignments": {"5x": "/path/HG002.5x.bam", "10x": "...", "30x": "...", "300x": "..."}
        },
        "HG003": {...}, "HG004": {...}
      }
    }

Every path must be visible inside the dev container: under the project, under `~/genomes`, or under
the directory given as `--data-root`, which run_suite.sh mounts. The reference must have its `.fai`
and GATK's `.dict` beside it. The arms and the scoring are the per-sample panel's.
"""

import json
import os

from lib import arms
from lib.steps import SuiteError

TITLE = "GIAB chr22 — SNPs and indels, one sample at a time, whole chromosome"
DESCRIPTION = """
Human (GRCh38). GIAB HG002, HG003 and HG004, each called alone over chromosome 22 inside its GIAB
confident regions, at the depths the configuration names. Truth: GIAB v4.2.1 small variants,
scored by RTG vcfeval, as for the per-sample panel.
"""
IN_FULL_RUN = False
ARMS = ["pvc_defaults", "pvc_fitted", "gatk", "freebayes"]


def steps(ctx):
    path = ctx.external_config
    if not path:
        raise SuiteError("external_chr22 needs --external-config <file.json>; see "
                         "benchmarks/suite/bench/external_chr22.py for its contents")
    if not os.path.exists(path):
        raise SuiteError(f"external_chr22: no configuration file at {path}")
    with open(path) as fh:
        config = json.load(fh)
    reference = config["reference"]
    missing = [p for p in [reference, reference + ".fai"] if not os.path.exists(p)]
    for s, c in config["samples"].items():
        missing += [p for p in [c["truth"], c["confident"], *c["alignments"].values()]
                    if not os.path.exists(p)]
    if missing:
        raise SuiteError("external_chr22: these files are not visible in the container: "
                         + ", ".join(missing[:5]) + (" ..." if len(missing) > 5 else ""))
    out, catalog, sdf = arms.reference_steps(ctx, "chr22_reference", reference)
    for s, c in config["samples"].items():
        for depth, alignment in c["alignments"].items():
            base = dict(benchmark="external_chr22", dataset=f"GIAB {s}, chromosome 22 confident",
                        samples=1, depth=depth, region_bp=arms.BedLength(c["confident"]))
            for arm in ARMS:
                step, callset = arms.arm_steps(ctx, arm, f"external_chr22/{s}/{depth}", reference,
                                               catalog, {s: alignment}, c["confident"], base)
                out.append(step)
                out.append(arms.score_step(ctx, callset, base, c["truth"], c["confident"], sdf,
                                           call_sample=s, truth_sample=s))
    return out
