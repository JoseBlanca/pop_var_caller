"""Agreement with GATK on the tomato cohort, which has no truth set (plan §4.7).

Uses the whole-cohort calls the speed benchmark makes (pvc both ways, and GATK's joint calls), so it
runs no caller of its own; it needs `speed` to have run at all 62 accessions (2 in the slice).
"""

import os

from lib import agreement, arms
from lib.steps import Step, file_identity

from . import speed

TITLE = "Agreement with GATK — tomato cohort (no truth set)"
DESCRIPTION = """
Tomato (S. lycopersicum, SL4.0), all 62 accessions over 80 regions of 100 kb, about 3 reads a
position. There is no truth set, so this is agreement, not accuracy. Both callers' VCFs are
left-aligned, split and decomposed the same way (bcftools norm -m -any -a); the unit is one sample's
genotype at one variant where at least one caller gives that sample the allele. same_genotype: both
give it the same number of copies; other_genotype: both give it the allele, in different numbers;
only_pvc / only_GATK: one caller gives it the allele and the other none; no_call_*: one caller gives
it the allele and the other leaves the sample uncalled.
"""


def steps(ctx):
    size = 2 if ctx.slice else speed.SIZES[-1]
    _, regions = speed.regions_step(ctx)
    samples = list(speed.cohort(ctx))[:size]
    gatk = os.path.join(ctx.results_root, "competitors", f"speed/tomato_{size}", "gatk",
                        "calls.vcf.gz")
    out = []
    for arm in ("pvc_defaults", "pvc_fitted"):
        calls = os.path.join(ctx.pvc_root(), f"speed/tomato_{size}", arm, "calls.vcf.gz")
        base = dict(benchmark="tomato_agreement",
                    dataset="tomato cohort, SL4.0, 80 regions of 100 kb", samples=size,
                    depth="about 3x", region_bp=arms.BedLength(regions),
                    caller=arms.ARM_LABELS[arm], caller_version=ctx.version("pvc"))

        def action(d, calls=calls, base=base):
            for path in (calls, gatk):
                if not os.path.exists(path):
                    raise FileNotFoundError(f"{path} is missing: run --only speed first")
            counts = agreement.compare(calls, gatk, speed.REFERENCE, regions, samples)
            return {"rows": agreement.rows(base, counts, "pvc", "GATK")}

        out.append(Step(f"tomato_agreement/{arm}", "score", action,
                        lambda calls=calls: {"a": file_identity(calls), "b": file_identity(gatk),
                                             "version": 1},
                        attrs={"caller": "pvc", "home": "pvc"}))
    return out
