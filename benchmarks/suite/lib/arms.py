"""The arms every SNP/indel benchmark runs, built once here so every benchmark runs them alike.

An *arm* is one way of calling the same data:

| arm             | label in the table | what runs                                                      |
|-----------------|--------------------|----------------------------------------------------------------|
| `pvc_defaults`  | pvc (defaults)     | `call-from-alignments --defaults`                              |
| `pvc_fitted`    | pvc (fitted)       | `generate-psps`, `estimate-parameters`, `call-from-psps`       |
| `gatk`          | GATK               | HaplotypeCaller (one sample), or GVCFs, CombineGVCFs, GenotypeGVCFs |
| `freebayes`     | freebayes          | one process per block of regions, all samples together         |

Each arm is one step, and its scoring another. pvc's arms live under the pvc commit's directory and
rerun for every commit; GATK's and freebayes' live under `competitors/` and run once.
"""

import os
import shutil
from collections import namedtuple

from . import callers, hostjobs, pvc_eval, vcfeval
from .measure import check_alignments, merge, run_commands, warm
from .results import row, time_rows
from .steps import Step, file_identity

GRCH38 = os.path.expanduser(
    "~/genomes/h_sapiens/gca_grch38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna")

ARM_LABELS = {"pvc_defaults": "pvc (defaults)", "pvc_fitted": "pvc (fitted)", "gatk": "GATK",
              "freebayes": "freebayes", "hipstr": "HipSTR", "deepvariant": "DeepVariant"}
ARM_CALLER = {"pvc_defaults": "pvc", "pvc_fitted": "pvc", "gatk": "gatk", "freebayes": "freebayes",
              "hipstr": "hipstr", "deepvariant": "deepvariant"}
ARM_KIND = {"pvc_defaults": "pvc", "pvc_fitted": "pvc", "gatk": "competitor",
            "freebayes": "competitor", "hipstr": "competitor", "deepvariant": "competitor"}
ARM_TOOL = {"pvc_defaults": "pvc", "pvc_fitted": "pvc", "gatk": "gatk", "freebayes": "freebayes",
            "hipstr": "hipstr", "deepvariant": "deepvariant"}
# freebayes writes every candidate it considers; its documentation asks for a quality filter.
MIN_QUAL = {"freebayes": 20}


# What the scorers need of a reference: the FASTA, vcfeval's copy of it, and pvc_eval's STR regions
# over it (None without pvc_eval).
ScoringReference = namedtuple("ScoringReference", "fasta sdf str_regions")


def reference_steps(ctx, name, reference):
    """The pvc repeat catalogue (per pvc commit), and what the scorers need of the reference (once).

    Returns (steps, catalogue path, ScoringReference)."""
    catalog_path = os.path.join(ctx.pvc_root(), "catalogs", name, f"{name}.repeats.parquet")

    def catalog(step_dir):
        log = os.path.join(step_dir, "repeat_catalog.log")
        m = run_commands([callers.pvc_catalog(ctx, reference, catalog_path, log)],
                         os.path.join(step_dir, "logs"))
        base = dict(benchmark="repeat_catalog", dataset=f"{name} whole genome", samples="-",
                    depth="-", region_bp="whole genome", variant_class="-", caller="pvc",
                    caller_version=ctx.version("pvc"))
        return {"measurement": m, "rows": time_rows(base, m, "repeat-catalog")}

    sdf_path = os.path.join(ctx.results_root, "inputs", "refs", f"{name}.sdf", "sdf")

    def sdf(step_dir):
        m = run_commands([vcfeval.sdf_command(ctx, reference, sdf_path)],
                         os.path.join(step_dir, "logs"))
        return {"measurement": m}

    fai = reference + ".fai"
    steps = [
        Step(f"catalogs/{name}", "pvc", catalog, {"reference": file_identity(fai)},
             attrs={"caller": "pvc"}),
        Step(f"refs/{name}.sdf", "input", sdf, {"reference": file_identity(fai)}),
    ]
    str_regions = None
    if pvc_eval.available(ctx):
        # pvc_eval's STR finder over the whole reference, once: a scan at every scoring step would
        # cost about 40 s a chromosome each time. Remade when pvc_eval changes, since its finder
        # may have.
        str_regions = os.path.join(ctx.results_root, "inputs", "refs", f"{name}.pvc_eval_str",
                                   "str_regions.bed")

        def find_strs(step_dir):
            m = run_commands([pvc_eval.str_regions_command(ctx, reference, str_regions)],
                             os.path.join(step_dir, "logs"))
            return {"measurement": m}

        steps.append(Step(f"refs/{name}.pvc_eval_str", "input", find_strs,
                          {"reference": file_identity(fai), "pvc_eval": ctx.version("pvc_eval")}))
    return steps, catalog_path, ScoringReference(reference, sdf_path, str_regions)


def bed_length(path):
    total = 0
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 3 and not line.startswith(("#", "track", "browser")):
                total += int(parts[2]) - int(parts[1])
    return total


class BedLength:
    """A BED's total length, read when first turned into text (the BED may not exist yet)."""

    def __init__(self, path):
        self.path = path

    def __str__(self):
        return str(bed_length(self.path))

    __repr__ = __str__


class CallSet:
    """Where an arm's calls end up, and what the scoring needs to know about them."""

    def __init__(self, arm, step_name, path, sample=None):
        self.arm = arm
        self.step_name = step_name
        self.path = path
        self.sample = sample


def arm_steps(ctx, arm, prefix, reference, catalog, alignments_by_sample, regions, base,
              joint=False):
    """One arm's calling step on one configuration; returns (step, CallSet).

    `alignments_by_sample` maps sample name to its alignment file(s); several samples means a
    cohort, called jointly. `regions` is the BED the callers are restricted to. `base` is the row
    configuration (benchmark, dataset, samples, depth, region_bp).
    """
    name = f"{prefix}/{arm}"
    kind = ARM_KIND[arm]
    step_dir = os.path.join(ctx.pvc_root() if kind == "pvc" else
                            os.path.join(ctx.results_root, "competitors"), name)
    calls = os.path.join(step_dir, "calls.vcf.gz")
    samples = list(alignments_by_sample)
    files = {s: (v if isinstance(v, list) else [v]) for s, v in alignments_by_sample.items()}
    flat = [f for s in samples for f in files[s]]
    row_base = dict(base, caller=ARM_LABELS[arm], caller_version=ctx.version(ARM_TOOL[arm]))

    def fingerprint():
        return {
            "tool": ctx.version(ARM_TOOL[arm]) if kind != "pvc" else "pvc",
            "alignments": [file_identity(f) for f in flat],
            "regions": file_identity(regions),
            "reference": file_identity(reference + ".fai"),
            "threads": ctx.threads,
        }

    if arm == "pvc_defaults":
        def action(d):
            check_alignments(flat + [reference])
            warm(flat + [reference])  # every tool reads from memory (measure.warm)
            log = os.path.join(d, "pvc.log")
            m = run_commands([callers.pvc_call_from_alignments(
                ctx, reference, catalog, flat, regions, calls, log)], os.path.join(d, "logs"))
            return {"measurement": m, "rows": time_rows(row_base, m, "call-from-alignments")}

    elif arm == "pvc_fitted":
        def action(d):
            check_alignments(flat + [reference])
            warm(flat + [reference])  # every tool reads from memory (measure.warm)
            psp_dir = os.path.join(d, "psps")
            logs = os.path.join(d, "logs")
            m_gen = run_commands(callers.pvc_generate_psps(
                ctx, reference, catalog, files, regions, psp_dir, logs), logs,
                max_concurrent=ctx.threads)
            psps = [os.path.join(psp_dir, f"{s}.psp") for s in samples]
            params = os.path.join(d, "fitted.parameters.toml")
            m_fit = run_commands([callers.pvc_estimate_parameters(
                ctx, reference, catalog, psps, params, os.path.join(d, "fit.log"))], logs)
            m_call = run_commands([callers.pvc_call_from_psps(
                ctx, reference, catalog, psps, calls, os.path.join(d, "call.log"),
                parameters=params)], logs)
            total = merge(m_gen, m_fit, m_call)
            rows = (time_rows(row_base, m_gen, "generate-psps")
                    + time_rows(row_base, m_fit, "estimate-parameters")
                    + time_rows(row_base, m_call, "call-from-psps")
                    + time_rows(row_base, total, "whole pipeline"))
            return {"measurement": total, "rows": rows}

    elif arm == "gatk":
        def action(d):
            check_alignments(flat + [reference])
            warm(flat + [reference])  # every tool reads from memory (measure.warm)
            logs = os.path.join(d, "logs")
            os.makedirs(os.path.join(ctx.results_root, "java_tmp"), exist_ok=True)
            if not joint:
                (s,) = samples
                m = run_commands([callers.gatk_haplotypecaller(
                    ctx, reference, files[s][0], regions, calls, ctx.threads)], logs)
                return {"measurement": m, "rows": time_rows(row_base, m, "HaplotypeCaller")}
            gvcf_dir = os.path.join(d, "gvcfs")
            os.makedirs(gvcf_dir)
            m_hc = run_commands(callers.gatk_gvcfs(
                ctx, reference, {s: files[s] for s in samples}, regions, gvcf_dir), logs,
                max_concurrent=ctx.threads)
            combined = os.path.join(d, "combined.g.vcf.gz")
            comb, gen = callers.gatk_joint(
                ctx, reference, [os.path.join(gvcf_dir, f"{s}.g.vcf.gz") for s in samples],
                regions, combined, calls)
            m_imp = run_commands([comb], logs)
            m_gen = run_commands([gen], logs)
            total = merge(m_hc, m_imp, m_gen)
            rows = (time_rows(row_base, m_hc, "HaplotypeCaller GVCFs")
                    + time_rows(row_base, m_imp, "CombineGVCFs")
                    + time_rows(row_base, m_gen, "GenotypeGVCFs")
                    + time_rows(row_base, total, "whole pipeline"))
            return {"measurement": total, "rows": rows}

    elif arm == "deepvariant":
        def action(d):
            # DeepVariant runs in Google's own image through the host job runner (hostjobs.py),
            # one sample at a time with its own parallelism (`--num_shards`, one per core). A
            # cohort is called the way DeepVariant documents: a GVCF per sample, joined by GLnexus.
            hostjobs.require(ctx)
            check_alignments(flat + [reference])
            warm(flat + [reference])
            measurements, rows = [], []
            for s in samples:
                tmp = os.path.join(d, f"intermediate_{s}")
                out = calls if not joint else os.path.join(d, f"{s}.vcf.gz")
                gvcf = f" --output_gvcf={os.path.join(d, f'{s}.g.vcf.gz')}" if joint else ""
                cmd = (f"/opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref={reference} "
                       f"--reads={files[s][0]} --regions={regions} --output_vcf={out}{gvcf} "
                       f"--num_shards={ctx.threads} --intermediate_results_dir={tmp} "
                       f"--vcf_stats_report=false")
                m = hostjobs.run(ctx, "deepvariant", cmd, d, f"deepvariant.{s}")
                measurements.append(m)
                shutil.rmtree(tmp, ignore_errors=True)
            if joint:
                gvcfs = " ".join(os.path.join(d, f"{s}.g.vcf.gz") for s in samples)
                m = hostjobs.run(ctx, "glnexus",
                                 f"glnexus_cli --config DeepVariantWGS --bed {regions} "
                                 f"--threads {ctx.threads} {gvcfs} > {os.path.join(d, 'joint.bcf')}",
                                 d, "glnexus")
                measurements.append(m)
                run_commands([("bcf_to_vcf", f"bcftools view --no-version {os.path.join(d, 'joint.bcf')} "
                                             f"-Oz -o {calls}")], os.path.join(d, "logs"))
                shutil.rmtree(os.path.join(d, "GLnexus.DB"), ignore_errors=True)
            run_commands([("index", f"bcftools index -t -f {calls}")], os.path.join(d, "logs"))
            measured = [m for m in measurements if m]
            total = merge(*measured) if len(measured) == len(measurements) else None
            stage = "DeepVariant + GLnexus" if joint else "DeepVariant"
            return {"measurement": total, "rows": time_rows(row_base, total, stage) if total else []}

    elif arm == "freebayes":
        def action(d):
            check_alignments(flat + [reference])
            warm(flat + [reference])  # every tool reads from memory (measure.warm)
            logs = os.path.join(d, "logs")
            blocks = callers.split_regions(regions, os.path.join(d, "blocks"), ctx.threads * 4)
            out_dir = os.path.join(d, "block_vcfs")
            os.makedirs(out_dir)
            m_fb = run_commands(callers.freebayes_blocks(reference, flat, blocks, out_dir), logs,
                                max_concurrent=ctx.threads)
            vcfs = sorted(os.path.join(out_dir, f) for f in os.listdir(out_dir))
            m_join = run_commands([callers.freebayes_join(vcfs, calls)], logs)
            shutil.rmtree(out_dir)
            total = merge(m_fb, m_join)
            return {"measurement": total, "rows": time_rows(row_base, total, "freebayes")}
    else:
        raise ValueError(arm)

    step = Step(name, kind, action, fingerprint,
                attrs={"caller": ARM_CALLER[arm], "depth": base.get("depth")})
    return step, CallSet(arm, name, calls)


def score_step(ctx, callset, base, truth, bed, scoring, call_sample=None, truth_sample=None):
    """One call set scored against the truth by vcfeval and, where there is a checkout, by
    pvc_eval (lib/pvc_eval.py); lives beside the calls it scores."""
    kind = ARM_KIND[callset.arm]
    name = f"{callset.step_name}/score" + (f"_{call_sample}" if call_sample else "")
    row_base = dict(base, caller=ARM_LABELS[callset.arm],
                    caller_version=ctx.version(ARM_TOOL[callset.arm]))
    with_pvc_eval = scoring.str_regions is not None

    def fingerprint():
        return {"calls": file_identity(callset.path), "truth": file_identity(truth),
                "bed": file_identity(bed), "scorer": ctx.version("vcfeval"),
                "second_scorer": ctx.version("pvc_eval") if with_pvc_eval else None,
                "min_qual": MIN_QUAL.get(callset.arm), "samples": [call_sample, truth_sample]}

    def action(d):
        hostjobs.skip_with(ctx, callset)
        m, rows = vcfeval.score(ctx, d, row_base, truth, callset.path, bed, scoring.sdf,
                                min_qual=MIN_QUAL.get(callset.arm), call_sample=call_sample,
                                truth_sample=truth_sample)
        if with_pvc_eval:
            _, more = pvc_eval.score(ctx, d, row_base, truth, callset.path, bed, scoring.fasta,
                                     scoring.str_regions, min_qual=MIN_QUAL.get(callset.arm),
                                     call_sample=call_sample, truth_sample=truth_sample)
            rows += more
        return {"rows": rows}

    return Step(name, "score", action, fingerprint,
                attrs={"caller": ARM_CALLER[callset.arm], "depth": base.get("depth"),
                       "home": "pvc" if kind == "pvc" else "competitor"})
