"""The simulated stuttering benchmark: STR genotypes with exact truth (plan §4.5).

A drawn 4 Mb genome with 3,000 STR segments and 500 STR bundles; reads drawn through PCR or
PCR-free stutter and mapped with minimap2 (doc/devel/ng/spec/str_stutter_benchmark.md). This runs
the stages of `scripts/ng_str_stutter_benchmark.sh`, split the suite's way:

| stage                                   | kind       | made                         |
|-----------------------------------------|------------|------------------------------|
| simulate, check the typing, map         | input      | once                         |
| GATK, freebayes, HipSTR (twice)         | competitor | once                         |
| repeat catalogue; pvc defaults, fitted, true stutter | pvc | every pvc commit          |
| scoring                                 | score      | competitors once, pvc per commit |

The scorer is the driver's, `scripts/ng_str_stutter_score.py`; each scoring pass gets a directory
of links holding the simulator's truth files and the callers' calls under the names it expects.
"""

import os
import shlex
import shutil

from lib import arms, callers
from lib.measure import check_alignments, run_commands, warm
from lib.results import row, time_rows
from lib.steps import Step, file_identity

TITLE = "Simulated stuttering benchmark — STR genotypes with exact truth"
DESCRIPTION = """
Simulated. A drawn 4 Mb genome with 3,000 STR segments (500 per motif length, 1 to 6 bases) and
500 STR bundles (two or three repeats a few bases apart), plus SNPs and indels; reads drawn through
PCR or PCR-free stutter, mapped with minimap2; 1 or 50 samples at 5× to 60×. The unit is one sample
at one STR segment (or bundle): `right` means both called sequences over the segment and 30 bases
each side are the truth's; `false_call` is a variant called where the truth has none. pvc (fitted)
fits its stutter from the run's own reads with `estimate-parameters`; pvc (true stutter) is given
the simulator's rates and is a reference row, not a user's result. freebayes is kept at QUAL ≥ 20.
HipSTR is given the segments only, not the bundles, once at its defaults and once with a fixed
stutter model and no read floor (its defaults call nothing below about 100 reads a locus).
"""

DESIGNS = ["pcr", "pcr_free"]
SAMPLE_COUNTS = [1, 50]
DEPTHS = [5, 15, 30, 60]
SLICE_SETTINGS = "segments_per_period=10 bundles=10 genome_bp=300000 contigs=1"
q = shlex.quote

PVC_ARMS = ["pvc_defaults", "pvc_fitted"]
COMPETITOR_ARMS = ["gatk", "freebayes"]
LABELS = {"ng_defaults": "pvc (defaults)", "ng_fitted": "pvc (fitted)",
          "ng_true_stutter": "pvc (true stutter)", "gatk": "GATK", "freebayes": "freebayes",
          "hipstr": "HipSTR", "hipstr_fixed_stutter": "HipSTR (fixed stutter)"}
TOOLS = {"ng_defaults": "pvc", "ng_fitted": "pvc", "ng_true_stutter": "pvc", "gatk": "gatk",
         "freebayes": "freebayes", "hipstr": "hipstr", "hipstr_fixed_stutter": "hipstr"}
OUTCOMES = ["right", "length_right", "wrong", "called_reference", "no_call", "false_call",
            "ref_right"]


def configurations(ctx):
    if ctx.slice:
        return [(d, 2, 5) for d in DESIGNS]
    return [(d, n, x) for d in DESIGNS for n in SAMPLE_COUNTS for x in DEPTHS]


def steps(ctx):
    out = []
    for design, samples, depth in configurations(ctx):
        out.extend(_configuration(ctx, design, samples, depth))
    return out


def _example(ctx, name):
    return os.path.join(ctx.project_dir, "target-container", "release", "examples", name)


def _configuration(ctx, design, n_samples, depth):
    cfg = f"{design}_{n_samples}s_{depth}x"
    name = f"stutter_sim/{cfg}"
    sim_dir = os.path.join(ctx.results_root, "inputs", name, "simulate")
    bam_dir = os.path.join(ctx.results_root, "inputs", name, "map")
    reference = os.path.join(sim_dir, "reference.fa")
    confident = os.path.join(sim_dir, "confident.bed")
    sample_names = [f"sim{i:03d}" for i in range(n_samples)]
    bams = {s: os.path.join(bam_dir, f"{s}.bam") for s in sample_names}
    settings = f"design={design} samples={n_samples} depth={depth} threads={ctx.threads}"
    if ctx.slice:
        settings += " " + SLICE_SETTINGS
    steps = []

    # ---- inputs: simulate (with the typing checked against what was planted), then map --------
    def simulate(d):
        logs = os.path.join(d, "logs")
        m = run_commands([("simulate", f"{q(_example(ctx, 'ng_str_stutter_benchmark'))} {q(d)} "
                                       f"{settings} > {q(os.path.join(d, 'simulate.log'))} "
                                       f"&& samtools faidx {q(reference)} "
                                       f"&& samtools dict {q(reference)} -o {q(d)}/reference.dict")],
                         logs)
        # The typed STR segments must be the planted ones, and each planted bundle typed as one
        # bundle with the caller's reference runs at the simulator's stutter sites; the driver's
        # check, unchanged. Needs a catalogue, built here only for this check.
        catalog = os.path.join(d, "check_catalog.parquet")
        check = (
            f"{q(ctx.pvc_bin)} repeat-catalog --reference {q(reference)} --output {q(catalog)} --force "
            f"&& ln -sf {q(catalog)} {q(reference)}.repeats.parquet "
            f"&& NG_REFERENCE_CHECK=skip {q(_example(ctx, 'ng_typed_region_dump'))} {q(reference)} "
            f"{q(confident)} calling > {q(d)}/typing_check.tsv "
            f"&& diff <(awk -F'\\t' '$4==\"ssr_locus\" {{print $1\"\\t\"($2-1)\"\\t\"$3\"\\t\"$6}}' "
            f"{q(d)}/typing_check.tsv | sort) <(sort {q(d)}/segments.bed) "
            f"&& NG_REFERENCE_CHECK=skip {q(_example(ctx, 'ng_str_bundle_runs'))} {q(reference)} "
            f"{q(confident)} > {q(d)}/bundle_runs_check.tsv "
            f"&& python3 {q(os.path.join(ctx.project_dir, 'scripts', 'ng_str_stutter_check_bundles.py'))} {q(d)} "
            f"&& rm -f {q(reference)}.repeats.parquet {q(catalog)}")
        run_commands([("check_typing", check)], logs)
        return {"measurement": m}

    steps.append(Step(f"{name}/simulate", "input", simulate,
                      {"settings": settings.replace(f"threads={ctx.threads}", ""),
                       "simulator": file_identity(os.path.join(
                           ctx.project_dir, "examples", "ng_str_stutter_benchmark.rs"))}))

    def map_reads(d):
        cmds = []
        for s in sample_names:
            r1, r2 = (os.path.join(sim_dir, f"{s}_R{i}.fastq.gz") for i in (1, 2))
            cmds.append((f"map.{s}",
                         f"minimap2 -ax sr -t {ctx.threads} -R '@RG\\tID:{s}\\tSM:{s}\\tLB:{s}\\tPL:ILLUMINA' "
                         f"{q(reference)} {q(r1)} {q(r2)} | samtools sort -@ 2 -o {q(bams[s])} - "
                         f"&& samtools index {q(bams[s])} && rm -f {q(r1)} {q(r2)}"))
        m = run_commands(cmds, os.path.join(d, "logs"))
        return {"measurement": m}

    steps.append(Step(f"{name}/map", "input", map_reads,
                      lambda: {"reference": file_identity(reference + ".fai"),
                               "minimap2": ctx.version("minimap2")}))

    base = dict(benchmark="stutter_sim", dataset=f"simulated 4 Mb, {design.replace('_', '-')} "
                                                 f"stutter", samples=n_samples, depth=f"{depth}x",
                region_bp=300000 if ctx.slice else 4000000)
    if ctx.slice:
        base["dataset"] = f"simulated 300 kb (slice), {design.replace('_', '-')} stutter"

    # ---- pvc, every commit --------------------------------------------------------------------
    catalog = os.path.join(ctx.pvc_root(), name, "catalog", "reference.repeats.parquet")

    def make_catalog(d):
        m = run_commands([callers.pvc_catalog(ctx, reference, catalog, os.path.join(d, "catalog.log"))],
                         os.path.join(d, "logs"))
        return {"measurement": m}

    steps.append(Step(f"{name}/catalog", "pvc", make_catalog,
                      lambda: {"reference": file_identity(reference + ".fai")},
                      attrs={"caller": "pvc", "depth": f"{depth}x"}))

    callsets = {}
    for arm in PVC_ARMS + COMPETITOR_ARMS:
        step, callset = arms.arm_steps(ctx, arm, name, reference, catalog, bams, confident, base,
                                       joint=True)
        steps.append(step)
        callsets[arm] = callset

    true_step, true_calls = _true_stutter_step(ctx, name, sim_dir, reference, catalog, callsets,
                                               bams, confident, base)
    steps.append(true_step)
    hipstr_steps = [_hipstr_step(ctx, name, sim_dir, reference, bams, base, fixed)
                    for fixed in (False, True)]
    steps.extend(s for s, _ in hipstr_steps)

    pvc_calls = {"ng_defaults": callsets["pvc_defaults"].path,
                 "ng_fitted": callsets["pvc_fitted"].path, "ng_true_stutter": true_calls}
    competitor_calls = {"gatk": callsets["gatk"].path, "freebayes": callsets["freebayes"].path,
                        "hipstr": hipstr_steps[0][1], "hipstr_fixed_stutter": hipstr_steps[1][1]}
    # Once pvc's calls are scored its psps are deleted: about 100 MB a sample at 30x, some 25 GB
    # over the grid for every pvc commit. The true-stutter step remakes them if it must rerun.
    steps.append(_score_step(ctx, name, "pvc", sim_dir, pvc_calls, base,
                             prune=os.path.join(ctx.pvc_root(), name, "pvc_fitted", "psps")))
    steps.append(_score_step(ctx, name, "competitor", sim_dir, competitor_calls, base))
    return steps


def _true_stutter_step(ctx, name, sim_dir, reference, catalog, callsets, bams, confident, base):
    """pvc given the simulator's own stutter: the defaults' parameters with the stutter tables
    replaced (scripts/ng_str_stutter_true_parameters.py), over the fitted arm's psps — or over
    psps it makes again, untimed, when scoring has already deleted those."""
    samples = list(bams)
    step_name = f"{name}/pvc_true_stutter"
    calls = os.path.join(ctx.pvc_root(), step_name, "calls.vcf.gz")
    defaults_log = os.path.join(ctx.pvc_root(), name, "pvc_defaults", "pvc.log")
    psp_dir = os.path.join(ctx.pvc_root(), name, "pvc_fitted", "psps")
    row_base = dict(base, caller="pvc (true stutter)", caller_version=ctx.version("pvc"))

    def action(d):
        view = os.path.join(d, "view")
        os.makedirs(view)
        for f in ("stutter_rates.tsv", "libraries.tsv"):
            os.symlink(os.path.join(sim_dir, f), os.path.join(view, f))
        script = os.path.join(ctx.project_dir, "scripts", "ng_str_stutter_true_parameters.py")
        params = os.path.join(view, "true_stutter.parameters.toml")
        prep = (f"sed -n '/^# ---- parameters file begins ----$/,/^# ---- parameters file ends ----$/p' "
                f"{q(defaults_log)} | sed '1d;$d' > {q(view)}/ng_defaults.parameters.toml "
                f"&& python3 {q(script)} {q(view)}")
        run_commands([("true_parameters", prep)], os.path.join(d, "logs"))
        source = psp_dir
        if not all(os.path.exists(os.path.join(psp_dir, f"{s}.psp")) for s in samples):
            source = os.path.join(d, "psps")
            run_commands(callers.pvc_generate_psps(ctx, reference, catalog,
                                                   {k: [v] for k, v in bams.items()}, confident, source,
                                                   os.path.join(d, "logs")),
                         os.path.join(d, "logs"), max_concurrent=ctx.threads)
        psps = [os.path.join(source, f"{s}.psp") for s in samples]
        m = run_commands([callers.pvc_call_from_psps(ctx, reference, catalog, psps, calls,
                                                     os.path.join(d, "call.log"), parameters=params)],
                         os.path.join(d, "logs"))
        if source != psp_dir:
            shutil.rmtree(source)
        return {"measurement": m, "rows": time_rows(row_base, m, "call-from-psps")}

    step = Step(step_name, "pvc", action,
                lambda: {"defaults": file_identity(callsets["pvc_defaults"].path),
                         "fitted": file_identity(callsets["pvc_fitted"].path)},
                attrs={"caller": "pvc", "depth": base["depth"]})
    return step, calls


def _hipstr_step(ctx, name, sim_dir, reference, bams, base, fixed):
    """HipSTR over the STR segments, all samples jointly; given no bundles (one period a region)."""
    arm = "hipstr_fixed_stutter" if fixed else "hipstr"
    step_name = f"{name}/{arm}"
    d_final = os.path.join(ctx.results_root, "competitors", step_name)
    calls = os.path.join(d_final, "calls.vcf.gz")
    extra = "--def-stutter-model --min-reads 1" if fixed else ""
    row_base = dict(base, caller=LABELS[arm], caller_version=ctx.version("hipstr"))

    def action(d):
        check_alignments(list(bams.values()))
        warm(list(bams.values()) + [reference])
        regions = os.path.join(d, "hipstr_regions.bed")
        segments = os.path.join(sim_dir, "segments.bed")
        names = ",".join(bams)
        cmd = (f"awk -F'\\t' -v OFS='\\t' '{{ print $1, $2 + 1, $3, $4, ($3 - $2) / $4, "
               f"sprintf(\"seg%05d\", NR - 1) }}' {q(segments)} > {q(regions)} "
               f"&& HipSTR --bams {q(','.join(bams.values()))} --bam-samps {names} "
               f"--bam-libs {names} --fasta {q(reference)} --regions {q(regions)} "
               f"--str-vcf {q(calls)} --log {q(os.path.join(d, 'hipstr.log'))} {extra} "
               f"&& tabix -p vcf {q(calls)}")
        m = run_commands([("hipstr", cmd)], os.path.join(d, "logs"))
        return {"measurement": m, "rows": time_rows(row_base, m, "HipSTR")}

    step = Step(step_name, "competitor", action,
                lambda: {"tool": ctx.version("hipstr"), "args": extra,
                         "bams": [file_identity(b) for b in bams.values()],
                         "segments": file_identity(os.path.join(sim_dir, "segments.bed"))},
                attrs={"caller": "hipstr", "depth": base["depth"]})
    return step, calls


def _score_step(ctx, name, home, sim_dir, calls, base, prune=None):
    """The driver's scorer over one side's callers, in a directory of links."""
    script = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib",
                          "stutter_score_callers.py")

    def action(d):
        view = os.path.join(d, "view")
        os.makedirs(view)
        for f in os.listdir(sim_dir):
            if not f.endswith(".fastq.gz") and f not in ("logs", ".done"):
                os.symlink(os.path.join(sim_dir, f), os.path.join(view, f))
        args = []
        for caller, path in calls.items():
            if not os.path.exists(path):
                continue
            target = os.path.join(view, f"{caller}.vcf.gz")
            if caller == "freebayes":
                # freebayes writes every candidate it considers; kept at QUAL >= 20 (plan D2).
                run_commands([("freebayes_qual", f"bcftools view --no-version -i 'QUAL>=20' "
                                                 f"{q(path)} -Oz -o {q(target)}")],
                             os.path.join(d, "logs"))
            else:
                os.symlink(path, target)
            args.append(f"{caller}={caller}.vcf.gz")
        run_commands([("score", f"python3 {q(script)} {q(view)} {' '.join(args)}")],
                     os.path.join(d, "logs"))
        rows = _rows(view, base, ctx)
        if prune and os.path.isdir(prune):
            shutil.rmtree(prune)
        return {"rows": rows}

    callers_for_attr = "pvc" if home == "pvc" else None
    attrs = {"depth": base["depth"], "home": home}
    if callers_for_attr:
        attrs["caller"] = callers_for_attr
    return Step(f"{name}/score_{home}", "score", action,
                lambda: {c: file_identity(p) for c, p in calls.items() if os.path.exists(p)}
                | {"version": 1},
                attrs=attrs)


def _rows(view, base, ctx):
    """Counts by caller and outcome, for segments (all, and by motif length) and for bundles."""
    out = []
    for cells_file, unit in (("score_cells.tsv", "STR segment"), ("score_bundle_cells.tsv", "STR bundle")):
        path = os.path.join(view, cells_file)
        if not os.path.exists(path):
            continue
        counts = {}
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                rec = dict(zip(header, line.rstrip("\n").split("\t")))
                period = rec.get("period") or rec.get("longest_member_period")
                for cls in (unit, f"{unit}, motif length {period}" if unit == "STR segment" else None):
                    if cls is None:
                        continue
                    key = (rec["caller"], cls)
                    c = counts.setdefault(key, {o: 0 for o in OUTCOMES} | {"variant": 0, "reference": 0})
                    c[rec["outcome"]] = c.get(rec["outcome"], 0) + 1
                    c[rec["truth"]] += 1
        for (caller, cls), c in counts.items():
            row_base = dict(base, caller=LABELS.get(caller, caller),
                            caller_version=ctx.version(TOOLS.get(caller, "pvc")),
                            variant_class=cls, stage="genotype")
            out.append(row(row_base, metric="variant_pairs", value=c["variant"]))
            out.append(row(row_base, metric="reference_pairs", value=c["reference"]))
            for o in OUTCOMES:
                out.append(row(row_base, metric=o, value=c.get(o, 0)))
    return out
