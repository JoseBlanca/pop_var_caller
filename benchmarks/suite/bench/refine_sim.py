"""The simulated cohorts for genotype refinement and phase (the neighbour-evidence pass).

Three cohorts, drawn by `ng_str_stutter_benchmark` over the same 4 Mb genome as `stutter_sim`
(STR segments, STR bundles, SNPs and indels), each called by pvc with refinement off and on:

| cohort         | how it is drawn                                                     | samples, depth |
|----------------|---------------------------------------------------------------------|----------------|
| `f2`           | two inbred parents and their F2, one crossover a contig a meiosis   | parents at 30×, 20 F2 at 1× |
| `natural_F0`   | msprime coalescent, outbred (`scripts/msprime_haplotypes.py`)      | 400: 100 each at 30×, 20×, 10×, 5× |
| `natural_F0.9` | the same, selfing to F = 0.9, with autozygous runs                  | the same |

One contig plays a whole chromosome in the F2, so its crossovers are about 80 times closer together
than tomato's: a harder test than a real cross (owner, 2026-11-05).

| stage                                          | kind  | made                 |
|------------------------------------------------|-------|----------------------|
| haplotypes (msprime), simulate, map            | input | once                 |
| pvc: psps, parameters, calls off and on        | pvc   | every pvc commit     |
| scoring, each call set                         | score | every pvc commit     |

Refinement is forced on (`--refine-genotypes on`) rather than left to `auto`: an F2's inbreeding
against its own allele frequencies is about 0, below auto's floor, so auto would not run it there.

msprime runs inside the dev container through `uv`, pinned, its interpreter and package cached
under `inputs/uv`; the first run needs the network.
"""

import os
import shlex

from lib import arms, callers, pvc_eval
from lib.measure import check_alignments, run_commands, warm
from lib.results import row, time_rows
from lib.steps import Step, file_identity

TITLE = "Simulated cohorts for genotype refinement and phase"
DESCRIPTION = """
Simulated. Three cohorts over one drawn 4 Mb genome (STR segments and bundles, SNPs, indels), reads
mapped with minimap2: an F2 (two inbred parents at 30×, 20 F2 at 1×, one crossover per contig per
meiosis, so recombination far denser than in a real chromosome), and two natural populations from a
coalescent simulation (msprime), 400 samples each, 100 at each of 30×, 20×, 10× and 5×, one outbred
and one selfing to F = 0.9. pvc fits its parameters, then calls with refinement off and on. The
`depth` column is the group of samples. Genotypes are scored by pvc_eval, one unit a sample at a
locus, split by kind and by whether the sample had reads there; phase over SNPs only, a wrong link
being one switch between consecutive phased heterozygotes of one phase set.
"""

GENOME = "genome_bp=4000000 contigs=4"
SLICE_GENOME = "genome_bp=300000 contigs=1 segments_per_period=10 bundles=10"
MSPRIME = "uv run --no-project --python 3.12 --with msprime==1.3.3 python"
ARMS = {"off": "pvc (refinement off)", "on": "pvc (refinement on)"}
q = shlex.quote


def configurations(ctx):
    """(name, simulator settings, inbreeding for msprime or None, [(group, samples)])."""
    if ctx.slice:
        return [
            ("f2", "cross=f2 samples=6 deep_samples=2 deep_depth=30 depth=1", None,
             [("parents 30x", 2), ("F2 1x", 4)]),
            ("natural_F0", "samples=8 depths=30:2,20:2,10:2,5:2", 0.0,
             [("30x", 2), ("20x", 2), ("10x", 2), ("5x", 2)]),
        ]
    tiers = [("30x", 100), ("20x", 100), ("10x", 100), ("5x", 100)]
    depths = "depths=" + ",".join(f"{g[:-1]}:{n}" for g, n in tiers)
    return [
        ("f2", "cross=f2 samples=22 deep_samples=2 deep_depth=30 depth=1", None,
         [("parents 30x", 2), ("F2 1x", 20)]),
        ("natural_F0", f"samples=400 {depths}", 0.0, tiers),
        ("natural_F0.9", f"samples=400 {depths}", 0.9, tiers),
    ]


def steps(ctx):
    out = []
    for config in configurations(ctx):
        out.extend(_configuration(ctx, *config))
    return out


def _example(ctx, name):
    return os.path.join(ctx.project_dir, "target-container", "release", "examples", name)


def _configuration(ctx, cfg, settings, inbreeding, groups):
    name = f"refine_sim/{cfg}"
    sim_dir = os.path.join(ctx.results_root, "inputs", name, "simulate")
    bam_dir = os.path.join(ctx.results_root, "inputs", name, "map")
    reference = os.path.join(sim_dir, "reference.fa")
    confident = os.path.join(sim_dir, "confident.bed")
    truth = os.path.join(sim_dir, "truth.vcf.gz")
    n_samples = sum(n for _, n in groups)
    sample_names = [f"sim{i:03d}" for i in range(n_samples)]
    group_of = [g for g, n in groups for _ in range(n)]
    bams = {s: os.path.join(bam_dir, f"{s}.bam") for s in sample_names}
    genome = SLICE_GENOME if ctx.slice else GENOME
    sim_settings = f"{genome} {settings} threads={ctx.threads}"
    if inbreeding is not None:
        sim_settings += f" haplotypes={os.path.join(sim_dir, 'haplotypes.tsv')} inbreeding={inbreeding}"
    genome_bp = int(genome.split("genome_bp=")[1].split()[0])
    contigs = int(genome.split("contigs=")[1].split()[0])
    steps = []

    # ---- inputs: haplotypes, simulate, map -----------------------------------------------------
    def simulate(d):
        logs = os.path.join(d, "logs")
        cmds = []
        if inbreeding is not None:
            uv = os.path.join(ctx.results_root, "inputs", "uv")
            # The simulator's contigs come out near genome_bp / contigs long; half as much again
            # leaves room, and the simulator refuses a file whose contigs are too short.
            cmds.append(("haplotypes",
                         f"UV_CACHE_DIR={q(uv)}/cache UV_PYTHON_INSTALL_DIR={q(uv)}/python "
                         f"{MSPRIME} {q(os.path.join(ctx.project_dir, 'scripts', 'msprime_haplotypes.py'))} "
                         f"{q(os.path.join(d, 'haplotypes.tsv'))} --samples {n_samples} "
                         f"--contigs {contigs} --contig-bp {genome_bp * 3 // (2 * contigs)} "
                         f"--inbreeding {inbreeding}"))
        cmds.append(("simulate",
                     f"{q(_example(ctx, 'ng_str_stutter_benchmark'))} {q(d)} {sim_settings} "
                     f"> {q(os.path.join(d, 'simulate.log'))} "
                     f"&& samtools faidx {q(reference)} "
                     f"&& bgzip -f {q(os.path.join(d, 'truth.vcf'))} && tabix -p vcf {q(truth)}"))
        m = run_commands(cmds, logs)
        with open(os.path.join(d, "groups.tsv"), "w") as fh:
            for s, g in zip(sample_names, group_of):
                fh.write(f"{s}\t{g}\n")
        return {"measurement": m}

    steps.append(Step(f"{name}/simulate", "input", simulate,
                      {"settings": sim_settings.replace(f"threads={ctx.threads}", ""),
                       "msprime": MSPRIME if inbreeding is not None else None,
                       "simulator": file_identity(os.path.join(
                           ctx.project_dir, "examples", "ng_str_stutter_benchmark.rs")),
                       "haplotype_script": file_identity(os.path.join(
                           ctx.project_dir, "scripts", "msprime_haplotypes.py"))}))

    def map_reads(d):
        cmds = []
        for s in sample_names:
            r1, r2 = (os.path.join(sim_dir, f"{s}_R{i}.fastq.gz") for i in (1, 2))
            cmds.append((f"map.{s}",
                         f"minimap2 -ax sr -t {ctx.threads} -R '@RG\\tID:{s}\\tSM:{s}\\tLB:{s}\\tPL:ILLUMINA' "
                         f"{q(reference)} {q(r1)} {q(r2)} | samtools sort -@ 2 -o {q(bams[s])} - "
                         f"&& samtools index {q(bams[s])} && rm -f {q(r1)} {q(r2)}"))
        return {"measurement": run_commands(cmds, os.path.join(d, "logs"))}

    steps.append(Step(f"{name}/map", "input", map_reads,
                      lambda: {"reference": file_identity(reference + ".fai"),
                               "minimap2": ctx.version("minimap2")}))

    region_bp = 300000 if ctx.slice else genome_bp
    dataset = {"f2": "simulated F2", "natural_F0": "simulated natural population, outbred",
               "natural_F0.9": "simulated natural population, selfing F = 0.9"}[cfg]
    if ctx.slice:
        dataset += " (slice)"
    base = dict(benchmark="refine_sim", dataset=dataset, samples=n_samples, depth="all",
                region_bp=region_bp)

    # ---- pvc: the catalogue and pvc_eval's STR regions, then one fit and two calls --------------
    ref_steps, catalog, scoring = arms.reference_steps(ctx, f"refine_sim_{cfg}", reference)
    steps.extend(s for s in ref_steps if not s.name.endswith(".sdf"))

    pvc_dir = os.path.join(ctx.pvc_root(), name, "pvc")
    calls = {arm: os.path.join(pvc_dir, f"calls_refine_{arm}.vcf.gz") for arm in ARMS}

    def run_pvc(d):
        flat = list(bams.values())
        check_alignments(flat + [reference])
        warm(flat + [reference])
        logs = os.path.join(d, "logs")
        psp_dir = os.path.join(d, "psps")
        m_gen = run_commands(callers.pvc_generate_psps(
            ctx, reference, catalog, {s: [b] for s, b in bams.items()}, confident, psp_dir, logs),
            logs, max_concurrent=ctx.threads)
        psps = [os.path.join(psp_dir, f"{s}.psp") for s in sample_names]
        params = os.path.join(d, "fitted.parameters.toml")
        m_fit = run_commands([callers.pvc_estimate_parameters(
            ctx, reference, catalog, psps, params, os.path.join(d, "fit.log"))], logs)
        rows = (time_rows(dict(base, caller="pvc", caller_version=ctx.version("pvc")), m_gen,
                          "generate-psps")
                + time_rows(dict(base, caller="pvc", caller_version=ctx.version("pvc")), m_fit,
                            "estimate-parameters"))
        for arm, label in ARMS.items():
            _, command = callers.pvc_call_from_psps(
                ctx, reference, catalog, psps, calls[arm], os.path.join(d, f"call_{arm}.log"),
                parameters=params, extra=f"--refine-genotypes {arm}")
            m_call = run_commands([(f"call_from_psps_refine_{arm}", command)], logs)
            rows += time_rows(dict(base, caller=label, caller_version=ctx.version("pvc")), m_call,
                              "call-from-psps")
        return {"measurement": m_fit, "rows": rows}

    steps.append(Step(f"{name}/pvc", "pvc", run_pvc,
                      lambda: {"bams": [file_identity(b) for b in bams.values()],
                               "catalog": file_identity(catalog)},
                      attrs={"caller": "pvc"}))

    for arm, label in ARMS.items():
        steps.append(_score_step(ctx, name, arm, label, calls[arm], truth, confident, scoring,
                                 os.path.join(sim_dir, "groups.tsv"), base))
    return steps


def _score_step(ctx, name, arm, label, calls, truth, confident, scoring, groups, base):
    """pvc_eval over every sample at once, then lib/refine_score.py splits its table by group."""
    script = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib",
                          "refine_score.py")
    with_pvc_eval = scoring.str_regions is not None

    def action(d):
        logs = os.path.join(d, "logs")
        loci = ""
        if with_pvc_eval:
            kept = os.path.join(d, "calls.kept.vcf.gz")
            sorted_truth = os.path.join(d, "truth.sorted.vcf.gz")
            filtered = os.path.join(d, "calls.filtered.vcf.gz")
            loci = os.path.join(d, "loci.tsv")
            prepare = (f"bcftools view --no-version -f PASS,. -Oz -o {q(filtered)} {q(calls)} && "
                       f"{pvc_eval.in_reference_order(scoring.fasta, filtered, kept)} && "
                       f"rm {q(filtered)} && "
                       f"{pvc_eval.in_reference_order(scoring.fasta, truth, sorted_truth)}")
            compare = pvc_eval._command(
                ctx, f"compare --reference {q(scoring.fasta)} --a {q(sorted_truth)} --b {q(kept)} "
                     f"--regions {q(confident)} --str-regions {q(scoring.str_regions)} "
                     f"--all-rows --output {q(loci)} --summary {q(os.path.join(d, 'summary.tsv'))}")
            run_commands([("pvc_eval_prepare", prepare), ("pvc_eval_compare", compare)], logs)
        table = os.path.join(d, "scores.tsv")
        run_commands([("score", f"python3 {q(script)} {q(calls)} {q(truth)} {q(groups)} "
                                f"{q(table)} {q(loci)}")], logs)
        rows = []
        with open(table) as fh:
            next(fh)
            for line in fh:
                group, cls, stage, metric, value = line.rstrip("\n").split("\t")
                rows.append(row(base, depth=group, variant_class=cls, stage=stage, metric=metric,
                                value=int(value), caller=label,
                                caller_version=ctx.version("pvc")))
        if loci and os.path.exists(loci):
            os.remove(loci)  # millions of rows; scores.tsv keeps what is reported
        return {"rows": rows}

    return Step(f"{name}/score_refine_{arm}", "score", action,
                lambda: {"calls": file_identity(calls), "truth": file_identity(truth),
                         "scorer": file_identity(script),
                         "pvc_eval": ctx.version("pvc_eval") if with_pvc_eval else None},
                attrs={"caller": "pvc", "home": "pvc"})
