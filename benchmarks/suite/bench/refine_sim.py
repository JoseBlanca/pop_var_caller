"""The simulated cohorts for genotype refinement and phase (the neighbour-evidence pass).

Three cohorts, drawn by `ng_str_stutter_benchmark` over the same 4 Mb genome as `stutter_sim`
(STR segments, STR bundles, SNPs and indels), each called by pvc with refinement off and on, and by
the tools that do the same job without a reference panel:

| cohort         | how it is drawn                                                     | samples, depth |
|----------------|---------------------------------------------------------------------|----------------|
| `f2`           | two inbred parents and their F2, one crossover a contig a meiosis   | parents at 30×, 100 F2 at 1× |
| `natural_F0`   | msprime coalescent, outbred (`scripts/msprime_haplotypes.py`)      | 400: 100 each at 30×, 10×, 5×, 2× |
| `natural_F0.9` | the same, selfing to F = 0.9, with autozygous runs                  | the same |

One contig plays a whole chromosome in the F2, so its crossovers are about 80 times closer together
than tomato's: a harder test than a real cross (owner, 2026-11-05).

| stage                                          | kind       | made                 |
|------------------------------------------------|------------|----------------------|
| haplotypes (msprime), simulate, map, coverage  | input      | once                 |
| pvc: psps, parameters, calls off and on        | pvc        | every pvc commit     |
| GATK; GATK then Beagle 4.1; STITCH             | competitor | once                 |
| scoring, each call set                         | score      | with its call set    |

**The competitors.** GATK joint calling as the other benchmarks run it (lib/arms.py), then Beagle
4.1 imputing and phasing from GATK's genotype likelihoods (`gl=`: Beagle 5 takes hard genotypes
only); and STITCH, which imputes from the reads alone at GATK's biallelic SNP sites, in its Bioconda
container through the host job runner (lib/hostjobs.py). Both are given the simulation's own
recombination rate, and for the natural populations its effective size: their best case, not a
user's guess. STITCH calls SNPs only, so only its SNP rows are kept.

**With or without reads** is read from the alignments (the coverage step: any mapped read over any
base of the locus), the same for every caller, since Beagle and STITCH write no DP.

Refinement is forced on (`--refine-genotypes on`) rather than left to `auto`: an F2's inbreeding
against its own allele frequencies is about 0, below auto's floor, so auto would not run it there.

msprime runs inside the dev container through `uv`, pinned, its interpreter and package cached
under `inputs/uv`; the first run needs the network.
"""

import os
import shlex
import time

from lib import arms, callers, hostjobs, pvc_eval
from lib.measure import check_alignments, merge, run_commands, warm
from lib.results import row, time_rows
from lib.steps import Step, file_identity

TITLE = "Simulated cohorts for genotype refinement and phase"
DESCRIPTION = """
Simulated. Three cohorts over one drawn 4 Mb genome (STR segments and bundles, SNPs, indels), reads
mapped with minimap2: an F2 (two inbred parents at 30×, 100 F2 at 1×, one crossover per contig per
meiosis, so recombination far denser than in a real chromosome), and two natural populations from a
coalescent simulation (msprime), 400 samples each, 100 at each of 30×, 10×, 5× and 2×, one outbred
and one selfing to F = 0.9. pvc fits its parameters, then calls with refinement off and on. The
`depth` column is the group of samples. Competitors: GATK; GATK then Beagle 4.1 from GATK's
genotype likelihoods; STITCH from the reads at GATK's SNP sites (SNPs only). Beagle and STITCH are
given the simulation's true recombination rate and effective size. Genotypes are scored by
pvc_eval, one unit a sample at a locus, split by kind and by whether the sample had any mapped read
over the locus; phase over SNPs only, a wrong link being one switch between consecutive phased
heterozygotes of one phase set.
"""

GENOME = "genome_bp=4000000 contigs=4"
SLICE_GENOME = "genome_bp=300000 contigs=1 segments_per_period=10 bundles=10"
MSPRIME = "uv run --no-project --python 3.12 --with msprime==1.3.3 python"
ARMS = {"off": "pvc (refinement off)", "on": "pvc (refinement on)"}
# The simulated reads are made once and kept. Their fingerprint is the settings and this number, not
# the simulator's or the haplotype script's source: a change to either does not remake them. Raise
# it when such a change should reach the benchmark. The sources' checksums are kept in the step's
# record (`meta`), for the record only.
SIMULATION_VERSION = 1
# msprime_haplotypes.py's defaults, which this module leaves: the outbred effective size, and
# recombination per base a generation (1 cM/Mb). Selfing to F scales them to Ne / (1 + F) and
# r (1 − F), as the script does.
MSPRIME_NE = 10_000
MSPRIME_RECOMBINATION = 1e-8
# STITCH's founder haplotypes: the two inbred parents in the F2; in a natural population a mosaic of
# K, with nGen = 4 Ne / K, the rule STITCH's own help gives.
STITCH_K_NATURAL = 20
COMPETITORS = {"gatk": "GATK", "beagle": "GATK + Beagle 4.1", "stitch": "STITCH"}
# What counts as reads over a locus: `samtools depth`'s default reads (mapped, primary, not
# duplicate, not QC-failed), at any base of it, a base a read deletes included (`-J`), so a sample
# carrying a deletion has reads over it. Raise when the rule changes.
COVERAGE_RULE = 1
q = shlex.quote


def configurations(ctx):
    """(name, simulator settings, inbreeding for msprime or None, [(group, samples)])."""
    if ctx.slice:
        return [
            ("f2", "cross=f2 samples=6 deep_samples=2 deep_depth=30 depth=1", None,
             [("parents 30x", 2), ("F2 1x", 4)]),
            ("natural_F0", "samples=8 depths=30:2,10:2,5:2,2:2", 0.0,
             [("30x", 2), ("10x", 2), ("5x", 2), ("2x", 2)]),
        ]
    # 30× checks that refinement does no harm to well-covered samples; 2× is where it has the most
    # to do (owner, 2026-11-05: 2× replaced 20×, which added least).
    tiers = [("30x", 100), ("10x", 100), ("5x", 100), ("2x", 100)]
    depths = "depths=" + ",".join(f"{g[:-1]}:{n}" for g, n in tiers)
    return [
        ("f2", "cross=f2 samples=102 deep_samples=2 deep_depth=30 depth=1", None,
         [("parents 30x", 2), ("F2 1x", 100)]),
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
                     f"&& samtools dict {q(reference)} -o {q(os.path.join(d, 'reference.dict'))} "
                     f"&& bgzip -f {q(os.path.join(d, 'truth.vcf'))} && tabix -p vcf {q(truth)}"))
        m = run_commands(cmds, logs)
        with open(os.path.join(d, "groups.tsv"), "w") as fh:
            for s, g in zip(sample_names, group_of):
                fh.write(f"{s}\t{g}\n")
        return {"measurement": m,
                "meta": {"simulator": file_identity(os.path.join(
                             ctx.project_dir, "examples", "ng_str_stutter_benchmark.rs")),
                         "haplotype_script": file_identity(os.path.join(
                             ctx.project_dir, "scripts", "msprime_haplotypes.py"))}}

    steps.append(Step(f"{name}/simulate", "input", simulate,
                      {"settings": sim_settings.replace(f"threads={ctx.threads}", ""),
                       "msprime": MSPRIME if inbreeding is not None else None,
                       "simulation": SIMULATION_VERSION}))

    def map_reads(d):
        cmds = []
        for s in sample_names:
            r1, r2 = (os.path.join(sim_dir, f"{s}_R{i}.fastq.gz") for i in (1, 2))
            cmds.append((f"map.{s}",
                         f"minimap2 -ax sr -t {ctx.threads} -R '@RG\\tID:{s}\\tSM:{s}\\tLB:{s}\\tPL:ILLUMINA' "
                         f"{q(reference)} {q(r1)} {q(r2)} | samtools sort -@ 2 -o {q(bams[s])} - "
                         f"&& samtools index {q(bams[s])} && rm -f {q(r1)} {q(r2)}"))
        return {"measurement": run_commands(cmds, os.path.join(d, "logs"))}

    # The settings too: reads drawn again over the same reference must be mapped again.
    steps.append(Step(f"{name}/map", "input", map_reads,
                      lambda: {"reference": file_identity(reference + ".fai"),
                               "minimap2": ctx.version("minimap2"),
                               "settings": sim_settings.replace(f"threads={ctx.threads}", ""),
                               "simulation": SIMULATION_VERSION}))

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

    # ---- which sample had reads where, from the alignments -------------------------------------
    cov_dir = os.path.join(ctx.results_root, "inputs", name, "coverage")

    def coverage(d):
        runs = ("awk -v OFS='\\t' '$1 == c && $2 == e + 1 { e = $2; next } "
                "{ if (c != \"\") print c, s - 1, e; c = $1; s = $2; e = $2 } "
                "END { if (c != \"\") print c, s - 1, e }'")
        cmds = [(f"coverage.{s}",
                 f"samtools depth -J {q(b)} | {runs} > {q(os.path.join(d, s + '.bed'))}")
                for s, b in bams.items()]
        return {"measurement": run_commands(cmds, os.path.join(d, "logs"),
                                            max_concurrent=ctx.threads)}

    steps.append(Step(f"{name}/coverage", "input", coverage,
                      lambda: {"bams": [file_identity(b) for b in bams.values()],
                               "rule": COVERAGE_RULE}))
    coverage_done = os.path.join(cov_dir, ".done")
    groups_tsv = os.path.join(sim_dir, "groups.tsv")

    for arm, label in ARMS.items():
        steps.append(_score_step(ctx, f"{name}/score_refine_{arm}", label, ctx.version("pvc"),
                                 calls[arm], truth, confident, scoring, groups_tsv, cov_dir,
                                 coverage_done, base, home="pvc", caller="pvc"))

    # ---- competitors, run once ------------------------------------------------------------------
    gatk_step, gatk_set = arms.arm_steps(ctx, "gatk", name, reference, catalog,
                                         {s: [b] for s, b in bams.items()}, confident, base,
                                         joint=True)
    steps.append(gatk_step)
    comp_dir = os.path.join(ctx.results_root, "competitors", name)
    contigs_fai = reference + ".fai"
    ne = MSPRIME_NE / (1 + inbreeding) if inbreeding is not None else None
    # cM per base a meiosis: the F2's one crossover per contig, so 100 cM over each contig's length;
    # the populations' effective rate.
    per_bp_cm = (None if inbreeding is None
                 else MSPRIME_RECOMBINATION * (1 - inbreeding) * 100)
    beagle_calls = os.path.join(comp_dir, "beagle", "calls.vcf.gz")
    stitch_calls = os.path.join(comp_dir, "stitch", "calls.vcf.gz")
    gatk_version = ctx.version("gatk")

    def contig_lengths():
        with open(contigs_fai) as fh:
            return [(f.split("\t")[0], int(f.split("\t")[1])) for f in fh]

    def cm_per_mb(length):
        return 100 / (length / 1e6) if per_bp_cm is None else per_bp_cm * 1e6

    def beagle(d):
        logs = os.path.join(d, "logs")
        genetic_map = os.path.join(d, "genetic.map")
        with open(genetic_map, "w") as fh:  # PLINK map: contig, id, cM, bp; Beagle interpolates
            for chrom, length in contig_lengths():
                fh.write(f"{chrom}\t.\t0\t1\n{chrom}\t.\t{cm_per_mb(length) * length / 1e6:.6f}"
                         f"\t{length}\n")
        gl = os.path.join(d, "gl.vcf.gz")
        out = os.path.join(d, "imputed")
        phased = os.path.join(d, "phased")
        ne_arg = f" ne={int(round(ne))}" if ne is not None else ""
        beagle_cmd = (f"java -Xmx{max(8, ctx.threads)}g -jar {q(ctx.beagle_jar)} "
                      f"map={q(genetic_map)}{ne_arg} nthreads={ctx.threads}")
        m_prep = run_commands([("gl_input", f"bcftools view --no-version -e 'ALT=\"*\"' "
                                            f"-Oz -o {q(gl)} {q(gatk_set.path)}")], logs)
        # Two passes, as Beagle 4.1 is used: genotypes from the likelihoods (`gl=`), which it
        # writes unphased, then phase over those genotypes (`gt=`).
        m = run_commands([("beagle", f"{beagle_cmd} gl={q(gl)} out={q(out)}")], logs)
        m_phase = run_commands([("beagle_phase",
                                 f"{beagle_cmd} gt={q(out + '.vcf.gz')} out={q(phased)}")], logs)
        run_commands([("index", f"mv {q(phased + '.vcf.gz')} {q(beagle_calls)} && "
                                f"bcftools index -t -f {q(beagle_calls)}")], logs)
        row_base = dict(base, caller=COMPETITORS["beagle"],
                        caller_version=f"{gatk_version} + {ctx.version('beagle')}")
        return {"measurement": merge(m, m_phase),
                "rows": time_rows(row_base, m_prep, "Beagle input")
                + time_rows(row_base, m, "Beagle genotypes")
                + time_rows(row_base, m_phase, "Beagle phase")
                + time_rows(row_base, merge(m, m_phase), "Beagle")}

    steps.append(Step(f"{name}/beagle", "competitor", beagle,
                      lambda: {"gatk_calls": file_identity(gatk_set.path),
                               "beagle": file_identity(ctx.beagle_jar), "ne": ne,
                               "cm_per_bp": per_bp_cm, "threads": ctx.threads},
                      attrs={"caller": "beagle"}))

    def stitch(d):
        logs = os.path.join(d, "logs")
        bamlist = os.path.join(d, "bams.txt")
        with open(bamlist, "w") as fh:
            fh.writelines(f"{bams[s]}\n" for s in sample_names)
        k = 2 if inbreeding is None else STITCH_K_NATURAL
        n_gen = 2 if inbreeding is None else max(1, round(4 * ne / k))
        vcfs, started = [], time.time()
        for chrom, length in contig_lengths():
            pos = os.path.join(d, f"{chrom}.pos")
            run_commands([(f"posfile.{chrom}",
                           f"bcftools view --no-version -r {q(chrom)} -m2 -M2 -v snps "
                           f"{q(gatk_set.path)} "
                           f"| bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' "
                           f"> {q(pos)}")], logs)
            if os.path.getsize(pos) == 0:
                continue
            out = os.path.join(d, chrom)
            os.makedirs(os.path.join(out, "tmp"))
            hostjobs.run(ctx, "stitch",
                         f"STITCH.R --chr={chrom} --bamlist={bamlist} --posfile={pos} "
                         f"--outputdir={out} --tempdir={os.path.join(out, 'tmp')} --K={k} "
                         f"--nGen={n_gen} --expRate={cm_per_mb(length):.6f} --method=diploid "
                         f"--nCores={ctx.threads}", d, f"stitch.{chrom}")
            vcfs.append(os.path.join(out, f"stitch.{chrom}.vcf.gz"))
        wall = time.time() - started
        # STITCH writes no contig lines, which bcftools needs to join the contigs' files.
        named = [v.replace(".vcf.gz", ".contigs.vcf.gz") for v in vcfs]
        run_commands([(f"contigs.{i}", f"bcftools reheader --fai {q(contigs_fai)} -o {q(n)} {q(v)} "
                                       f"&& bcftools index -t -f {q(n)}")
                      for i, (v, n) in enumerate(zip(vcfs, named))], logs)
        vcfs = named
        run_commands([("concat", f"bcftools concat --no-version -Oz -o {q(stitch_calls)} "
                                 + " ".join(q(v) for v in vcfs)
                                 + f" && bcftools index -t -f {q(stitch_calls)}")], logs)
        row_base = dict(base, caller=COMPETITORS["stitch"], caller_version=ctx.version("stitch"))
        # STITCH's image has no Python, so the suite cannot measure it inside; its wall time only.
        return {"rows": [row(row_base, stage="STITCH", metric="wall_seconds",
                             value=round(wall, 3))],
                "meta": {"K": k, "nGen": n_gen}}

    steps.append(Step(f"{name}/stitch", "competitor", stitch,
                      lambda: {"gatk_calls": file_identity(gatk_set.path),
                               "bams": [file_identity(b) for b in bams.values()],
                               "stitch": ctx.version("stitch"), "ne": ne, "cm_per_bp": per_bp_cm,
                               "K_natural": STITCH_K_NATURAL, "threads": ctx.threads},
                      attrs={"caller": "stitch"}))

    for key, path in (("gatk", gatk_set.path), ("beagle", beagle_calls), ("stitch", stitch_calls)):
        version = {"gatk": gatk_version, "beagle": f"{gatk_version} + {ctx.version('beagle')}",
                   "stitch": ctx.version("stitch")}[key]
        steps.append(_score_step(ctx, f"{name}/score_{key}", COMPETITORS[key], version, path,
                                 truth, confident, scoring, groups_tsv, cov_dir, coverage_done,
                                 base, home="competitor", caller=key, snps_only=key == "stitch"))
    return steps


def _score_step(ctx, step_name, label, version, calls, truth, confident, scoring, groups, cov_dir,
                coverage_done, base, home, caller, snps_only=False):
    """pvc_eval over every sample at once, then lib/refine_score.py splits its table by group.

    `snps_only` keeps the SNP rows alone, for a caller that calls nothing else (STITCH): its
    indels and repeats would all count as missed."""
    script = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib",
                          "refine_score.py")
    with_pvc_eval = scoring.str_regions is not None

    def action(d):
        if caller == "stitch":
            hostjobs.require(ctx)
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
                                f"{q(cov_dir)} {q(table)} {q(loci)}")], logs)
        rows = []
        with open(table) as fh:
            next(fh)
            for line in fh:
                group, cls, stage, metric, value = line.rstrip("\n").split("\t")
                if snps_only and not cls.startswith("SNP"):
                    continue
                rows.append(row(base, depth=group, variant_class=cls, stage=stage, metric=metric,
                                value=int(value), caller=label, caller_version=version))
        if loci and os.path.exists(loci):
            os.remove(loci)  # millions of rows; scores.tsv keeps what is reported
        return {"rows": rows}

    return Step(step_name, "score", action,
                lambda: {"calls": file_identity(calls), "truth": file_identity(truth),
                         "scorer": file_identity(script), "coverage": file_identity(coverage_done),
                         "snps_only": snps_only,
                         "pvc_eval": ctx.version("pvc_eval") if with_pvc_eval else None},
                attrs={"caller": caller, "home": home})
