"""Simulated tomato population with known inbreeding: does a wrong F cost genotypes?

50 diploid samples, 10 at each of F = 0, 0.25, 0.5, 0.75 and 1, simulated by
genome-population-simulator (lib/popsim.py) over 5 Mb of tomato chromosome 1, with the same SNPs at
the same frequencies in every group. Reads are simulated once at 30x, mapped with minimap2 and
downsampled to each depth, so every depth has the same true genotypes.

| stage                                     | kind       | made                              |
|-------------------------------------------|------------|-----------------------------------|
| simulate, map, downsample                 | input      | once                              |
| GATK (GVCFs, joint), DeepVariant + GLnexus | competitor | once                              |
| catalogue; pvc defaults, fitted, F stated | pvc        | every pvc commit                  |
| genotype scoring; pvc's F against truth   | score      | beside what they score            |

pvc runs three ways: its defaults, which take F = 0 for every sample; fitted, which estimates each
sample's F; and with F stated as 0.9 for every sample (`estimate-parameters --inbreeding 0.9`),
the opposite mistake for the outbred groups. GATK and DeepVariant take no F at all.
"""

import os
import shlex
import subprocess

from lib import arms, callers, genotype_concordance, hostjobs, inputs, popsim
from lib.measure import run_commands
from lib.results import row, time_rows
from lib.steps import Step, file_identity

from bench import speed

TITLE = "Simulated tomato population — genotypes under known inbreeding"
DESCRIPTION = """
Simulated (genome-population-simulator). 50 diploid samples over 5 Mb of tomato chromosome 1
(SL4.0ch01:5,000,001-10,000,000, no gaps), 10 samples at each inbreeding coefficient F = 0, 0.25,
0.5, 0.75 and 1, all called together. One SNP a kilobase, the same SNPs at the same frequencies in
every group (ALT frequency uniform between 0.05 and 0.95), so F is the only difference between the
groups; F acts site by site, with no runs of homozygosity. 150 bp PCR-free pairs simulated at 30x,
mapped with minimap2, downsampled to each depth. pvc runs with --call-satellites, since the
simulation plants no stutter and no indels, so a SNP inside a repeat too long for pvc's repeat model
is called, not left out. pvc runs at its defaults (F = 0 for every sample),
fitted (each sample's F estimated), and with F stated as 0.9 for every sample; GATK and DeepVariant
take no F. Each sample's genotype is compared with the truth at every planted SNP: `het_called_het`
is the share of true heterozygotes called heterozygous, `homozygote_called_het` the homozygotes
called heterozygous, `nonref_concordance` the non-reference genotypes right over all non-reference
genotypes in the truth or the calls. "pvc fitted F" rows set each group's estimated F against the
F realised in the simulation.
"""

FS = [0.0, 0.25, 0.5, 0.75, 1.0]
PER_GROUP = 10
CONTIG_SOURCE, START, END = "SL4.0ch01", 5_000_001, 10_000_000
SNP_DENSITY, MIN_AF, MAX_AF = 0.001, 0.05, 0.95
COVERAGE = 30
DEPTHS = ["3x", "5x", "10x", "15x", "30x"]
SEED = 42
STATED_F = 0.9
# pvc genotypes satellite tracts too (off by default). Safe here: GPS plants SNPs only, with no
# stutter and no indels, so a long repeat holds nothing pvc's SNP/indel model would get wrong; left
# off, pvc wrote nothing at the SNPs inside them (17 positions at 30x, 2026-10-31).
CALL_SATELLITES = True
ARMS = ["pvc_defaults", "pvc_fitted", "gatk", "deepvariant"]
STATED_LABEL = f"pvc (F stated {STATED_F})"
# The slice: 100 kb, 2 samples a group, 10x, called at 5x.
SLICE = dict(end=START + 99_999, per_group=2, coverage=10, depths=["5x"])
q = shlex.quote


def _settings(ctx):
    end = SLICE["end"] if ctx.slice else END
    return dict(contig_source=CONTIG_SOURCE, start=START, end=end, fs=FS,
                per_group=SLICE["per_group"] if ctx.slice else PER_GROUP,
                coverage=SLICE["coverage"] if ctx.slice else COVERAGE, density=SNP_DENSITY,
                min_af=MIN_AF, max_af=MAX_AF, seed=SEED)


def _gps_commit():
    repo = os.path.expanduser("~/genome_population_simulator")
    out = subprocess.run(["git", "-C", repo, "rev-parse", "--short", "HEAD"], capture_output=True,
                         text=True)
    return out.stdout.strip() or "unknown"


def steps(ctx):
    cfg = _settings(ctx)
    contig = f"{cfg['contig_source']}_{cfg['start']}_{cfg['end']}"
    name = "inbreeding_sim"
    sim_dir = os.path.join(ctx.results_root, "inputs", name, "simulate")
    reference = os.path.join(sim_dir, "reference.fa")
    regions = os.path.join(sim_dir, "regions.bed")
    truth = os.path.join(sim_dir, "truth.vcf.gz")
    snvs = os.path.join(sim_dir, "snvs.tsv")
    samples = popsim.sample_names(cfg["fs"], cfg["per_group"])
    full = {s: os.path.join(ctx.results_root, "inputs", name, "map", f"{s}.bam") for s in samples}
    out = []

    # ---- inputs: the segment, the SNPs, the simulation, then the reads mapped ----------------
    def simulate(d):
        with open(os.path.join(d, "reference.fa"), "w") as fh:
            seq = subprocess.run(["samtools", "faidx", speed.reference(ctx),
                                  f"{cfg['contig_source']}:{cfg['start']}-{cfg['end']}"],
                                 capture_output=True, text=True, check=True).stdout
            fh.write(f">{contig}\n" + "".join(seq.splitlines(keepends=True)[1:]))
        with open(reference) as fh:
            seq = popsim.read_segment(fh)
        n_snvs = popsim.write_snvs(snvs, contig, seq, cfg["fs"], cfg["density"], cfg["min_af"],
                                   cfg["max_af"], cfg["seed"])
        config = os.path.join(d, "general.toml")
        popsim.write_config(config, "reference.fa", contig, len(seq), "snvs.tsv", cfg["fs"],
                            cfg["per_group"], cfg["coverage"], cfg["seed"])
        with open(regions, "w") as fh:
            fh.write(f"{contig}\t0\t{len(seq)}\n")
        gps = q(popsim.GPS_BIN)
        m = run_commands([("simulate",
                           f"{gps} validate {q(config)} && {gps} run {q(config)} "
                           f"--output {q(os.path.join(d, 'gps'))} "
                           f"&& samtools faidx {q(reference)} "
                           f"&& samtools dict {q(reference)} -o {q(os.path.join(d, 'reference.dict'))} "
                           f"&& bcftools view --no-version {q(os.path.join(d, 'gps', 'outputs', 'truth.vcf'))} "
                           f"-Oz -o {q(truth)} && bcftools index -t {q(truth)}")],
                         os.path.join(d, "logs"))
        return {"measurement": m, "meta": {"snvs": n_snvs, "contig": contig}}

    out.append(Step(f"{name}/simulate", "input", simulate,
                    lambda: {"settings": cfg, "simulator": file_identity(popsim.GPS_BIN),
                             "simulator_commit": _gps_commit(), "popsim": 1}))

    def map_reads(d):
        fastq = os.path.join(sim_dir, "gps", "outputs", "fastq", "reads")
        cmds = []
        for s in samples:
            r1, r2 = (os.path.join(fastq, f"{s}_R{i}.fastq.gz") for i in (1, 2))
            cmds.append((f"map.{s}",
                         f"minimap2 -ax sr -t {ctx.threads} -R '@RG\\tID:{s}\\tSM:{s}\\tLB:{s}\\tPL:ILLUMINA' "
                         f"{q(reference)} {q(r1)} {q(r2)} | samtools sort -@ 4 -o {q(full[s])} - "
                         f"&& samtools index {q(full[s])}"))
        m = run_commands(cmds, os.path.join(d, "logs"))
        return {"measurement": m}

    out.append(Step(f"{name}/map", "input", map_reads,
                    lambda: {"reference": file_identity(reference + ".fai"),
                             "minimap2": ctx.version("minimap2")}))

    depths = SLICE["depths"] if ctx.slice else DEPTHS
    ladder, alignment = inputs.depth_ladder(ctx, name, reference, full, regions, depths)
    out.extend(ladder)

    # ---- pvc's repeat catalogue, every commit -------------------------------------------------
    catalog = os.path.join(ctx.pvc_root(), name, "catalog", "reference.repeats.parquet")

    def make_catalog(d):
        m = run_commands([callers.pvc_catalog(ctx, reference, catalog, os.path.join(d, "catalog.log"))],
                         os.path.join(d, "logs"))
        return {"measurement": m}

    out.append(Step(f"{name}/catalog", "pvc", make_catalog,
                    lambda: {"reference": file_identity(reference + ".fai")},
                    attrs={"caller": "pvc"}))

    for depth in depths:
        base = dict(benchmark=name,
                    dataset=f"simulated tomato, {contig}, {len(samples)} samples in 5 F groups",
                    samples=len(samples), depth=depth, region_bp=arms.BedLength(regions))
        bams = {s: alignment(s, depth) for s in samples}
        prefix = f"{name}/{depth}"
        callsets = {}
        for arm in ARMS:
            step, callset = arms.arm_steps(ctx, arm, prefix, reference, catalog, bams, regions,
                                           base, joint=True, call_satellites=CALL_SATELLITES)
            out.append(step)
            callsets[arm] = callset
        step, stated_calls, stated_params = _stated_f_step(ctx, prefix, reference, catalog, bams,
                                                          regions, base, callsets["pvc_fitted"])
        out.append(step)
        scored = [(arms.ARM_LABELS[a], callsets[a], arms.ARM_KIND[a], a) for a in ARMS]
        scored.append((STATED_LABEL, arms.CallSet("pvc_fitted", step.name, stated_calls), "pvc",
                       "pvc_stated_f"))
        for label, callset, kind, arm in scored:
            out.append(_score_step(ctx, label, callset, kind, arm, truth, reference, samples, base))
        fitted_params = os.path.join(ctx.pvc_root(), callsets["pvc_fitted"].step_name,
                                     "fitted.parameters.toml")
        out.append(_f_step(ctx, prefix, fitted_params, truth, snvs, samples, base))
    return out


def _stated_f_step(ctx, prefix, reference, catalog, bams, regions, base, fitted):
    """pvc told every sample's F (STATED_F): `estimate-parameters --inbreeding`, then
    `call-from-psps`, over the fitted arm's psps (made again when they are gone). The psps record
    CALL_SATELLITES, and the commands that read them follow it."""
    step_name = f"{prefix}/pvc_stated_f"
    d_final = os.path.join(ctx.pvc_root(), step_name)
    calls = os.path.join(d_final, "calls.vcf.gz")
    params = os.path.join(d_final, "stated.parameters.toml")
    psp_dir = os.path.join(ctx.pvc_root(), fitted.step_name, "psps")
    samples = list(bams)
    row_base = dict(base, caller=STATED_LABEL, caller_version=ctx.version("pvc"))

    def action(d):
        logs = os.path.join(d, "logs")
        source = psp_dir
        if not all(os.path.exists(os.path.join(psp_dir, f"{s}.psp")) for s in samples):
            source = os.path.join(d, "psps")
            run_commands(callers.pvc_generate_psps(ctx, reference, catalog,
                                                   {k: [v] for k, v in bams.items()}, regions,
                                                   source, logs, call_satellites=CALL_SATELLITES),
                         logs, max_concurrent=ctx.threads)
        psps = [os.path.join(source, f"{s}.psp") for s in samples]
        m_fit = run_commands([callers.pvc_estimate_parameters(
            ctx, reference, catalog, psps, params, os.path.join(d, "fit.log"),
            inbreeding=STATED_F)], logs)
        m_call = run_commands([callers.pvc_call_from_psps(
            ctx, reference, catalog, psps, calls, os.path.join(d, "call.log"), parameters=params)],
            logs)
        return {"measurement": m_call,
                "rows": time_rows(row_base, m_fit, "estimate-parameters")
                + time_rows(row_base, m_call, "call-from-psps")}

    step = Step(step_name, "pvc", action,
                lambda: {"alignments": [file_identity(b) for b in bams.values()],
                         "regions": file_identity(regions), "stated_f": STATED_F,
                         "threads": ctx.threads, "call_satellites": CALL_SATELLITES},
                attrs={"caller": "pvc", "depth": base["depth"]})
    return step, calls, params


def _score_step(ctx, label, callset, kind, arm, truth, reference, samples, base):
    row_base = dict(base, caller=label,
                    caller_version=ctx.version("pvc" if kind == "pvc" else arms.ARM_TOOL[arm]))
    min_qual = arms.MIN_QUAL.get(arm)

    def action(d):
        hostjobs.skip_with(ctx, callset)
        rows = genotype_concordance.score(truth, callset.path, reference, samples,
                                          popsim.group_of, row_base,
                                          os.path.join(d, "per_sample.tsv"), min_qual=min_qual)
        return {"rows": rows}

    return Step(f"{callset.step_name}/score", "score", action,
                lambda: {"calls": file_identity(callset.path), "truth": file_identity(truth),
                         "version": 1},
                attrs={"caller": "pvc" if kind == "pvc" else arms.ARM_CALLER[arm],
                       "depth": base["depth"], "home": "pvc" if kind == "pvc" else "competitor"})


def _f_step(ctx, prefix, params, truth, snvs, samples, base):
    """Each group's F as pvc fitted it, against the F the simulation realised."""
    row_base = dict(base, caller="pvc (fitted)", caller_version=ctx.version("pvc"),
                    stage="inbreeding coefficient")

    def action(d):
        fitted = genotype_concordance.fitted_inbreeding(params)
        realised = genotype_concordance.realised_inbreeding(truth, snvs, samples, popsim.group_of)
        rows = []
        with open(os.path.join(d, "per_sample.tsv"), "w") as fh:
            fh.write("sample\tsimulated_f\tfitted_f\twarrant\n")
            for s in samples:
                f, warrant = fitted.get(s, ("", "missing"))
                fh.write(f"{s}\t{popsim.group_of(s)}\t{f}\t{warrant}\n")
        for g in sorted(realised):
            values = [fitted[s][0] for s in samples if popsim.group_of(s) == g and s in fitted]
            cls = dict(row_base, variant_class=f"F = {g:.2f}")
            rows.append(row(cls, metric="realised_f", value=realised[g][0]))
            if values:
                rows.append(row(cls, metric="fitted_f_mean", value=round(sum(values) / len(values), 4)))
                rows.append(row(cls, metric="fitted_f_min", value=round(min(values), 4)))
                rows.append(row(cls, metric="fitted_f_max", value=round(max(values), 4)))
        return {"rows": rows}

    return Step(f"{prefix}/pvc_fitted/inbreeding", "score", action,
                lambda: {"parameters": file_identity(params), "truth": file_identity(truth),
                         "version": 1},
                attrs={"caller": "pvc", "depth": base["depth"], "home": "pvc"})
