"""Speed and peak memory of each caller's whole pipeline on the tomato cohort (plan §4.6).

62 tomato accessions at about 3 reads a position over 80 regions of 100 kb (8 Mb of SL4.0), called
at 1, 8 and 62 samples so the table shows how each caller grows with the cohort. The one-human-sample
figures come from the HG002 5 Mb benchmark (`hg002_5mb`), whose arms are timed the same way.

Each caller runs as its documentation tells a user to (plan §3.4), at 6 threads: one process per
sample for the one-sample stages (pvc generate-psps, GATK HaplotypeCaller GVCFs), one threaded
process for the cohort stages, freebayes as 6 processes over blocks of the regions. GATK's cohort
stages are single-threaded and run as GATK's joint-calling workflow runs them, one process per
shard of the regions.

GATK's GVCFs are made once, for the largest cohort (lib/arms.py, `gvcf_store_step`): a sample's
GVCF is the same file in every cohort, and making it again for each would add hours that say
nothing new. A cohort's per-sample stage figures are derived from those runs and marked so.
"""

import os
import subprocess

from lib import arms, readfloor
from lib.steps import Step, file_identity

TITLE = "Speed and peak memory — tomato cohort"
DESCRIPTION = """
Tomato (S. lycopersicum, SL4.0). 62 accessions at about 3 reads a position, over 80 regions of
100 kb (8 Mb); the first 1, 8 and all 62 accessions in name order. (The data directory holds 63
files; SRR7279543 is left out because its CRAM is truncated.) Each caller's whole pipeline at
6 threads, inside a Linux virtual machine given 6 cores and 32 GB on the machine named below:
wall_seconds is elapsed time, cpu_seconds user plus system time summed over every process,
peak_rss_mb_all_processes the largest total resident memory of the processes running at once (the
real cost of a stage that runs one process per sample), peak_rss_mb_largest_process the largest
single process. pvc is timed both ways it can be run: call-from-alignments in one step, and
generate-psps, estimate-parameters and call-from-psps ("whole pipeline" is their sum). GATK makes
each accession's GVCF once, in a run of the largest cohort, and each cohort's "HaplotypeCaller
GVCFs" figures are derived from those runs (marked †); it then joins the GVCFs (GenomicsDBImport on
x86-64, CombineGVCFs elsewhere) and genotypes them, one process per shard of the regions. The
repeat catalogue (repeat-catalog), built once per reference, is timed separately under
"repeat_catalog".
"""

SIZES = [1, 8, 62]
# The slice's tomato ground: the first region, cut to this many bases (tomato1's regions are
# 100 kb, so there it is the whole first region).
SLICE_REGION_BP = 100_000

# Files of the tomato data directory left out of every tomato benchmark, and why.
EXCLUDED = {
    "SRR7279543.p1.bench.cram": "truncated: no end-of-file block, and its reads stop at "
                                "SL4.0ch08:54,312,875, so it has none in the regions on chromosomes "
                                "9 to 12 (found 2026-10-22; the owner chose to drop it)",
}
ARMS = ["pvc_defaults", "pvc_fitted", "gatk", "freebayes"]
LOCAL_REFERENCE = os.path.expanduser("~/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa")
TOMATO1 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "tomato",
                       "tomato1_cohort.tsv")


def reference(ctx):
    """SL4.0: the machine profile's copy (`--tomato-reference`), else the laptop's."""
    return ctx.tomato_reference or LOCAL_REFERENCE


def sizes(ctx):
    return [2] if ctx.slice else (ctx.tomato_sizes or SIZES)


def cohort(ctx):
    """The tomato cohort, sample name -> alignment file(s), in the order samples are added.

    On the laptop, tomato1's CRAMs (sliced to the regions). On a machine whose profile names a
    directory of tomato CRAMs (`--tomato-crams`; kimura's whole-genome collection), the same 62
    accessions, found by sample name, followed by further accessions of the same project
    (PRJNA454805) and then the others, by name, for sample counts above 62.
    """
    if ctx.tomato_crams:
        return _server_cohort(ctx)
    crams = ctx.data("benchmarks", "tomato1", "crams")
    out = {}
    for f in sorted(os.listdir(crams)):
        if not f.endswith(".cram") or f in EXCLUDED:
            continue
        path = os.path.join(crams, f)
        header = subprocess.run(["samtools", "view", "-H", path], capture_output=True, text=True,
                                check=True).stdout
        sample = next(field[3:] for line in header.splitlines() if line.startswith("@RG")
                      for field in line.split("\t") if field.startswith("SM:"))
        out[sample] = path
    return out


def _server_cohort(ctx):
    """tomato1's accessions first, then the rest, from every CRAM under ctx.tomato_crams.

    Every CRAM's header is read once (its read groups' SM and LB), and the table is kept in
    `inputs/tomato/catalogue.tsv`; delete it to read the directory again.
    """
    import csv

    table = os.path.join(ctx.results_root, "inputs", "tomato", "catalogue.tsv")
    if not os.path.exists(table):
        os.makedirs(os.path.dirname(table), exist_ok=True)
        rows = []
        for root, _, files in os.walk(ctx.tomato_crams):
            for f in sorted(files):
                if not f.endswith(".cram") or f in EXCLUDED:
                    continue
                path = os.path.join(root, f)
                header = subprocess.run(["samtools", "view", "-H", path], capture_output=True,
                                        text=True).stdout
                for line in header.splitlines():
                    if line.startswith("@RG"):
                        tags = dict(t.split(":", 1) for t in line.split("\t")[1:] if ":" in t)
                        rows.append((tags.get("SM", ""), tags.get("LB", "").split("_")[0], path))
                        break
        with open(table + ".tmp", "w") as fh:
            fh.writelines(f"{s}\t{p}\t{f}\n" for s, p, f in sorted(set(rows)))
        os.replace(table + ".tmp", table)
    files, project = {}, {}
    with open(table) as fh:
        for line in fh:
            sample, proj, path = line.rstrip("\n").split("\t")
            if sample:
                files.setdefault(sample, []).append(path)
                project[sample] = proj
    with open(TOMATO1) as fh:
        wanted = [r["sample"] for r in csv.DictReader(fh, delimiter="\t")]
    order = [s for s in wanted if s in files]
    missing = [s for s in wanted if s not in files]
    if missing:
        ctx.log(f"tomato cohort: {len(missing)} of tomato1's 62 accessions are not in "
                f"{ctx.tomato_crams}; others take their places: {', '.join(missing[:10])}")
    rest = sorted((s for s in files if s not in order),
                  key=lambda s: (project.get(s) != "PRJNA454805", s))
    return {s: files[s] for s in order + rest}


def regions_step(ctx):
    # tomato1's 80 regions of 100 kb (8 Mb), the only ground its sliced CRAMs cover; a machine with
    # whole-genome CRAMs can set larger regions in its profile (BENCH_TOMATO_REGIONS; kimura: all
    # of chromosome 1).
    src = ctx.tomato_regions or ctx.data("benchmarks", "tomato1", "regions.bed")
    path = os.path.join(ctx.results_root, "inputs", "tomato", "regions", "regions.bed")

    def action(d):
        with open(src) as fh:
            lines = [l for l in fh if l.strip() and not l.startswith(("#", "track"))]
        if ctx.slice:
            # The slice is the first region, cut to SLICE_REGION_BP: a profile's regions can be one
            # whole chromosome (kimura's and rick's are SL4.0ch01, 90.9 Mb), and a slice of that is
            # no slice.
            contig, start, end = lines[0].split("\t")[:3]
            start, end = int(start), int(end)
            lines = [f"{contig}\t{start}\t{min(end, start + SLICE_REGION_BP)}\n"]
        with open(path, "w") as fh:
            fh.writelines(lines)
        return {}

    return Step("tomato/regions", "input", action,
                {"bed": file_identity(src), "slice_regions": f"1:{SLICE_REGION_BP}" if ctx.slice else "all"}), path


def steps(ctx):
    REFERENCE = reference(ctx)
    out, catalog, _ = arms.reference_steps(ctx, "SL4.0", REFERENCE)
    out = [s for s in out if not s.name.startswith("refs/")]  # no truth here, so no vcfeval copy
    step, regions = regions_step(ctx)
    out.append(step)
    samples = cohort(ctx)
    names = list(samples)
    for size in sizes(ctx):
        if size > len(names):
            ctx.log(f"speed: {size} accessions asked for, the cohort has {len(names)}; skipped")
    run_sizes = [size for size in sizes(ctx) if size <= len(names)]
    # The cohorts are nested (the first 1, 2, 5, ... accessions), so GATK makes each accession's
    # GVCF once, for the largest, and every cohort's GATK step joins its share of them.
    gvcfs = None
    if run_sizes:
        largest = {s: samples[s] for s in names[:max(run_sizes)]}
        gvcf_step, gvcfs = arms.gvcf_store_step(ctx, "speed/tomato", REFERENCE, largest, regions)
    for size in run_sizes:
        chosen = {s: samples[s] for s in names[:size]}
        base = dict(benchmark="speed", dataset="tomato cohort, SL4.0, 80 regions of 100 kb",
                    samples=size, depth="about 3x", region_bp=arms.BedLength(regions))
        for arm in ARMS:
            if arm == "gatk" and gvcf_step not in out:
                out.append(gvcf_step)
            step, _ = arms.arm_steps(ctx, arm, f"speed/tomato_{size}", REFERENCE, catalog, chosen,
                                     regions, base, joint=True, gvcfs=gvcfs)
            out.append(step)
        out.extend(readfloor.floor_steps(ctx, f"speed/tomato_{size}", REFERENCE, chosen, regions,
                                         base))
        out.extend(readfloor.startup_steps(ctx, f"speed/tomato_{size}", REFERENCE, catalog,
                                           chosen, regions, base))
    return out
