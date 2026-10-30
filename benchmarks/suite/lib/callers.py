"""How each caller is run, so that every benchmark runs it the same way.

Each function returns shell commands, as (label, command) pairs for measure.run_commands, and the
caller of the function decides how many run at once. The parallelism follows each caller's own
documentation (plan §3.4):

- one-sample tools run one process per sample, as many at a time as the suite has cores: GATK
  HaplotypeCaller making GVCFs, and pvc `generate-psps`;
- cohort tools are one process using every core: pvc `estimate-parameters`, `call-from-psps`,
  `call-from-alignments`;
- GATK's cohort stages (GenomicsDBImport or CombineGVCFs, then GenotypeGVCFs) are single-threaded,
  so they run as GATK's joint-calling workflow runs them: one process per shard of the regions,
  as many at a time as the suite has cores, and the shards' VCFs joined;
- freebayes has no threads, so it runs one process per block of regions, as `freebayes-parallel`
  does, and the blocks' VCFs are joined.
"""

import os
import platform
import shlex

q = shlex.quote


# ---- pvc --------------------------------------------------------------------------------------

def pvc_catalog(ctx, reference, out, log):
    return ("repeat_catalog",
            f"{q(ctx.pvc_bin)} repeat-catalog --reference {q(reference)} --output {q(out)} "
            f"--log {q(log)} --force")


def pvc_call_from_alignments(ctx, reference, catalog, alignments, regions, out, log,
                             parameters=None):
    alignment_args = " ".join(f"--alignment {q(a)}" for a in alignments)
    source = f"--parameters {q(parameters)}" if parameters else "--defaults"
    regions_arg = f"--regions {q(regions)}" if regions else ""
    return ("call_from_alignments",
            f"{q(ctx.pvc_bin)} call-from-alignments --reference {q(reference)} "
            f"--catalog {q(catalog)} {alignment_args} {regions_arg} {source} "
            f"--threads {ctx.threads} --output {q(out)} --log {q(log)} --force")


def pvc_generate_psps(ctx, reference, catalog, alignments_by_sample, regions, out_dir, log_dir):
    """One command per sample. `alignments_by_sample` maps sample name to its alignment files."""
    cmds = []
    regions_arg = f"--regions {q(regions)}" if regions else ""
    for sample, files in alignments_by_sample.items():
        alignment_args = " ".join(f"--alignment {q(a)}" for a in files)
        cmds.append((f"generate_psps.{sample}",
                     f"{q(ctx.pvc_bin)} generate-psps --reference {q(reference)} "
                     f"--catalog {q(catalog)} {alignment_args} {regions_arg} "
                     f"--output-dir {q(out_dir)} --force"))
    return cmds


def pvc_estimate_parameters(ctx, reference, catalog, psps, out, log):
    psp_args = " ".join(f"--psp {q(p)}" for p in psps)
    return ("estimate_parameters",
            f"{q(ctx.pvc_bin)} estimate-parameters --reference {q(reference)} "
            f"--catalog {q(catalog)} {psp_args} --threads {ctx.threads} --output {q(out)} "
            f"--log {q(log)} --force")


def pvc_call_from_psps(ctx, reference, catalog, psps, out, log, parameters=None):
    psp_args = " ".join(f"--psp {q(p)}" for p in psps)
    source = f"--parameters {q(parameters)}" if parameters else "--defaults"
    return ("call_from_psps",
            f"{q(ctx.pvc_bin)} call-from-psps --reference {q(reference)} --catalog {q(catalog)} "
            f"{psp_args} {source} --threads {ctx.threads} --output {q(out)} --log {q(log)} "
            f"--force")


# ---- GATK -------------------------------------------------------------------------------------

# The Java heap of each process in a stage that runs one process per core: 30 at once on a
# 30-thread run stay within 60 GB. HaplotypeCaller making a GVCF, and each shard's join.
SHARD_HEAP = "2g"


def _gatk(ctx, tool, args, heap=None):
    return (f"{q(ctx.gatk)} --java-options {q('-Xmx' + (heap or ctx.java_heap))} {tool} {args} "
            f"--tmp-dir {q(os.path.join(ctx.results_root, 'java_tmp'))}")


def gatk_haplotypecaller(ctx, reference, alignment, intervals, out, threads):
    """One sample straight to a VCF, using `threads` pair-HMM threads (GATK's native threading)."""
    return ("haplotypecaller",
            _gatk(ctx, "HaplotypeCaller",
                  f"-R {q(reference)} -I {q(alignment)} -L {q(intervals)} -O {q(out)} "
                  f"--native-pair-hmm-threads {threads}"))


def gatk_gvcfs(ctx, reference, alignment_by_sample, intervals, out_dir):
    """One single-threaded HaplotypeCaller process per sample, each writing `<sample>.g.vcf.gz`.

    A sample's value is one alignment file or a list of them (several runs of one sample).
    """
    cmds = []
    for sample, alignment in alignment_by_sample.items():
        out = os.path.join(out_dir, f"{sample}.g.vcf.gz")
        files = alignment if isinstance(alignment, list) else [alignment]
        inputs = " ".join(f"-I {q(f)}" for f in files)
        cmds.append((f"haplotypecaller.{sample}",
                     _gatk(ctx, "HaplotypeCaller",
                           f"-R {q(reference)} {inputs} -L {q(intervals)} -O {q(out)} "
                           f"-ERC GVCF --native-pair-hmm-threads 1", heap=SHARD_HEAP)))
    return cmds


def gatk_joint_route():
    """How GATK joins the samples' GVCFs on this machine.

    GATK's Best Practices name GenomicsDBImport for cohorts, but its native library exists only for
    x86-64 ("can't load AMD 64 .so on a AARCH64 platform" on the arm64 laptop). There the route is
    CombineGVCFs, GATK's other supported one. The stage rows carry the route's name, so a table
    says which one each machine ran.
    """
    return "GenomicsDBImport" if platform.machine() in ("x86_64", "AMD64") else "CombineGVCFs"


def gatk_join_shards(ctx, reference, gvcfs, shard_beds, work_dir, route):
    """One process per shard of the regions joining every sample's GVCF over that shard.

    GATK's joint-calling workflow scatters the import and the genotyping over intervals this way;
    each of these commands is single-threaded. Returns (commands, one GenotypeGVCFs input per
    shard): a GenomicsDB workspace (`gendb://`) or a combined GVCF.
    """
    vs = " ".join(f"-V {q(g)}" for g in gvcfs)
    cmds, sources = [], []
    for bed in shard_beds:
        name = os.path.splitext(os.path.basename(bed))[0]
        if route == "GenomicsDBImport":
            workspace = os.path.join(work_dir, f"{name}.genomicsdb")
            # --merge-input-intervals: one partition over the shard instead of one per interval;
            # there are no GVCF records between the intervals, which HaplotypeCaller kept to.
            cmds.append((f"genomicsdbimport.{name}",
                         _gatk(ctx, "GenomicsDBImport",
                               f"{vs} -L {q(bed)} --merge-input-intervals "
                               f"--genomicsdb-workspace-path {q(workspace)} --batch-size 50",
                               heap=SHARD_HEAP)))
            sources.append(f"gendb://{workspace}")
        else:
            combined = os.path.join(work_dir, f"{name}.g.vcf.gz")
            cmds.append((f"combinegvcfs.{name}",
                         _gatk(ctx, "CombineGVCFs",
                               f"-R {q(reference)} {vs} -L {q(bed)} -O {q(combined)}",
                               heap=SHARD_HEAP)))
            sources.append(combined)
    return cmds, sources


def gatk_genotype_shards(ctx, reference, sources, shard_beds, out_dir):
    """GenotypeGVCFs, one process per shard, each writing `<shard>.vcf.gz`."""
    cmds = []
    for source, bed in zip(sources, shard_beds):
        name = os.path.splitext(os.path.basename(bed))[0]
        out = os.path.join(out_dir, f"{name}.vcf.gz")
        cmds.append((f"genotypegvcfs.{name}",
                     _gatk(ctx, "GenotypeGVCFs",
                           f"-R {q(reference)} -V {q(source)} -L {q(bed)} -O {q(out)}",
                           heap=SHARD_HEAP)))
    return cmds


def gatk_gather(shard_vcfs, out):
    """The shards' VCFs joined into one sorted, indexed VCF (GATK's GatherVcfs).

    Sorted, not concatenated in the shards' order: GATK writes each shard in the reference's order,
    so where the regions put one chromosome in two shards with others between them, a plain
    concatenation stops ("the chromosome block chr18 is not contiguous") and leaves a partial file.
    """
    inputs = " ".join(q(v) for v in shard_vcfs)
    return ("gather_shards",
            f"bcftools concat --no-version -a -Ou {inputs} | bcftools sort -Oz -o {q(out)} "
            f"&& bcftools index -t {q(out)}")


# ---- freebayes --------------------------------------------------------------------------------

def split_regions(regions_bed, out_dir, blocks):
    """Split a BED into about `blocks` files of similar total length, keeping intervals whole.

    Intervals longer than 200 kb are cut into 200 kb pieces first, so one long interval cannot
    become one slow block.
    """
    pieces = []
    with open(regions_bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = line.split()[:3]
            start, end = int(start), int(end)
            while start < end:
                stop = min(end, start + 200_000)
                pieces.append((chrom, start, stop))
                start = stop
    total = sum(e - s for _, s, e in pieces)
    target = max(1, total // max(1, blocks))
    os.makedirs(out_dir, exist_ok=True)
    files, current, size = [], [], 0
    for piece in pieces:
        current.append(piece)
        size += piece[2] - piece[1]
        if size >= target:
            files.append(current)
            current, size = [], 0
    if current:
        files.append(current)
    paths = []
    for i, block in enumerate(files):
        path = os.path.join(out_dir, f"block{i:04d}.bed")
        with open(path, "w") as fh:
            for chrom, s, e in block:
                fh.write(f"{chrom}\t{s}\t{e}\n")
        paths.append(path)
    return paths


def freebayes_blocks(reference, alignments, block_beds, out_dir):
    """One freebayes process per block of regions, all samples together."""
    alignment_args = " ".join(q(a) for a in alignments)
    cmds = []
    for bed in block_beds:
        name = os.path.splitext(os.path.basename(bed))[0]
        out = os.path.join(out_dir, f"{name}.vcf")
        cmds.append((f"freebayes.{name}",
                     f"freebayes -f {q(reference)} -t {q(bed)} --ploidy 2 {alignment_args} "
                     f"> {q(out)}"))
    return cmds


def freebayes_join(block_vcfs, out):
    """Join the blocks' VCFs into one sorted, compressed, indexed VCF."""
    inputs = " ".join(q(v) for v in block_vcfs)
    return ("freebayes_join",
            f"bcftools concat --no-version {inputs} | bcftools sort -Oz -o {q(out)} "
            f"&& bcftools index -t {q(out)}")


# ---- inputs -----------------------------------------------------------------------------------

def downsample(reference, alignment, fraction, out, seed=42):
    """`samtools view -s`, the way the GIAB panel's depths were made, to a BAM with an index."""
    sub = f"{seed + fraction:.10f}" if fraction < 1 else None
    sample_arg = f"-s {sub}" if sub else ""
    return ("downsample",
            f"samtools view -@ 4 -T {q(reference)} -b {sample_arg} -o {q(out)} {q(alignment)} "
            f"&& samtools index {q(out)}")
