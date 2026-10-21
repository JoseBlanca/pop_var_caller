"""GIAB HG002, all of chromosome 20: accuracy, and the human speed figures on a realistic size.

The alignments are GIAB's public HG002 60× novoalign BAM (`HG002.GRCh38.60x.1.bam`, the same
alignment the per-sample panel's 300× files come from), of which the suite fetches chromosome 20
alone over HTTPS, through the BAM's index, and adds the read group the file lacks. The truth is
GIAB's HG002 v4.2.1 small-variant benchmark, already in the data checkout
(`benchmarks/giab/all_bench_regions/`, the same bytes GIAB publishes), cut to chromosome 20.

The callers run over the whole chromosome, as a user would; the scoring is restricted to GIAB's
confident regions. Depths 5×, 10× and 30× are downsampled by the suite (seed 42); "full" is the
file's own depth, about 60×.
"""

import os
import shlex

from lib import arms, inputs, readfloor
from lib.measure import run_commands
from lib.steps import Step, file_identity

TITLE = "GIAB HG002, chromosome 20 — SNPs and indels over a whole chromosome"
DESCRIPTION = """
Human (GRCh38). GIAB HG002, all of chromosome 20 (64 Mb), from GIAB's public 60× Illumina
alignments (novoalign, 2×148 bp); 5×, 10× and 30× downsampled by the suite (seed 42), "full" is
about 60×. Each caller runs over the whole chromosome. Truth: GIAB v4.2.1 small variants, scored by
RTG vcfeval inside GIAB's confident regions of chromosome 20. The time and memory rows are the
human speed figures: one sample, a whole chromosome, at 6 threads.
"""

URL = ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
       "HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/"
       "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam")
CONTIG = "chr20"
SLICE = ("chr20", 10_000_000, 10_500_000)
DEPTHS = ["5x", "10x", "30x", "full"]
ARMS = ["pvc_defaults", "pvc_fitted", "gatk", "freebayes"]
READ_GROUP = r"@RG\tID:HG002\tSM:HG002\tLB:HG002\tPL:ILLUMINA"
q = shlex.quote


def steps(ctx):
    giab = ctx.data("benchmarks", "giab", "all_bench_regions")
    whole_truth = os.path.join(giab, "vcfs", "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")
    whole_bed = os.path.join(giab, "bed", "HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed")
    region = f"{SLICE[0]}:{SLICE[1] + 1}-{SLICE[2]}" if ctx.slice else CONTIG
    start, end = (SLICE[1], SLICE[2]) if ctx.slice else (0, None)
    depths = ["10x"] if ctx.slice else DEPTHS
    inputs_dir = os.path.join(ctx.results_root, "inputs", "giab_chr20")
    out, catalog, sdf = arms.reference_steps(ctx, "GRCh38", arms.GRCH38)

    # ---- the alignment: chromosome 20 fetched through the index, with a read group added -------
    # GIAB's header says SO:unsorted (as its other HG002 downloads do, benchmarks/ssr_hg002/README),
    # which pvc refuses; sorting the fetched reads makes the header true.
    alignment = os.path.join(inputs_dir, "alignment", "HG002.chr20.60x.bam")

    def fetch(d):
        cmd = (f"samtools view -u -h {q(URL)} {region} "
               f"| samtools addreplacerg -m overwrite_all -r {q(READ_GROUP)} -O bam -u - "
               f"| samtools sort -@ 4 -m 1G -T {q(os.path.join(d, 'sort'))} -o {q(alignment)} - "
               f"&& samtools index {q(alignment)}")
        m = run_commands([("fetch", cmd)], os.path.join(d, "logs"), cwd=d)
        return {"measurement": m}

    out.append(Step("giab_chr20/alignment", "input", fetch, {"url": URL, "region": region, "sorted": True}))

    # ---- the regions the callers run over, and the truth and confident regions cut to them ----
    calling = os.path.join(inputs_dir, "truth", "calling.bed")
    confident = os.path.join(inputs_dir, "truth", "confident.bed")
    truth = os.path.join(inputs_dir, "truth", "truth.vcf.gz")

    def cut_truth(d):
        fai = arms.GRCH38 + ".fai"
        contig_len = next(int(l.split()[1]) for l in open(fai) if l.split()[0] == CONTIG)
        stop = end if end is not None else contig_len
        with open(calling, "w") as fh:
            fh.write(f"{CONTIG}\t{start}\t{stop}\n")
        with open(whole_bed) as src, open(confident, "w") as dst:
            for line in src:
                f = line.split()
                if f[0] != CONTIG:
                    continue
                s, e = max(int(f[1]), start), min(int(f[2]), stop)
                if s < e:
                    dst.write(f"{CONTIG}\t{s}\t{e}\n")
        m = run_commands([("truth", f"bcftools view --no-version -r {region} {q(whole_truth)} "
                                    f"-Oz -o {q(truth)} && bcftools index -t {q(truth)}")],
                         os.path.join(d, "logs"))
        return {"measurement": m}

    out.append(Step("giab_chr20/truth", "input", cut_truth,
                    {"truth": file_identity(whole_truth), "bed": file_identity(whole_bed),
                     "region": region}))

    ladder, aligned = inputs.depth_ladder(ctx, "giab_chr20", arms.GRCH38, {"HG002": alignment},
                                          calling, depths)
    out.extend(ladder)

    for depth in depths:
        base = dict(benchmark="giab_chr20",
                    dataset="GIAB HG002, chromosome 20" + (" (slice, 500 kb)" if ctx.slice else ""),
                    samples=1, depth=depth, region_bp=arms.BedLength(calling))
        for arm in ARMS:
            step, callset = arms.arm_steps(ctx, arm, f"giab_chr20/{depth}", arms.GRCH38, catalog,
                                           {"HG002": aligned("HG002", depth)}, calling, base)
            out.append(step)
            out.append(arms.score_step(ctx, callset, base, truth, confident, sdf,
                                       call_sample="HG002", truth_sample="HG002"))
        out.extend(readfloor.floor_steps(ctx, f"giab_chr20/{depth}", arms.GRCH38,
                                         {"HG002": aligned("HG002", depth)}, calling, base))
        out.extend(readfloor.startup_steps(ctx, f"giab_chr20/{depth}", arms.GRCH38, catalog,
                                           {"HG002": aligned("HG002", depth)}, calling, base))
    return out
