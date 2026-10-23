"""Where the human benchmarks' data come from: the data checkout if it has them, GIAB otherwise.

On the laptop the data checkout holds every file, and they are used as they are: nothing here
changes a path there, so nothing already measured reruns. On a machine without them (a server), the
suite fetches the same data from their public sources as input steps, made once:

- **GRCh38** (the analysis set with decoys, no alts) from NCBI; its index is byte-identical to the
  laptop's copy.
- **Alignments**: GIAB's 300× novoalign BAMs of HG002, HG003 and HG004, fetched through each BAM's
  index over a benchmark's regions only, given the read group the files lack, and sorted. The
  panel's depths are downsampled with the fractions recorded when the laptop's files were made
  (`data/giab_panel/fractions.tsv`); `samtools view -s` keeps a read by its name and the seed, so
  the same fraction keeps the same reads.
- **Truth**: GIAB's v4.2.1 small-variant VCF and confident regions per sample, and the HG002
  tandem-repeat benchmark v1.0.1, cut to each benchmark's regions.

The region lists this project chose (the panel's 100 regions per sample, the trio's 100, the
tandem-repeat benchmark's 50,000) are in `benchmarks/suite/data/`, so every machine uses the same
ground; the 5 Mb set's list is `benchmarks/human_genome_bottle/HG002_bench_azar_merged_1000.bed`.
"""

import csv
import os
import shlex

from .measure import run_commands
from .steps import Step, file_identity

q = shlex.quote
SUITE_DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

FTP = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab"
SAMPLE_DIRS = {
    "HG002": ("HG002_NA24385_son", "NIST_HiSeq_HG002_Homogeneity-10953946"),
    "HG003": ("HG003_NA24149_father", "NIST_HiSeq_HG003_Homogeneity-12389378"),
    "HG004": ("HG004_NA24143_mother", "NIST_HiSeq_HG004_Homogeneity-14572558"),
}
REFERENCE_URL = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/"
                 "seqs_for_alignment_pipelines.ucsc_ids/"
                 "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz")
REFERENCE_GZ_MD5 = "a056c57649f3c9964c68aead3849bbf8"
REFERENCE_NAME = "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
LOCAL_REFERENCE = os.path.expanduser(f"~/genomes/h_sapiens/gca_grch38/{REFERENCE_NAME}")
TR_TRUTH_URL = (f"{FTP}/release/AshkenazimTrio/HG002_NA24385_son/TandemRepeats_v1.0/GRCh38/"
                "HG002_GRCh38_TandemRepeats_v1.0.1.vcf.gz")


def bam_url(sample):
    family, homogeneity = SAMPLE_DIRS[sample]
    return (f"{FTP}/data/AshkenazimTrio/{family}/{homogeneity}/"
            f"NHGRI_Illumina300X_AJtrio_novoalign_bams/{sample}.GRCh38.300x.bam")


def truth_urls(sample):
    family, _ = SAMPLE_DIRS[sample]
    base = f"{FTP}/release/AshkenazimTrio/{family}/NISTv4.2.1/GRCh38/{sample}_GRCh38_1_22_v4.2.1"
    return base + "_benchmark.vcf.gz", base + "_benchmark_noinconsistent.bed"


def data_file(name):
    """A file the suite carries (benchmarks/suite/data/)."""
    return os.path.join(SUITE_DATA, name)


def _inputs(ctx, *parts):
    return os.path.join(ctx.results_root, "inputs", "giab", *parts)


# ---- the reference ------------------------------------------------------------------------------

def grch38(ctx):
    """(steps, path) of the GRCh38 analysis set: the local copy, or NCBI's, downloaded once."""
    if os.path.exists(LOCAL_REFERENCE + ".fai"):
        return [], LOCAL_REFERENCE
    path = _inputs(ctx, "reference", REFERENCE_NAME)

    def action(d):
        gz = os.path.join(d, REFERENCE_NAME + ".gz")
        m = run_commands([("download",
                           f"curl -fsSL -o {q(gz)} {q(REFERENCE_URL)} "
                           f"&& test \"$(md5sum {q(gz)} | cut -d' ' -f1)\" = {REFERENCE_GZ_MD5} "
                           f"&& gunzip {q(gz)} && samtools faidx {q(path)} "
                           f"&& samtools dict {q(path)} -o {q(path[:-4] + '.dict')}")],
                         os.path.join(d, "logs"))
        return {"measurement": m}

    return [Step("giab/reference", "input", action, {"url": REFERENCE_URL, "md5": REFERENCE_GZ_MD5})], path


# ---- truth ----------------------------------------------------------------------------------------

def whole_truth(ctx, sample):
    """(steps, vcf, bed): GIAB v4.2.1 for one sample, whole genome."""
    local = ctx.data("benchmarks", "giab", "all_bench_regions")
    vcf = os.path.join(local, "vcfs", f"{sample}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")
    bed = os.path.join(local, "bed", f"{sample}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed")
    if os.path.exists(vcf) and os.path.exists(bed):
        return [], vcf, bed
    vcf_url, bed_url = truth_urls(sample)
    d_final = _inputs(ctx, "truth", sample)
    vcf, bed = os.path.join(d_final, os.path.basename(vcf_url)), os.path.join(d_final, os.path.basename(bed_url))

    def action(d):
        m = run_commands([("download", f"curl -fsSL -o {q(vcf)} {q(vcf_url)} "
                                       f"&& curl -fsSL -o {q(vcf)}.tbi {q(vcf_url)}.tbi "
                                       f"&& curl -fsSL -o {q(bed)} {q(bed_url)}")],
                         os.path.join(d, "logs"))
        return {"measurement": m}

    return [Step(f"giab/truth/{sample}", "input", action, {"vcf": vcf_url, "bed": bed_url})], vcf, bed


def cut_truth(ctx, name, vcf, regions):
    """(step, path): a truth VCF cut to a benchmark's regions, indexed."""
    path = _inputs(ctx, "cut", name, "truth.vcf.gz")

    def action(d):
        m = run_commands([("cut", f"bcftools view --no-version -R {q(regions)} {q(vcf)} "
                                  f"| bcftools sort -Oz -o {q(path)} && bcftools index -t {q(path)}")],
                         os.path.join(d, "logs"))
        return {"measurement": m}

    return Step(f"giab/cut/{name}", "input", action,
                lambda: {"vcf": file_identity(vcf), "regions": file_identity(regions)}), path


def tr_truth(ctx):
    """(steps, path): the HG002 tandem-repeat benchmark v1.0.1, whole."""
    local = ctx.data("benchmarks", "ssr_hg002", "truth", "HG002_GRCh38_TandemRepeats_v1.0.1_50000.vcf.gz")
    if os.path.exists(local):
        return [], local
    path = _inputs(ctx, "truth", "tandem_repeats", os.path.basename(TR_TRUTH_URL))

    def action(d):
        m = run_commands([("download", f"curl -fsSL -o {q(path)} {q(TR_TRUTH_URL)} "
                                       f"&& curl -fsSL -o {q(path)}.tbi {q(TR_TRUTH_URL)}.tbi")],
                         os.path.join(d, "logs"))
        return {"measurement": m}

    return [Step("giab/truth/tandem_repeats", "input", action, {"url": TR_TRUTH_URL})], path


# ---- alignments ---------------------------------------------------------------------------------

def fetched(ctx, name, sample, regions, reference):
    """(step, path): a sample's 300× GIAB alignment over `regions`, with a read group, sorted.

    Read from GIAB's server, or, where the machine's profile names local copies of the same
    alignments (`--giab-alignments`, a path with `{sample}`; rick's are GIAB's novoalign 300× BAMs
    as CRAMs of chromosomes 1 to 22), from those: the same reads, at disk speed instead of hours of
    HTTP range requests.
    """
    path = _inputs(ctx, "alignments", name, sample, f"{sample}.300x.bam")
    rg = f"@RG\\tID:{sample}\\tSM:{sample}\\tLB:{sample}\\tPL:ILLUMINA"
    local = ctx.giab_alignments.format(sample=sample) if ctx.giab_alignments else None

    def action(d):
        source = f"-T {q(reference)} {q(local)}" if local else q(bam_url(sample))
        m = run_commands([("fetch",
                           f"samtools view -u -h -M -L {q(regions)} {source} "
                           f"| samtools addreplacerg -m overwrite_all -w -r {q(rg)} -O bam -u - "
                           f"| samtools sort -@ 4 -m 1G -T {q(os.path.join(d, 'sort'))} -o {q(path)} - "
                           f"&& samtools index {q(path)}")],
                         os.path.join(d, "logs"), cwd=d)
        return {"measurement": m}

    def fingerprint():
        source = {"local": file_identity(local)} if local else {"url": bam_url(sample)}
        return {**source, "regions": file_identity(regions)}

    return Step(f"giab/alignments/{name}/{sample}", "input", action, fingerprint), path


def downsampled(ctx, name, sample, source, reference, depth, fraction):
    """(step, path): `source` downsampled with `samtools view -s 42.<fraction>`."""
    path = _inputs(ctx, "alignments", name, sample, f"{sample}.{depth}.seed42.bam")

    def action(d):
        m = run_commands([("downsample",
                           f"samtools view -@ 4 -T {q(reference)} -b -s {42 + fraction:.12f} "
                           f"-o {q(path)} {q(source)} && samtools index {q(path)}")],
                         os.path.join(d, "logs"))
        return {"measurement": m}

    return Step(f"giab/alignments/{name}/{sample}/{depth}", "input", action,
                lambda: {"source": file_identity(source), "fraction": fraction}), path


def panel_fractions():
    """{(sample, depth): fraction} recorded when the laptop's panel files were made."""
    with open(data_file("giab_panel/fractions.tsv")) as fh:
        return {(r["sample"], r["depth"]): float(r["fraction"])
                for r in csv.DictReader(fh, delimiter="\t")}
