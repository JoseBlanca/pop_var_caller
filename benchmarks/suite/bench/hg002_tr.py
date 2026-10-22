"""HG002 tandem-repeat benchmark: STR genotypes against HipSTR and GATK (plan §4.4).

GIAB's HG002 tandem-repeat benchmark v1.0.1, 50,000 regions, at seven depths. Each caller's STR
genotypes are judged by `benchmarks/lib/tract_qual_experiment.py`, which rebuilds both copies of the
DNA over each repeat tract from the caller's records and from the truth's, and compares them; and
the small variants in the same regions by vcfeval.

**The repeat tracts are fixed at the first run.** The scorer needs a list of tracts, and HipSTR a
list of loci with their motif. Both come from pvc's typing of the regions (`pvc dev type-regions`),
made once into `inputs/` so that every pvc commit and every competitor is scored over the same
tracts. A later change to pvc's typing does not move the tract list; delete
`inputs/hg002_tr/tracts` to remake it, which also means rerunning HipSTR.
"""

import os
import shlex

from lib import arms, hostjobs, inputs, vcfeval
from lib.measure import check_alignments, run_commands, warm
from lib.results import row, time_rows
from lib.steps import Step, file_identity

TITLE = "HG002 tandem repeats — STR genotypes, one human sample"
DESCRIPTION = """
Human (GRCh38). GIAB HG002 tandem-repeat benchmark v1.0.1: 50,000 regions (6.09 Mb), 36,497 truth
records, called alone at 5× to 50× (downsampled from the 288× source, seed 42) and at the full
depth ("300x"). STR rows: per repeat tract, split into homopolymers and repeats of period 2 or more,
how many tracts the truth calls variant, how many the caller also calls, and how often the caller's
genotype is right — compared as the two rebuilt DNA sequences over the tract, and as their repeat
lengths. Small-variant rows: RTG vcfeval inside the regions. HipSTR runs as it was run on this
benchmark before (`--use-unpaired --min-reads 5`) over pvc's typed repeat loci; GATK HaplotypeCaller
with 6 pair-HMM threads. freebayes is not run here.
"""

DEPTHS = ["5x", "10x", "15x", "20x", "30x", "50x", "300x"]
ARMS = ["pvc_defaults", "pvc_fitted", "gatk", "hipstr", "deepvariant"]
GENOTYPE_METRICS = ["tracts_truth_calls", "tracts_both_call", "genotype_right", "no_call",
                    "genotype_accuracy", "right_as_lengths", "length_accuracy",
                    "called_hom_truth_het", "called_het_truth_hom", "wrong_some_other_way"]
HIPSTR_ARGS = "--use-unpaired --min-reads 5"
HIPSTR_QUAL = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib",
                           "hipstr_qual.py")
q = shlex.quote


def steps(ctx):
    data = ctx.data("benchmarks", "ssr_hg002")
    truth = os.path.join(data, "truth", "HG002_GRCh38_TandemRepeats_v1.0.1_50000.vcf.gz")
    tier = os.path.join(data, "regions", "HG002_GRCh38_TandemRepeats_v1.0.1_Tier_50000.bed")
    depths = ["10x"] if ctx.slice else DEPTHS
    out, catalog, sdf = arms.reference_steps(ctx, "GRCh38", arms.GRCH38)

    # The tier BED is unsorted; sort it, then take the slice's first 200 regions.
    inputs_dir = os.path.join(ctx.results_root, "inputs", "hg002_tr")
    regions = os.path.join(inputs_dir, "regions", "regions.bed")

    def make_regions(step_dir):
        keep = "| awk 'NR<=200' " if ctx.slice else ""
        run_commands([("sort", f"sort -k1,1V -k2,2n {q(tier)} | cut -f1-3 {keep}> {q(regions)}")],
                     os.path.join(step_dir, "logs"))
        return {}

    out.append(Step("hg002_tr/regions", "input", make_regions,
                    {"bed": file_identity(tier), "slice_regions": 200 if ctx.slice else "all"}))

    tracts_dir = os.path.join(inputs_dir, "tracts")
    tracts = os.path.join(tracts_dir, "tracts.bed")
    hipstr_regions = os.path.join(tracts_dir, "hipstr_regions.bed")

    def make_tracts(step_dir):
        # `pvc dev type-regions` reads the catalogue from beside the reference, which is read-only
        # here, so it gets a reference made of links with the catalogue beside them.
        view = os.path.join(step_dir, "reference")
        os.makedirs(view)
        ref = os.path.join(view, "GRCh38.fna")
        os.symlink(arms.GRCH38, ref)
        os.symlink(arms.GRCH38 + ".fai", ref + ".fai")
        os.symlink(catalog, ref + ".repeats.parquet")
        typed = os.path.join(step_dir, "typed_regions.tsv")
        m = run_commands([("type_regions",
                           f"{q(ctx.pvc_bin)} dev type-regions --reference {q(ref)} "
                           f"--regions {q(regions)} --output {q(typed)}")],
                         os.path.join(step_dir, "logs"))
        n = 0
        with open(typed) as fh, open(tracts, "w") as tr, open(hipstr_regions, "w") as hs:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 7 or f[3] != "ssr_locus":
                    continue
                chrom, start, end, motif, period, copies = f[0], int(f[1]), int(f[2]), f[4], f[5], f[6]
                tr.write(f"{chrom}\t{start - 1}\t{end}\t{period}\n")
                hs.write(f"{chrom}\t{start}\t{end}\t{period}\t{copies}\t{chrom}_{start}_{motif}\n")
                n += 1
        return {"measurement": m, "meta": {"tracts": n, "pvc": ctx.commit}}

    out.append(Step("hg002_tr/tracts", "input", make_tracts,
                    lambda: {"regions": file_identity(regions)}))

    for depth in depths:
        bam = os.path.join(data, "bam", depth, f"HG002_TR_v1.0.1_Tier_{depth}.bam")
        base = dict(benchmark="hg002_tr", dataset="GIAB HG002 tandem-repeat benchmark v1.0.1",
                    samples=1, depth=depth, region_bp=arms.BedLength(regions))
        for arm in ARMS:
            if arm == "hipstr":
                step, callset = _hipstr_step(ctx, depth, bam, hipstr_regions, base)
            else:
                step, callset = arms.arm_steps(ctx, arm, f"hg002_tr/{depth}", arms.GRCH38,
                                               catalog, {"HG002": bam}, regions, base)
            out.append(step)
            out.append(_str_score_step(ctx, callset, base, truth, regions, tracts, depth))
            out.append(arms.score_step(ctx, callset, base, truth, regions, sdf,
                                       call_sample="HG002", truth_sample="HG002"))
    return out


def _hipstr_step(ctx, depth, bam, hipstr_regions, base):
    name = f"hg002_tr/{depth}/hipstr"
    step_dir = os.path.join(ctx.results_root, "competitors", name)
    calls = os.path.join(step_dir, "calls.vcf.gz")
    row_base = dict(base, caller="HipSTR", caller_version=ctx.version("hipstr"))

    def action(d):
        check_alignments([bam, arms.GRCH38])
        warm([bam, arms.GRCH38])
        raw = os.path.join(d, "hipstr.vcf.gz")
        cmd = (f"HipSTR --bams {q(bam)} --bam-samps HG002 --bam-libs HG002 "
               f"--fasta {q(arms.GRCH38)} --regions {q(hipstr_regions)} --str-vcf {q(raw)} "
               f"--log {q(os.path.join(d, 'hipstr.log'))} {HIPSTR_ARGS}")
        m = run_commands([("hipstr", cmd)], os.path.join(d, "logs"))
        # HipSTR leaves QUAL empty; both scorers need it (lib/hipstr_qual.py). Not timed.
        plain = os.path.join(d, "calls.vcf")
        run_commands([("qual", f"python3 {q(HIPSTR_QUAL)} {q(raw)} {q(plain)} && "
                               f"bgzip -f {q(plain)} && tabix -p vcf {q(calls)}")],
                     os.path.join(d, "logs"))
        return {"measurement": m, "rows": time_rows(row_base, m, "HipSTR")}

    fingerprint = lambda: {"tool": ctx.version("hipstr"), "bam": file_identity(bam),
                           "regions": file_identity(hipstr_regions), "args": HIPSTR_ARGS,
                           "qual": file_identity(HIPSTR_QUAL)}
    step = Step(name, "competitor", action, fingerprint, attrs={"caller": "hipstr", "depth": depth})
    return step, arms.CallSet("hipstr", name, calls)


def _str_score_step(ctx, callset, base, truth, regions, tracts, depth):
    """tract_qual_experiment.py's genotype comparison, one row per period class and metric."""
    scorer = os.path.join(ctx.project_dir, "benchmarks", "lib", "tract_qual_experiment.py")
    row_base = dict(base, caller=arms.ARM_LABELS[callset.arm],
                    caller_version=ctx.version(arms.ARM_TOOL[callset.arm]))
    kind = arms.ARM_KIND[callset.arm]

    def action(d):
        hostjobs.skip_with(ctx, callset)
        genotype = os.path.join(d, "genotype.tsv")
        cmd = (f"python3 {q(scorer)} --reference {q(arms.GRCH38)} --truth {q(truth)} "
               f"--query {q(callset.path)} --confident-bed {q(regions)} --tract-bed {q(tracts)} "
               f"--arm {callset.arm} --ground hg002_tr --depth {depth} --sample HG002 "
               f"--calibration-out {q(os.path.join(d, 'calibration.tsv'))} "
               f"--sweep-out {q(os.path.join(d, 'sweep.tsv'))} "
               f"--genotype-out {q(genotype)} --genotype-sample HG002")
        run_commands([("tract_score", cmd)], os.path.join(d, "logs"))
        rows = []
        with open(genotype) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                rec = dict(zip(header, line.rstrip("\n").split("\t")))
                cls = rec["period_class"]
                for metric in GENOTYPE_METRICS:
                    value = rec.get(metric, "")
                    if value in (".", ""):
                        continue
                    rows.append(row(row_base, variant_class=f"STR {cls}", stage="tract genotype",
                                    metric=metric, value=value))
        return {"rows": rows}

    return Step(f"{callset.step_name}/score_tracts", "score", action,
                lambda: {"calls": file_identity(callset.path), "truth": file_identity(truth),
                         "tracts": file_identity(tracts), "scorer": file_identity(scorer)},
                attrs={"caller": arms.ARM_CALLER[callset.arm], "depth": depth,
                       "home": "pvc" if kind == "pvc" else "competitor"})
