"""SNP and indel scoring with RTG vcfeval, the same way for every caller (plan §5.1).

vcfeval compares haplotypes rather than records, so a variant written differently by a caller and
by the truth (an indel and a nearby SNP as one record, say) still matches. Each call set is scored
twice:

- **genotype match**, vcfeval's default: the alleles and the genotype must both agree;
- **allele match** (`--squash-ploidy`): the right variant, whatever the genotype.

The gap between the two is the genotype errors. Counts are split into SNPs and everything else
("indel") from vcfeval's own per-class tables, which also give each caller's best F1 over a sweep
of its QUAL threshold.
"""

import gzip
import os
import shlex

from .measure import run_commands
from .results import row

q = shlex.quote


def sdf_command(ctx, reference, out):
    return ("rtg_format", f"java -Xmx4g -jar {q(ctx.rtg_jar)} format -o {q(out)} {q(reference)}")


def prepare_calls(calls, out, sample=None, min_qual=None):
    """The records vcfeval should see: one sample, passing the caller's filters, carrying a variant.

    `min_qual` is the extra threshold freebayes needs (plan D2); other callers get none.
    """
    sample_arg = f"-s {q(sample)}" if sample else ""
    qual = f"-i 'QUAL>={min_qual}'" if min_qual is not None else ""
    return ("prepare_calls",
            f"bcftools view --no-version {sample_arg} -f PASS,. {qual} {q(calls)} "
            f"| bcftools view --no-version -c 1 -Oz -o {q(out)} && bcftools index -t {q(out)}")


def vcfeval_command(ctx, truth, calls, bed, sdf, out_dir, squash, label, truth_sample=None,
                    call_sample=None):
    sample = ""
    if truth_sample or call_sample:
        sample = f"--sample {q((truth_sample or call_sample) + ',' + (call_sample or truth_sample))}"
    return (label,
            f"java -Xmx8g -jar {q(ctx.rtg_jar)} vcfeval -b {q(truth)} -c {q(calls)} -t {q(sdf)} "
            f"-e {q(bed)} -o {q(out_dir)} --ref-overlap --vcf-score-field QUAL "
            f"-T {ctx.threads} {sample} {'--squash-ploidy' if squash else ''}")


def read_roc(path):
    """vcfeval's per-class table: one row per QUAL threshold, from the highest to the lowest."""
    total_baseline = None
    rows = []
    with gzip.open(path, "rt") as fh:
        header = None
        for line in fh:
            if line.startswith("#total baseline variants:"):
                total_baseline = int(line.split(":")[1])
            elif line.startswith("#score"):
                header = line[1:].split()
            elif not line.startswith("#") and line.strip():
                rows.append(dict(zip(header, (float(x) for x in line.split()))))
    return total_baseline, rows


def class_metrics(roc_path):
    """As-written counts (every call, the table's last row) and the best F1 over thresholds."""
    total, roc = read_roc(roc_path)
    if not roc:
        return {"truth": total or 0, "truth_found": 0, "missed": total or 0, "calls_matching": 0,
                "false_calls": 0, "recall": 0.0, "precision": 0.0, "f1": 0.0}, None
    last = roc[-1]
    out = {
        "truth": total,
        "truth_found": int(last["true_positives_baseline"]),
        "missed": int(total - last["true_positives_baseline"]),
        "calls_matching": int(last["true_positives_call"]),
        "false_calls": int(last["false_positives"]),
        "recall": round(last["sensitivity"], 5),
        "precision": round(last["precision"], 5),
        "f1": round(last["f_measure"], 5),
    }
    best = max(roc, key=lambda r: r["f_measure"])
    best_out = {"best_f1": round(best["f_measure"], 5), "best_f1_qual_threshold": best["score"],
                "best_f1_recall": round(best["sensitivity"], 5),
                "best_f1_precision": round(best["precision"], 5)}
    return out, best_out


def score(ctx, step_dir, base, truth, calls, bed, sdf, min_qual=None, call_sample=None,
          truth_sample=None):
    """Score one call set both ways; return (measurement, rows)."""
    prepared = os.path.join(step_dir, "calls.prepared.vcf.gz")
    cmds = [prepare_calls(calls, prepared, sample=call_sample, min_qual=min_qual)]
    m1 = run_commands(cmds, os.path.join(step_dir, "logs"))
    runs = []
    for squash, stage in ((False, "genotype match"), (True, "allele match")):
        out = os.path.join(step_dir, "vcfeval_" + ("allele" if squash else "genotype"))
        runs.append((out, stage, vcfeval_command(
            ctx, truth, prepared, bed, sdf, out, squash, os.path.basename(out),
            truth_sample=truth_sample, call_sample=call_sample)))
    m2 = run_commands([r[2] for r in runs], os.path.join(step_dir, "logs"), max_concurrent=1)
    rows = []
    for out, stage, _ in runs:
        for cls, roc in (("SNP", "snp_roc.tsv.gz"), ("indel", "non_snp_roc.tsv.gz")):
            metrics, best = class_metrics(os.path.join(out, roc))
            for k, v in metrics.items():
                rows.append(row(base, variant_class=cls, stage=stage, metric=k, value=v))
            if best:
                for k, v in best.items():
                    rows.append(row(base, variant_class=cls, stage=stage, metric=k, value=v))
    return m2, rows
