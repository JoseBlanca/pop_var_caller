"""SNP and indel scoring with RTG vcfeval, the same way for every caller (plan §5.1).

vcfeval compares haplotypes rather than records, so a variant written differently by a caller and
by the truth (an indel and a nearby SNP as one record, say) still matches. Each call set is scored
twice:

- **genotype match**, vcfeval's default: the alleles and the genotype must both agree;
- **allele match** (`--squash-ploidy`): the right variant, whatever the genotype.

The gap between the two is the genotype errors. Counts are split into SNPs and everything else
("indel") from vcfeval's own per-class tables, which also give each caller's best F1 over a sweep
of its QUAL threshold.

**A missed variant is split into "called with the wrong genotype" and "not called".** The two
counts are not the difference of the two scorings' `missed`: the two match variants differently,
and allele scoring can miss a variant genotype scoring found (pvc's indels on chromosome 20 at full
depth: 210 against 177). Instead each truth record genotype scoring missed is looked up among the
truth records allele scoring found: found there, its allele was called and its genotype was wrong;
not found, it was not called. "Not called" is taken as `missed` less the wrong genotypes, so the
two always add up to vcfeval's own count; the class of a record is judged here by its alleles'
lengths, which can differ from vcfeval's on a rare record.
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
                # A caller with no QUAL gives one row whose score is "None".
                rows.append(dict(zip(header, (_number(x) for x in line.split()))))
    return total_baseline, rows


def _number(text):
    return None if text == "None" else float(text)


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
    scored = [r for r in roc if r["score"] is not None]
    if not scored:
        return out, None
    best = max(scored, key=lambda r: r["f_measure"])
    best_out = {"best_f1": round(best["f_measure"], 5), "best_f1_qual_threshold": best["score"],
                "best_f1_recall": round(best["sensitivity"], 5),
                "best_f1_precision": round(best["precision"], 5)}
    return out, best_out


def _truth_records(path):
    """class -> the truth records in one of vcfeval's per-record files, by position and alleles."""
    out = {"SNP": set(), "indel": set()}
    if not os.path.exists(path):
        return out
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt = line.split("\t", 5)[:5]
            snp = len(ref) == 1 and all(len(a) == 1 for a in alt.split(","))
            out["SNP" if snp else "indel"].add((chrom, pos, ref, alt))
    return out


def wrong_genotypes(step_dir):
    """class -> truth records genotype scoring missed whose allele allele scoring found."""
    missed = _truth_records(os.path.join(step_dir, "vcfeval_genotype", "fn.vcf.gz"))
    found = _truth_records(os.path.join(step_dir, "vcfeval_allele", "tp-baseline.vcf.gz"))
    return {cls: len(missed[cls] & found[cls]) for cls in missed}


def miss_split_rows(step_dir, base, missed_by_class):
    """The `missed` of genotype scoring split into wrong genotype and not called (module doc)."""
    rows = []
    for cls, wrong in wrong_genotypes(step_dir).items():
        if cls not in missed_by_class:
            continue
        wrong = min(wrong, missed_by_class[cls])
        rows.append(row(base, variant_class=cls, stage="genotype match",
                        metric="missed_wrong_genotype", value=wrong))
        rows.append(row(base, variant_class=cls, stage="genotype match",
                        metric="missed_not_called", value=missed_by_class[cls] - wrong))
    return rows


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
    genotype_missed = {}
    for out, stage, _ in runs:
        for cls, roc in (("SNP", "snp_roc.tsv.gz"), ("indel", "non_snp_roc.tsv.gz")):
            metrics, best = class_metrics(os.path.join(out, roc))
            if stage == "genotype match":
                genotype_missed[cls] = metrics["missed"]
            for k, v in metrics.items():
                rows.append(row(base, variant_class=cls, stage=stage, metric=k, value=v))
            if best:
                for k, v in best.items():
                    rows.append(row(base, variant_class=cls, stage=stage, metric=k, value=v))
    rows.extend(miss_split_rows(step_dir, base, genotype_missed))
    return m2, rows
