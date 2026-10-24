"""The second scorer: pvc_eval, which sorts every error into one of three types.

vcfeval (lib/vcfeval.py) counts false calls and missed variants. It cannot count genotypes a
caller leaves missing (`./.`): it marks such a call ignored and never pairs it with the truth, so
a `./.` at a true variant becomes a missed variant and one where the truth has none vanishes
(checked on a three-site test, 24 October 2026). pvc_eval (github.com/JoseBlanca/pvc_eval) gives
every sample at every compared locus one outcome, the truth being its file A and the calls its file
B, and those outcomes fall into the three error types the figures draw (the owner's names):

| error type | pvc_eval outcomes |
|---|---|
| wrong GT (truth ref) | only in B |
| wrong GT (truth var) | only in A; same alleles, other genotype; different |
| no call (./.) | missing in B |

and the rest are right (same genotype; reference in both) or not scored (missing in A or in both,
too complex, invalid, partly outside the regions). Every outcome's count is kept as it is too, in
rows of their own, so nothing pvc_eval says is lost to the grouping.

pvc_eval counts loci, where vcfeval counts records, and it puts the loci at STRs in a class of
their own ("STR") beside SNP, indel and complex; so its counts sit close to vcfeval's but are not
the same, and vcfeval stays as the check (on the GIAB panel, HG002 at 30x, pvc fitted: 23 wrong at
true variants against vcfeval's 26 missed; 3 against 6 false calls).

It is given the calls before vcfeval's filter for records carrying a variant, so the `./.` and `0/0`
genotypes reach it; the caller's filters and freebayes' QUAL floor still apply. STR regions come
from `pvc_eval str-regions`, run once per reference (lib/arms.py), not from a scan at every step.
pvc_eval runs on the dev image's free-threaded Python 3.15 (Containerfile), from the checkout
run_suite.sh found and mounted.
"""

import csv
import os
import re
import shlex
import subprocess

from .measure import run_commands
from .results import row

q = shlex.quote

# pvc_eval's outcomes (its doc/spec.md §5.5), grouped into the three error types and the rest.
ERROR_TYPE = {
    "only in B": "wrong_gt_truth_ref",
    "only in A": "wrong_gt_truth_var",
    "same alleles, other genotype": "wrong_gt_truth_var",
    "different": "wrong_gt_truth_var",
    "missing in B": "no_call",
    "same genotype": "right",
    "reference in both": "right",
}
GROUPS = ["wrong_gt_truth_ref", "wrong_gt_truth_var", "no_call", "right", "not_scored"]
KINDS = ["SNP", "indel", "STR", "complex"]


def available(ctx):
    return bool(ctx.pvc_eval_dir)


def _metadata_dir(ctx):
    """A directory holding the package metadata of the checkout, for PYTHONPATH.

    pvc_eval reads its own version from its installed package's metadata (`str-regions` writes it
    in its BED's header), and a checkout run from PYTHONPATH has none, so the call fails
    ("No package metadata was found for pvc-eval"). Installing it would need the network at every
    run, to fetch its build backend; a minimal `.dist-info` with its name and the version in its
    pyproject.toml is all the lookup reads. Written once per pvc_eval commit, under the results."""
    with open(os.path.join(ctx.pvc_eval_dir, "pyproject.toml")) as fh:
        found = re.search(r'^version\s*=\s*"([^"]+)"', fh.read(), re.MULTILINE)
    version = found.group(1) if found else "0"
    path = os.path.join(ctx.results_root, "inputs", "pvc_eval_metadata", ctx.pvc_eval_commit)
    info = os.path.join(path, f"pvc_eval-{version}.dist-info")
    if not os.path.exists(os.path.join(info, "METADATA")):
        os.makedirs(info, exist_ok=True)
        with open(os.path.join(info, "METADATA"), "w") as fh:
            fh.write(f"Metadata-Version: 2.1\nName: pvc-eval\nVersion: {version}\n")
    return path


def _command(ctx, args):
    """A pvc_eval command line, run from the mounted checkout on the image's Python 3.15t."""
    path = f"{os.path.join(ctx.pvc_eval_dir, 'src')}:{_metadata_dir(ctx)}"
    return (f"PYTHONPATH={q(path)} PYTHONDONTWRITEBYTECODE=1 {q(ctx.pvc_eval_python)} -c "
            f"{q('import sys; from pvc_eval.cli import main; main(sys.argv[1:])')} {args}")


def str_regions_command(ctx, reference, out):
    return ("pvc_eval_str_regions",
            _command(ctx, f"str-regions --reference {q(reference)} --output {q(out)} "
                          f"--threads {ctx.threads}"))


def _first_sample(path):
    out = subprocess.run(["bcftools", "query", "-l", path], capture_output=True, text=True,
                         check=True)
    names = out.stdout.split()
    if not names:
        raise ValueError(f"{path} has no sample")
    return names[0]


def in_reference_order(reference, vcf, out):
    """A shell command writing `vcf` to `out` in the reference's contig order, which pvc_eval
    requires and vcfeval does not: HipSTR writes its contigs by name (chr19 before chr2), and two
    of GIAB's panel truth files have chr7 before chr6. The header's contigs are replaced by the
    `.fai`'s, then the records sorted by them; a file already in order comes out the same."""
    reheaded = out + ".reheaded.vcf.gz"
    return (f"bcftools reheader --fai {q(reference + '.fai')} -o {q(reheaded)} {q(vcf)} && "
            f"bcftools sort -Oz -o {q(out)} {q(reheaded)} && rm {q(reheaded)} && "
            f"bcftools index -t {q(out)}")


def score(ctx, step_dir, base, truth, calls, bed, reference, str_regions, min_qual=None,
          call_sample=None, truth_sample=None):
    """Compare one call set with the truth; return (measurement, rows)."""
    out = os.path.join(step_dir, "pvc_eval")
    os.makedirs(out, exist_ok=True)
    filtered = os.path.join(out, "calls.filtered.vcf.gz")
    kept = os.path.join(out, "calls.kept.vcf.gz")
    sorted_truth = os.path.join(out, "truth.sorted.vcf.gz")
    sample_arg = f"-s {q(call_sample)}" if call_sample else ""
    qual = f"-i 'QUAL>={min_qual}'" if min_qual is not None else ""
    prepare = ("pvc_eval_prepare",
               f"bcftools view --no-version {sample_arg} -f PASS,. {qual} -Oz -o {q(filtered)} "
               f"{q(calls)} && {in_reference_order(reference, filtered, kept)} && "
               f"rm {q(filtered)} && {in_reference_order(reference, truth, sorted_truth)}")
    run_commands([prepare], os.path.join(step_dir, "logs"))
    truth = sorted_truth
    sample_map = os.path.join(out, "samples.tsv")
    with open(sample_map, "w") as fh:
        fh.write(f"{truth_sample or _first_sample(truth)}\t{call_sample or _first_sample(kept)}\n")
    summary = os.path.join(out, "summary.tsv")
    compare = ("pvc_eval_compare", _command(
        ctx, f"compare --reference {q(reference)} --a {q(truth)} --b {q(kept)} "
             f"--regions {q(bed)} --str-regions {q(str_regions)} --sample-map {q(sample_map)} "
             f"--output {q(os.path.join(out, 'loci.tsv'))} --summary {q(summary)}"))
    m = run_commands([compare], os.path.join(step_dir, "logs"))
    return m, summary_rows(summary, base)


def summary_rows(path, base):
    """pvc_eval's summary as rows: each outcome's count per kind (stage "pvc_eval outcome"), and
    the counts grouped into the error types (stage "error type"), summed over strata."""
    raw = {}
    with open(path) as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            for kind in KINDS:
                key = (rec["outcome"], kind)
                raw[key] = raw.get(key, 0) + int(rec[kind])
    rows = []
    grouped = {(g, k): 0 for g in GROUPS for k in KINDS}
    for (outcome, kind), n in sorted(raw.items()):
        rows.append(row(base, variant_class=kind, stage="pvc_eval outcome", metric=outcome,
                        value=n))
        grouped[(ERROR_TYPE.get(outcome, "not_scored"), kind)] += n
    for (group, kind), n in grouped.items():
        rows.append(row(base, variant_class=kind, stage="error type", metric=group, value=n))
    return rows
