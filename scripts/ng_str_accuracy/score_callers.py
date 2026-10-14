"""The stuttering benchmark's scorer, over GATK and the callers named in <run>/callers.list.

    python3 scripts/ng_str_accuracy/score_callers.py <run-dir>

`scripts/ng_str_stutter_score.py` with its caller list replaced, so a run directory holding only
GATK's calls and a few builds' (`call_and_score.sh` writes the list) is scored in seconds rather
than over the first pass's nine callers. Rewrites <run-dir>/score_cells.tsv (one row per caller,
segment and sample) and score_summary.tsv, and, when the run has STR bundles,
score_bundle_cells.tsv and score_bundle_summary.tsv. Set SKIP_VCFEVAL=1 to skip vcfeval.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import ng_str_stutter_score as scorer  # noqa: E402

run = sys.argv[1]
names = [line.strip() for line in open(os.path.join(run, "callers.list")) if line.strip()]
scorer.CALLERS = [("gatk", "gatk.vcf.gz")] + [(name, name + ".vcf") for name in names]
scorer.main()
