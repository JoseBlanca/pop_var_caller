"""Run the stuttering benchmark's scorer over a chosen set of callers.

    python3 stutter_score_callers.py <view-dir> <caller>=<file> ...

`scripts/ng_str_stutter_score.py` with its caller list replaced, as
`scripts/ng_str_accuracy/score_callers.py` does, so the suite can score the competitors once and
pvc once per commit, each in a directory of its own. vcfeval is skipped (SKIP_VCFEVAL=1): this
benchmark is quoted for its STR genotypes, and the human benchmarks carry the SNP and indel figures.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "scripts"))
import ng_str_stutter_score as scorer  # noqa: E402

view = sys.argv[1]
scorer.CALLERS = [tuple(arg.split("=", 1)) for arg in sys.argv[2:]]
os.environ["SKIP_VCFEVAL"] = "1"
sys.argv = [sys.argv[0], view]
scorer.main()
