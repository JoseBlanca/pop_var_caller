"""A frame test: steps that take a few seconds, used to check stopping and resuming.

Not part of a full run. `--only smoke` runs it.
"""

import os

from lib.measure import run_commands
from lib.results import row, time_rows
from lib.steps import Step

TITLE = "Frame test"
DESCRIPTION = "Not a benchmark: three short steps that exercise the step runner."
IN_FULL_RUN = False


def steps(ctx):
    out = []
    for i in range(3):
        def action(step_dir, i=i):
            m = run_commands([(f"sleep{i}", f"sleep 4; echo {i} > {step_dir}/out.txt")],
                             os.path.join(step_dir, "logs"))
            base = dict(benchmark="smoke", dataset="none", samples=1, depth="-", region_bp="-",
                        variant_class="-", caller="pvc" if i % 2 else "gatk",
                        caller_version="test")
            return {"measurement": m, "rows": time_rows(base, m, f"sleep{i}")}
        kind = "pvc" if i % 2 else "competitor"
        out.append(Step(f"smoke/step{i}", kind, action, {"i": i}, attrs={"caller": "x"}))
    return out
