"""Measure one shell command the suite's way, from inside another tool's image.

    python3 measure_cli.py <result.json> <log-dir> <label> -- <shell command>

Used by the host job runner (host_jobs.sh) to run DeepVariant and GLnexus in their own images while
measuring them exactly as the suite measures every other caller (lib/measure.py: wall and CPU time,
the largest process's peak memory, the summed peak of all processes). Standard library only.
"""

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from measure import CommandFailed, run_commands  # noqa: E402


def main():
    out, log_dir, label = sys.argv[1:4]
    assert sys.argv[4] == "--", "usage: measure_cli.py <result.json> <log-dir> <label> -- <command>"
    command = " ".join(sys.argv[5:])
    try:
        m = run_commands([(label, command)], log_dir)
        result = {"ok": True, "measurement": m}
    except CommandFailed as err:
        result = {"ok": False, "error": str(err)}
    with open(out + ".tmp", "w") as fh:
        json.dump(result, fh)
    os.replace(out + ".tmp", out)
    sys.exit(0 if result["ok"] else 1)


if __name__ == "__main__":
    main()
