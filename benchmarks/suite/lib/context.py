"""Settings shared by every benchmark in one run, and where each step's directory is."""

import os
import subprocess
import sys

SUITE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJECT_DIR = os.path.dirname(os.path.dirname(SUITE_DIR))


class Context:
    def __init__(self, args):
        self.project_dir = PROJECT_DIR
        self.slice = args.slice
        base = os.path.join(PROJECT_DIR, "benchmark_results")
        self.results_root = os.path.join(base, "slice") if args.slice else base
        self.data_root = os.path.abspath(args.data_root)
        self.threads = args.threads
        self.commit = args.pvc_commit
        self.machine = args.machine
        self.rerun_competitors = args.rerun_competitors
        self.rescore = args.rescore
        self.redo_steps = set(args.redo_step or [])
        self.only = _split(args.only)
        self.depths = _split(args.depth)
        self.callers = _split(args.caller)
        self.external_config = getattr(args, "external_config", None)
        self.tomato_crams = getattr(args, "tomato_crams", None)
        self.tomato_reference = getattr(args, "tomato_reference", None)
        self.tomato_regions = getattr(args, "tomato_regions", None)
        self.giab_alignments = getattr(args, "giab_alignments", None)
        sizes = getattr(args, "tomato_sizes", None)
        self.tomato_sizes = [int(x) for x in sizes.split(",")] if sizes else None
        self.wait_for_quiet = not args.no_wait
        self.quiet_cores = args.quiet_cores
        self.quiet_minutes = args.quiet_minutes
        self.stop_at = args.stop_at_epoch
        self.redo_busy = args.redo_busy
        # The fastest calibration among finished steps, the yardstick of --redo-busy.
        self.fastest_calibration = _fastest_calibration(self.results_root)
        self.pvc_bin = os.path.join(PROJECT_DIR, "target-container", "release", "pvc")
        self.gatk = os.environ.get("GATK_BIN", "/opt/gatk/gatk")
        self.rtg_jar = os.path.join(PROJECT_DIR, "tmp", "tools", "rtg-tools-3.13", "RTG.jar")
        self.java_heap = "4g"
        os.makedirs(self.results_root, exist_ok=True)
        self._log = open(os.path.join(self.results_root, "suite.log"), "a")
        self._versions = {}

    # ---- where things go ------------------------------------------------------------------

    def step_dir(self, step):
        home = step.kind
        if home == "score":
            home = step.attrs.get("home", "pvc")
        if home == "input":
            return os.path.join(self.results_root, "inputs", step.name)
        if home == "competitor":
            return os.path.join(self.results_root, "competitors", step.name)
        return os.path.join(self.pvc_root(), step.name)

    def pvc_root(self):
        return os.path.join(self.results_root, "pvc", self.commit)

    def data(self, *parts):
        """A path inside the data checkout (the one holding benchmarks/ and its data)."""
        return os.path.join(self.data_root, *parts)

    # ---- selection ------------------------------------------------------------------------

    def selected(self, step):
        attrs = step.attrs
        if self.depths and "depth" in attrs and step.kind != "input":
            if str(attrs["depth"]) not in self.depths:
                return False
        if self.callers and "caller" in attrs and step.kind in ("competitor", "pvc", "score"):
            if attrs["caller"] not in self.callers:
                return False
        return True

    # ---- tool versions --------------------------------------------------------------------

    def version(self, tool):
        if tool not in self._versions:
            self._versions[tool] = _read_version(self, tool)
        return self._versions[tool]

    def log(self, message):
        print(message, flush=True)
        self._log.write(message + "\n")
        self._log.flush()


def _fastest_calibration(root):
    import glob
    import json

    best = None
    for path in glob.glob(os.path.join(root, "**", ".done"), recursive=True):
        try:
            with open(path) as fh:
                cal = json.load(fh).get("calibration_s")
        except (OSError, ValueError):
            continue
        if cal and (best is None or cal < best):
            best = cal
    return best


def _split(value):
    if not value:
        return set()
    return {v.strip() for v in value.split(",") if v.strip()}


def _first_line(cmd):
    try:
        out = subprocess.run(cmd, shell=True, capture_output=True, text=True, errors="replace",
                             timeout=120)
    except subprocess.TimeoutExpired:
        return "unknown"
    text = (out.stdout + out.stderr).strip().splitlines()
    return text[0].strip() if text else "unknown"


def _read_version(ctx, tool):
    if tool == "pvc":
        return f"pvc {ctx.commit}"
    if tool == "gatk":
        out = subprocess.run([ctx.gatk, "--version"], capture_output=True, text=True)
        for line in (out.stdout + out.stderr).splitlines():
            if "The Genome Analysis Toolkit (GATK)" in line:
                return "GATK " + line.split()[-1]
        return "GATK unknown"
    if tool == "freebayes":
        return "freebayes " + _first_line("freebayes --version").split()[-1]
    if tool == "hipstr":
        return "HipSTR " + _first_line("HipSTR --version 2>&1").split()[-1]
    if tool == "bcftools":
        return _first_line("bcftools --version")
    if tool == "samtools":
        return _first_line("samtools --version")
    if tool == "minimap2":
        return "minimap2 " + _first_line("minimap2 --version")
    if tool == "deepvariant":
        from .hostjobs import DEEPVARIANT_VERSION
        return f"DeepVariant {DEEPVARIANT_VERSION}"
    if tool == "vcfeval":
        return "RTG Tools " + _first_line(f"java -jar {ctx.rtg_jar} version 2>&1 | grep -i 'product'").split()[-1]
    sys.exit(f"no version reader for {tool}")
