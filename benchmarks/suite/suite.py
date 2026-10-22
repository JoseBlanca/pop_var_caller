"""The benchmark suite, run inside the dev container. Start it with run_suite.sh, not directly.

run_suite.sh, on the host, records the machine and the host's load, then calls this inside the
container. This builds pvc, runs every selected benchmark's steps in order (skipping finished
ones), and writes the results table and the summary for the current pvc commit.
"""

import argparse
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bench import BENCHMARKS  # noqa: E402
from lib.context import Context  # noqa: E402
from lib.results import build_table, write_summary  # noqa: E402
from lib.steps import Runner, SkipStep, StopRequested, SuiteError  # noqa: E402
from lib.measure import CommandFailed  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--only", help="comma-separated benchmarks: " + ",".join(BENCHMARKS))
    p.add_argument("--depth", help="comma-separated depths, e.g. 10x,30x")
    p.add_argument("--caller", help="comma-separated callers: pvc,gatk,freebayes,hipstr")
    p.add_argument("--threads", type=int, default=6)
    p.add_argument("--slice", action="store_true", help="the small test configuration")
    p.add_argument("--rerun-competitors", action="store_true")
    p.add_argument("--rescore", action="store_true", help="redo every scoring step")
    p.add_argument("--redo-step", action="append", help="redo this step (its name), any kind")
    p.add_argument("--external-config", help="configuration file of an external benchmark")
    p.add_argument("--tomato-crams", help="directory of tomato CRAMs (a machine profile sets it)")
    p.add_argument("--tomato-reference", help="SL4.0 FASTA (a machine profile sets it)")
    p.add_argument("--tomato-sizes", help="comma-separated sample counts for the speed benchmark")
    p.add_argument("--no-wait", action="store_true",
                   help="start timed steps even when the machine is busy (their times are marked)")
    p.add_argument("--quiet-cores", type=float, default=2.0)
    p.add_argument("--quiet-minutes", type=float, default=2.0)
    p.add_argument("--stop-at-epoch", type=float, default=None, help="set by run_suite.sh --stop-at")
    p.add_argument("--redo-busy", action="store_true",
                   help="rerun pvc and competitor steps whose times were taken on a busy machine")
    p.add_argument("--summary-only", action="store_true")
    p.add_argument("--status", action="store_true")
    p.add_argument("--no-build", action="store_true")
    p.add_argument("--data-root", required=True)
    p.add_argument("--pvc-commit", required=True)
    p.add_argument("--machine", required=True)
    return p.parse_args()


def selected_benchmarks(ctx):
    names = list(BENCHMARKS)
    if ctx.only:
        unknown = ctx.only - set(names)
        if unknown:
            sys.exit(f"unknown benchmark(s): {', '.join(sorted(unknown))}; known: {', '.join(names)}")
        names = [n for n in names if n in ctx.only]
    else:
        names = [n for n in names if getattr(BENCHMARKS[n], "IN_FULL_RUN", True)]
    return names


def build_pvc(ctx):
    ctx.log("building pvc (release)")
    # The stuttering benchmark's simulator and typing checks are examples of this crate.
    examples = ["ng_str_stutter_benchmark", "ng_typed_region_dump", "ng_str_bundle_runs"]
    cmd = ["cargo", "build", "--release", "--quiet", "--bin", "pvc"]
    for e in examples:
        cmd += ["--example", e]
    out = subprocess.run(cmd,
                         cwd=ctx.project_dir, capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit("pvc did not build:\n" + out.stderr[-4000:])


def _steps(ctx, name):
    try:
        return BENCHMARKS[name].steps(ctx)
    except SuiteError:
        raise
    except FileNotFoundError as err:
        raise SuiteError(f"{name}: {err}")


def status(ctx, names):
    runner = Runner(ctx)
    for name in names:
        steps = [s for s in BENCHMARKS[name].steps(ctx) if ctx.selected(s)]
        done = [s for s in steps if runner.is_done(s)]
        seconds = sum((runner.load_done(s).get("measurement") or {}).get("wall_s", 0) for s in done)
        print(f"{name:20s} {len(done):4d} of {len(steps):4d} steps done "
              f"({seconds / 3600:.1f} h of measured work so far)")
        for s in steps:
            if not runner.is_done(s):
                print(f"    next: {s.name}")
                break


def main():
    args = parse_args()
    ctx = Context(args)
    names = selected_benchmarks(ctx)
    if args.status:
        status(ctx, names)
        return 0
    code = 0
    if not args.summary_only:
        if not args.no_build:
            build_pvc(ctx)
        runner = Runner(ctx)
        failures = []
        skipped = []
        try:
            # A failed step is logged and the run goes on, so one broken step does not cost the
            # rest of a night's quiet hours; steps that need its output fail in turn and are
            # listed too. Everything finished is kept, and the failed steps rerun next time.
            for name in names:
                ctx.log(f"== {name}")
                try:
                    steps = _steps(ctx, name)
                except SuiteError as err:
                    ctx.log(f"error: {err}")
                    failures.append((name, str(err)))
                    continue
                for step in steps:
                    if not ctx.selected(step):
                        continue
                    try:
                        runner.run(step)
                    except SkipStep as why:
                        ctx.log(f"skipped: {step.name}: {why}")
                        skipped.append((step.name, str(why)))
                    except (SuiteError, CommandFailed, FileNotFoundError) as err:
                        ctx.log(f"error: {step.name}: {err}")
                        failures.append((step.name, str(err)))
            if skipped:
                reasons = sorted({why for _, why in skipped})
                ctx.log(f"{len(skipped)} step(s) skipped, this machine cannot run them: "
                        + "; ".join(reasons))
            if failures:
                ctx.log(f"{len(failures)} step(s) failed; they rerun next time:")
                for name, why in failures:
                    ctx.log(f"    {name}: {why[:300]}")
                code = 2
        except StopRequested as why:
            ctx.log(f"{why}; run the same command to resume")
            code = 3
        except (SuiteError, CommandFailed) as err:
            ctx.log(f"error: {err}")
            code = 2
        except Exception:
            import traceback
            ctx.log("error:\n" + traceback.format_exc())
            code = 2
        except KeyboardInterrupt:
            ctx.log("interrupted; the step that was running will be redone next time")
            code = 130
        ctx.log(f"{runner.done_count} steps run, {runner.skipped_count} already done")
    rows = build_table(ctx)
    latest = os.path.join(ctx.results_root, "latest")
    if os.path.islink(latest):
        os.remove(latest)
    if not os.path.exists(latest):
        os.symlink(os.path.relpath(ctx.pvc_root(), ctx.results_root), latest)
    path = write_summary(ctx, rows, BENCHMARKS)
    ctx.log(f"{len(rows)} rows in {os.path.join(ctx.pvc_root(), 'results.tsv')}; summary {path}")
    return code


if __name__ == "__main__":
    sys.exit(main())
